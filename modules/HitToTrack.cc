/** \class HitToTrack
 *
 * Takes hits produced by the HitFinder class
 * Performs track reconstruction
 * Returns Delphes tracks 
 *
 *  \authors W. Fawcett
 *
 */


#include "modules/HitToTrack.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/LineParameters.h"

// cstdlib 
#include <stdexcept>
#include <iostream>

// testing
#include <ctime>


// sign: god knows why this isn't a standard function!
  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//------------------------------------------------------------------------------

inline float quickDeltaPhi(float phi1, float phi2){
  if(phi1<0) phi1+= 2*M_PI;
  if(phi2<0) phi2+= 2*M_PI;
  float dPhi= fabs(phi1 - phi2);
  if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
  return dPhi; 
}


//------------------------------------------------------------------------------

HitToTrack::HitToTrack() 
{
}

//------------------------------------------------------------------------------

HitToTrack::~HitToTrack()
{
}


//------------------------------------------------------------------------------

void HitToTrack::Init()
{

  m_debug = GetBool("debug", false); 

  if(m_debug) std::cout << "HitToTrack::Init()" << std::endl;

  ////////////////////////////////////////////////////
  // Input parameters (copied from ParticlePropagator) 
  ////////////////////////////////////////////////////
  fInputArray = ImportArray(GetString("InputArray", "Delphes/hit"));
  //fInputArray   = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // window size for the seeding pattern algorithm 
  fPhiSeedingWindowSize = GetDouble("PhiSeedingWindowSize", 0.06); // radians  
  fEtaSeedingWindowSize = GetDouble("EtaSeedingWindowSize", 0.1); // 

  fSeedingTrackPtMinGeV = GetDouble("SeedingTrackPtMinGeV", 2.0); // GeV
  fBField = GetDouble("Bz", 4.0); // solenoid magnetic field [T]

  fLuminousRegionSize = GetDouble("LuminousRegionSize", 53.0); // length of the region covering 1 sigma of the luminous region
  fNVertxSigma = GetInt("NVertexSigma", 4); // number of sigma of luminout region to consider

  fZResiduumTolerance = GetDouble("ZResiduumTolerance", 0.5); // 0.5

  // Output arrays 
  fTrackOutputArray = ExportArray(GetString("OutputArray", "tracks"));

  // Print configuration
  std::cout << "HitToTrack::Init(): Initialized HitToTrack with parameters:" << std::endl; 
  std::cout << "HitToTrack::Init(): PhiSeedingWindowSize: " << fPhiSeedingWindowSize << std::endl;
  std::cout << "HitToTrack::Init(): EtaSeedingWindowSize: " << fEtaSeedingWindowSize << std::endl;
  std::cout << "HitToTrack::Init(): SeedingTrackPtMinGeV : " << fSeedingTrackPtMinGeV  << std::endl;
  std::cout << "HitToTrack::Init(): BField: " << fBField << std::endl;
  std::cout << "HitToTrack::Init(): LuminousRegionSize : " << fLuminousRegionSize  << std::endl;
  std::cout << "HitToTrack::Init(): NVertxSigma : " << fNVertxSigma  << std::endl;
  std::cout << "HitToTrack::Init(): ZResiduumTolerance: " << fZResiduumTolerance << std::endl;

  if(m_debug) std::cout << "Init(): Defined input and output arrays" << std::endl;

}

//------------------------------------------------------------------------------

void HitToTrack::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

std::vector<Candidate*> concatenateVector(std::vector<Candidate*>& A, std::vector<Candidate*>& B){
  // https://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
  std::vector<Candidate*> AB;
  AB.reserve( A.size() + B.size() ); // preallocate memory
  AB.insert( AB.end(), A.begin(), A.end() );
  AB.insert( AB.end(), B.begin(), B.end() );
  return AB;
}

//------------------------------------------------------------------------------

std::vector<Candidate*> concatenateHitsBasedOnLocations(std::map<std::string, std::vector<Candidate*>>& hitMap, std::vector<std::string>& locations){
  // Return a vector of all hits in all of the regions selected by locations 
  // Can probably make this more efficient ? 
  std::vector<Candidate*> newVec;
  for(const auto& location : locations){
    /******************
    try{
      newVec = concatenateVector(newVec, hitMap.at(location));
    }
    catch(const std::out_of_range& oor){
      std::cout << "out f range error with: " << location << std::endl;
      std::cout << oor.what() << std::endl;
    }
    ********************/
      newVec = concatenateVector(newVec, hitMap[location]);
    
  }
  return newVec; 
}



//------------------------------------------------------------------------------

std::vector< std::vector<Candidate*> >  HitToTrack::FindSeedsTriplet(
    std::map<int, std::vector<Candidate*> >& hc,
    std::map<std::string, std::vector<Candidate*> >& hitMap, 
    Location& loc, 
    std::vector<int> layerIDs
    ) const { // previously called associateHitsSimple


    //////////////////////////////////////////
    // Simplest possible algorithm to match seed tracks from a triplet 
    //////////////////////////////////////////

    // Calculate the min and max Z displacement of the track based on the vertex selection constraints
    float maxZ = fNVertxSigma * fLuminousRegionSize;
    float minZ = -1 * maxZ; 


    // get the inner and outer barrel radii
    int innerLayerID = layerIDs.at(0);
    int outerLayerID = layerIDs.back();
    int middleLayerID = 1; 
    double rInner = hc[innerLayerID].at(0)->Position.Perp();
    double rOuter = hc[outerLayerID].at(0)->Position.Perp(); 

    // reserve some space for the seeds
    std::vector<std::vector<Candidate*>> trackSeeds; 
    trackSeeds.reserve( hc[outerLayerID].size() ); 

    // Calculate the phi window in which the hits in the outer layer must match
    double bendingRadius = 1000 * 3.35 * fSeedingTrackPtMinGeV / fBField ; // [mm]
    double phiWindow = fabs( acos(rInner / (2*bendingRadius)) - acos(rOuter / (2*bendingRadius)) );
    phiWindow *= 2; // multiply by two, to have the deviation travelling in either direction. 



    // Draw a line between the hit in the innermost and outermost layer
    // See if there is a hit on the line in the intermediate layer (within some tolerance)
    for(const auto& innerHit : hc[innerLayerID]){

      const float zInner = innerHit->Position.Z(); // Z
      const float phiInner = innerHit->Phi;

      // get locations (areas) for other hits
      std::string innerHitLocation = loc.locationFromHit(innerHit); 
      std::vector<std::string> outerHitLocations  = loc.listOfLocationsInLayer(innerHitLocation, outerLayerID );
      std::vector<std::string> middleHitLocations = loc.listOfLocationsInLayer(innerHitLocation, middleLayerID);
      //std::cout << "There are " << outerHitLocations.size() << " outer hit locations: " << std::endl;
      //for(auto l : outerHitLocations) std::cout << "\t " <<  l << std::endl;

      // get vector of hits defined by list of locations
      std::vector<Candidate*> outerHitVector  = concatenateHitsBasedOnLocations(hitMap, outerHitLocations);
      std::vector<Candidate*> middleHitVector = concatenateHitsBasedOnLocations(hitMap, middleHitLocations);


      //for(const auto& outerHit : hc[outerLayerID]){
      for(const auto& outerHit : outerHitVector){

        // must be within phi criteria  
        if( quickDeltaPhi(phiInner, outerHit->Phi) > phiWindow) continue; 

        // calculate parameters of line from inner hit to outer hit 
        const float zOuter = outerHit->Position.Z(); // Z 
        LineParameters params;
        params.calculateLineParameters(zInner, rInner, zOuter, rOuter);

        // reject if line does not point to within 3 sigma of the luminous region
        float beamlineIntersect = params.x_intercept() ;
        //std::cout << "beamlineIntersect: " << beamlineIntersect << ", min: " << minZ << ", max: " << maxZ << std::endl;
        if(beamlineIntersect > maxZ || beamlineIntersect < minZ) continue;

        // intersection of the line with the intermediate layer
        float intersect = (582.0 - params.y_intercept())/params.gradient();

        //for(const auto& intermediateHit : hc[middleLayerID]){
        for(const auto& intermediateHit : middleHitVector){

          float middleZresiduum = intermediateHit->Position.Z() - intersect;
          // only select if intermediate hit matches within tolerance along Z  
          if(fabs(middleZresiduum) < fZResiduumTolerance){

            // reject the intermediate hit if it is also outside the phi window
            if( quickDeltaPhi(phiInner, intermediateHit->Phi) > phiWindow) continue; 

            // Three hits are matched -> a track 
            std::vector<Candidate*> matchedHits;
            matchedHits.reserve(4); // prevent vector from having to grow 
            matchedHits.push_back(innerHit);
            matchedHits.push_back(intermediateHit);
            matchedHits.push_back(outerHit);

            trackSeeds.push_back(matchedHits); 

          }
        }
      }
    }

    return trackSeeds;
}





TrackParameterSet HitToTrack::CalculateTrackParametersTriplet(std::vector<Candidate*>& seeds) const{

  /****************
   * First calculate the track parameters assuming a beamline constraint (and using the first and third hit)
   * Then calculate the track parameters using only the three hits
   * Good tracks will have consistent parameters
   * **************/

  ///////////////////////////////
  // Calculate track parameters with beamline constraint 
  // Assumes the track originates from (0, 0, z0)
  // Calculates the track parameters using the origin and two other points, corresponding to the innermost and outermost hit
  
  //if(m_debug) std::cout << "HitToTrack::CalculateTrackParametersTriplet()" << std::endl;

  if(seeds.size() != 3){
    std::cerr << "ERROR: more than three hits associate to this track. Algorith is not compatible." << std::endl;
    exit(1);
  }

  Candidate* hit1 = seeds.at(0);
  Candidate* hit2 = seeds.at(1);
  Candidate* hit3 = seeds.at(2);

  // Check hits are properly orderd radially. Shouldn't really happen that they are not 
  if(hit1->Position.Perp() > hit2->Position.Perp() || hit2->Position.Perp() > hit3->Position.Perp()){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    exit(1);
  }

  // copy of hit coordinates
  float x1 = hit1->Position.X();
  float y1 = hit1->Position.Y(); 
  float z1 = hit1->Position.Z(); 

  float x2 = hit2->Position.X();
  float y2 = hit2->Position.Y(); 
  float z2 = hit2->Position.Z(); 

  float x3 = hit3->Position.X();
  float y3 = hit3->Position.Y(); 
  float z3 = hit3->Position.Z(); 

  float r01 = hit1->Position.Perp(); // = sqrt(x^2 + y^2) 
  float r02 = hit2->Position.Perp();
  float r03 = hit3->Position.Perp(); 

  // Calculate the parameters in the longitudinal plane
  // Follows calculations by A. Schoning

  float r13 = hypotf( fabs(x3-x1) , fabs(y3-y1) ); // not sure if fabs actuall needed ...  
  float chord13 = x1*y3 - y1*x3; 
  float PHI1 = 2 * asin( chord13 / (r13*r03) );
  float PHI3 = 2 * asin( chord13 / (r13*r01) );

  // Calculate the track parameters in the transverse plane
  // Given the two points and the beamline constraint
  // This can be solved analytically ... 
   
  // calculate the centers of the circle touching coordinates (0,0) (x1, y1) (x3, y3)
  // described by (x-a)^2 + (y-b)^2 = R^2
  float a = (y3*r01*r01 - y1*r03*r03) / (2*(y3*x1 - y1*x3));
  float b = (x3*r01*r01 - x1*r03*r03) / (2*(x3*y1 - x1*y3));

  // radius of trajectory  
  float radius_OLD = sqrt(a*a + b*b); 
  float radius = (r01 * r03 * r13) / (2*chord13); 

  // kappa (1/radius) 
  float kappa_013 = 1/radius; 

  // transverse momentum
  float pT = 1.199 * fabs(radius/1000); // for a 4T magnetic field and radius in [m], divide by 1000 as length units in [mm]

  // arc lengths 
  double s1 = radius * PHI1;
  double s3 = radius * PHI3; 


  // phi angle given by line tangent to the circle at (0,0) 
  // atan2(y, x)
  float phi = atan2f(-b, a); // will return [-pi, pi] 

  // Assign the track parameters
  float z0 = z1 - s1* ( z3 - z1 ) / (s3 - s1);
  double theta = atan2( (s3-s1), (z3-z1) );  // atan2f returns theta in [-pi, pi] 


  // check eta calculation is correct
  double eta = -1*log( tan( fabs(theta)/2.0 )); // take fabs(theta), want -pi and pi to be treated the same
  if(isnan(eta)){
    std::cerr << "trackParametersBeamlineConstraint(): ERROR: Eta calculation performed incorrectly." << std::endl; 
    std::cerr << "theta: " << theta << std::endl;
    std::cerr << "tan(fabs(theta)/2.0): " << tan(fabs(theta)/2.0) << std::endl;
    std::cerr << "log(tan(fabs(theta))): " << log( tan(fabs(theta)/2.0) ) << std::endl;
    std::cerr << "Curious and curiouser ... " << std::endl;
    eta = -100.0;
    //return false;
  }




  /////////////////////////////////////////////
  // Calculate the track parameters without the beamline constraint
  /////////////////////////////////////////////
  
  // circle described by (x-a)^2 + (y-b)^2 = r^2 
  float b_nbc = ( (x3 - x1)*(r02*r02 - r01*r01) - (x2 - x1)*(r03*r03 - r01*r01) ) / ( 2*( (y2-y1)*(x3-x1) - (y3-y1)*(x2-x1)  ) );
  float a_nbc = ( r02*r02 - r01*r01 - 2*b*(y2 - y1) ) / ( 2*(x2 - x1));
  //float radius_nbc = hypotf( (x1-a_nbc), (y1-b_nbc) );

  // Andre's calculation
  float chord_123 = x2*y3 + x1*y2 + x3*y1 - x3*y2 -x2*y1 - x1*y3;
  float r12 = hypotf( (x2-x1), (y2-y1) );
  float r23 = hypotf( (x3-x2), (y3-y2) );
  float kappa_123 = 2*chord_123 / (r12 * r13 * r23 ); 
  float radius_nbc = fabs(1/kappa_123); 


  // Write the track parameter struct
  TrackParameterSet trackParameters;
  trackParameters.pT = pT;
  trackParameters.d0 = this->CalculateD0(a_nbc, b_nbc, radius_nbc);
  trackParameters.z0 = z0;
  trackParameters.eta = eta; 
  trackParameters.theta = theta;
  trackParameters.phi = phi;
  trackParameters.kappa_013 = kappa_013;
  trackParameters.kappa_123 = kappa_123;
  trackParameters.isFake = isFake(seeds);
  trackParameters.Charge = sgn(radius);  // use Andre's radius 


  return trackParameters; 

}


float HitToTrack::CalculateD0(float x0, float y0, float radius) const{
  // D0 calculated with respect to the detector coordinate origin
  // with track circle parameters (x, y, r)
  
  float distanceToOrigin = hypotf(x0, y0);

  // Sign convention
  // if origin is inside the circle (distanceToOrigin < radius)
  // d0 will be -ve
  //
  // if origin is outside the circle (distanceToOrigin > radius)
  // d0 will be +vq
  return distanceToOrigin - radius;
}

bool HitToTrack::isFake(std::vector<Candidate*>& seeds) const{
  // if any of the hits do not originate from the same particle, then the track is fake
  for(int i=0; i < seeds.size()-1; ++i){
    auto p1 = seeds.at(i)->GetCandidates()->At(0);
    auto p2 = seeds.at(i+1)->GetCandidates()->At(0);
    if( p1->GetUniqueID() != p2->GetUniqueID() ) return true; 
  }
  return false;
}


void HitToTrack::Process()
{

  if(m_debug) std::cout << "HitToTrack::Process()" << std::endl;


  // Define a location object, that defines the size of the eta and phi windows used for seeding 
  //Location loc(0.06, 0.1);
  Location loc(fPhiSeedingWindowSize, fEtaSeedingWindowSize);

  // Two maps for hit groupings 
  std::map<std::string, std::vector<Candidate*>> hitMap;
  std::map<int, std::vector<Candidate*> > hitContainer; 

  // Loop over the hits in the event
  // Group hits into layer-eta-phi regions (for faster matching)
  int nHitsEvent(0);
  int nHitsOuterPt2(0);
  fItInputArray->Reset();
  Candidate* hit; 
  while((hit = static_cast<Candidate*>(fItInputArray->Next())))
  {

    // Map for locations to hits 
    std::string locationString = loc.locationFromHit(hit);
    hitMap[ locationString ].push_back(hit); 

    // Map for layers to hits 
    hitContainer[hit->SurfaceID].push_back(hit); 

    if(hit->SurfaceID == 2){
      if(hit->PT > 2) nHitsOuterPt2++;
    }

    nHitsEvent++;
  }
  if(m_debug) std::cout << "HitToTrack::Process(): event has " << nHitsEvent << " hits in total, and " << nHitsOuterPt2 << " in the outermost layer with pT>2 GeV" << std::endl;


  // Only continue if there is at least 1 hit in the event
  if(nHitsEvent > 0){

    // get a list of the unique layer IDs 
    std::vector<int> layerIDs;
    for(const auto& id : hitContainer){
      layerIDs.push_back(id.first);
    }
    std::sort(layerIDs.begin(), layerIDs.end()); // sort into ascending order

    // Get the seeds 
    std::vector< std::vector<Candidate*> > seedSetVector = this->FindSeedsTriplet(hitContainer, hitMap, loc, layerIDs); 
    
    if(m_debug) std::cout << "HitToTrack::Process(): event has " << seedSetVector.size() << " sets of seeds" << std::endl;

    // Reconstruct the seeds into tracks, apply tighter constraints on the track selection
    for(auto& seeds : seedSetVector){

      // Calculate track parameters 
      TrackParameterSet parameters = this->CalculateTrackParametersTriplet(seeds);

      // Create a new (track) candidate using the outermost hit, 
      // Use the outermost hit so the Xd and L variables are set correctly 
      Candidate* track = static_cast<Candidate*>(seeds.back()->Clone());

      // Assign track parameters 
      track->D0 = parameters.d0;  
      track->DZ = parameters.z0;
      track->PT = parameters.pT;
      track->Phi = parameters.phi;
      track->Eta = parameters.eta; 
      track->Charge = parameters.Charge;
      track->kappa_123 = parameters.kappa_123;
      track->kappa_013 = parameters.kappa_013; 
      track->IsFake = parameters.isFake;

      // Set track momentum TLV
      TLorentzVector momentum;
      float pion_mass = 139.570 / 1000; // GeV 
      momentum.SetPtEtaPhiM(parameters.pT, parameters.eta, parameters.phi, pion_mass);
      track->P = momentum.P();
      track->Momentum = momentum; 
      track->CtgTheta = 1.0/tan( parameters.theta ); 

      // Set track position TLV
      track->Position = momentum; // might not make sense to set position = momemtum, but track position doesn't really mean anything. Would want the correct Eta and Phi coordinates if this was called, though


      // Add the other seeds to the track
      for(auto& seed : seeds){
        track->AddCandidate(seed);
      }
      fTrackOutputArray->Add(track); 

    } // end of loop over seeds
  } // end of if(hasHits)
} // end of Process()

//------------------------------------------------------------------------------
