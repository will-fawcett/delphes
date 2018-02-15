/** \class HitFinder
 *
 *  Uses the particle propagator to extract the location of a particle when intersecting with a specified surface (i.e. a hit)
 *  Two types of surface can be defined, barrels (which all must have the same length), and endcaps (which all must have the same radius).
 *  If a particle does not have sufficient pT to reach the surface, then the hit is not stored.
 *
 *  \authors W. Fawcett
 *
 */


#include "modules/HitFinder.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "classes/Barrel.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include "TVector3.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

// testing
#include <ctime>

using namespace std;




//------------------------------------------------------------------------------

HitFinder::HitFinder() :
  fBarrelLength(0), fEndCapRadius(0)
{
}

//------------------------------------------------------------------------------

HitFinder::~HitFinder()
{
}


//------------------------------------------------------------------------------

void HitFinder::Init()
{

  m_debug=false;

  if(m_debug) std::cout << "HitFinder::Init()" << std::endl;

  ////////////////////////////////////////////////////
  // Input parameters (copied from ParticlePropagator) 
  ////////////////////////////////////////////////////
  fInputArray   = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // import beamspot
  try{
    fBeamSpotInputArray = ImportArray(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"));
  }
  catch(runtime_error &e){
    fBeamSpotInputArray = 0;
  }
  
  fBz = GetDouble("Bz", 0.0);

  // Name of module useful for debugging 

  // Output arrays 
  fHitOutputArray = ExportArray(GetString("OutputArray", "hits"));

  ///////////////////////////////
  // parameters for barrel layers
  ///////////////////////////////
  fBarrelLength = GetDouble("BarrelLength", 1.0); // barrel length [m]
  ExRootConfParam barrelLayersParam = GetParam("BarrelLayerRadii"); // barrel layer radii [m]

  Long_t size = barrelLayersParam.GetSize();
  for(int i = 0; i < size; ++i){
    //std::cout << "parameter " << i << ": " << barrelLayersParam[i].GetDouble() << std::endl;
    fBarrelLayerRadii.push_back( barrelLayersParam[i].GetDouble() ); 
  }

  ///////////////////////////////
  // parameters for endcap layers
  ///////////////////////////////
  ExRootConfParam endcapZPositions = GetParam("EndCapZ"); // z positions of the endcaps. Only needed in the positive sense (symmetry assumed) 
  fEndCapRadius = GetDouble("EndCapRadius", 1.0); // end cap radius [m] 

  // Make sure the user only put in positive values for endcap
  // If an asymetrical tracker design was needed, then a hack would be required
  for(int i=0; i<endcapZPositions.GetSize(); ++i){
    if( endcapZPositions[i].GetDouble() < 0.0){
      std::cerr << "WARNIGN: HitFinder: Negative endcap z positions detected, which will be ignored. The user need only input positive values, and a symmetrical detector design is assumed." << std::endl;
      continue;
    }
    fEndcapZPositions.push_back( endcapZPositions[i].GetDouble() ); 
  }

  if(m_debug) std::cout << "Init(): Defined input and output arrays" << std::endl;

}

//------------------------------------------------------------------------------

void HitFinder::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void HitFinder::ParticlePropagator(float RADIUS_MAX, float HalfLengthMax, int SurfaceID, bool removeEndcaps, bool removeBarrel)
{

  if(m_debug) std::cout << "ParticlePropagator()" << std::endl;

  Candidate *candidate, *mother;
  TLorentzVector candidatePosition, candidateMomentum, beamSpotPosition;
  Double_t px, py, pz, pt, pt2, e, charge;
  Double_t x, y, z, t, r, phi;
  Double_t x_c, y_c, r_c, phi_c, phi_0;
  Double_t x_t, y_t, z_t, r_t;
  Double_t t1, t2, t3, t4, t5, t6;
  Double_t t_z, t_r, t_ra, t_rb;
  Double_t tmp, discr, discr2;
  Double_t delta, gammam, omega, asinrho;
  Double_t rcu, rc2, xd, yd, zd;
  Double_t pathLength, d0, dz, momentum, ctgTheta, phip, etap, alpha;
  //Double_t bsx, bsy, bsz;


  const Double_t c_light = 2.99792458E8;

  /***************
  if (!fBeamSpotInputArray || fBeamSpotInputArray->GetSize () == 0)
    beamSpotPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  else
  {
    Candidate &beamSpotCandidate = *((Candidate *) fBeamSpotInputArray->At(0));
    beamSpotPosition = beamSpotCandidate.Position;
  }
  ******************/

  // process all candidates 
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;
    x = candidatePosition.X()*1.0E-3;
    y = candidatePosition.Y()*1.0E-3;
    z = candidatePosition.Z()*1.0E-3;

    /**************
     * WJF optimize 
    bsx = beamSpotPosition.X()*1.0E-3;
    bsy = beamSpotPosition.Y()*1.0E-3;
    bsz = beamSpotPosition.Z()*1.0E-3;
    ***************/


    // Check that particle position is inside the cylinder that defines the detector volume 
    // NB TMath::Hypot(x, y) calculates hypotenuse see https://en.wikipedia.org/wiki/Hypot 
    if(TMath::Hypot(x, y) > RADIUS_MAX || TMath::Abs(z) > HalfLengthMax){
      continue;
    }

    charge = candidate->Charge;

    // Were only interested in charged particles, since we're looking for hits 
    if(TMath::Abs(charge) < 1.0E-9) continue; 

    px = candidateMomentum.Px();
    py = candidateMomentum.Py();
    pz = candidateMomentum.Pz();
    pt = candidateMomentum.Pt();
    pt2 = candidateMomentum.Perp2();
    e = candidateMomentum.E();

    if(pt2 < 1.0E-9){
      continue;
    }

    if(TMath::Abs(charge) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
    {
      // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
      tmp = px*y - py*x;
      discr2 = pt2*RADIUS_MAX*RADIUS_MAX - tmp*tmp;

      if(discr2 < 0.0)
      {
        // no solutions
        continue;
      }

      tmp = px*x + py*y;
      discr = TMath::Sqrt(discr2);
      t1 = (-tmp + discr)/pt2;
      t2 = (-tmp - discr)/pt2;
      t = (t1 < 0.0) ? t2 : t1;

      z_t = z + pz*t;
      if(TMath::Abs(z_t) > HalfLengthMax)
      {
        t3 = (+HalfLengthMax - z) / pz;
        t4 = (-HalfLengthMax - z) / pz;
        t = (t3 < 0.0) ? t4 : t3;
      }

      x_t = x + px*t;
      y_t = y + py*t;
      z_t = z + pz*t;


      mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      candidate->InitialPosition = candidatePosition;
      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*e*1.0E3);
      /*********
       * WJF optimize
      pathLength = TMath::Sqrt( (x_t - x)*(x_t - x) + (y_t - y)*(y_t - y) + (z_t - z)*(z_t - z));
      candidate->L = pathLength*1.0E3;
      *
      * ************/

      candidate->Momentum = candidateMomentum;
      candidate->AddCandidate(mother);

      if(TMath::Abs(charge) > 1.0E-9)
      {
        // This section should never be entered, relic from original ParticlePropagator
        switch(TMath::Abs(candidate->PID))
        {
          case 11:
            fHitOutputArray->Add(candidate);
            break;
          case 13:
            fHitOutputArray->Add(candidate);
            break;
          default:
            fHitOutputArray->Add(candidate);
        }
      }
      else
      {
        //fNeutralOutputArray->Add(candidate);
      }
    }
    else
    {

      // 1.  initial transverse momentum p_{T0}: Part->pt
      //     initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
      //     relativistic gamma: gamma = E/mc^2; gammam = gamma * m
      //     gyration frequency omega = q/(gamma m) fBz
      //     helix radius r = p_{T0} / (omega gamma m)

      gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c^2]
      omega = charge * fBz / (gammam);                // omega is here in [89875518/s]
      r = pt / (charge * fBz) * 1.0E9/c_light;        // in [m]

      phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

      // 2. helix axis coordinates
      x_c = x + r*TMath::Sin(phi_0);
      y_c = y - r*TMath::Cos(phi_0);
      r_c = TMath::Hypot(x_c, y_c);
      phi_c = TMath::ATan2(y_c, x_c);
      phi = phi_c;
      if(x_c < 0.0) phi += TMath::Pi();

      rcu = TMath::Abs(r);
      rc2 = r_c*r_c;

      // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
      xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
      xd = (rc2 > 0.0) ? xd / rc2 : -999;
      yd = y_c*(-rcu*r_c + rc2);
      yd = (rc2 > 0.0) ? yd / rc2 : -999;
      zd = z + (TMath::Sqrt(xd*xd + yd*yd) - TMath::Sqrt(x*x + y*y))*pz/pt;

      // use perigee momentum rather than original particle
      // momentum, since the orignal particle momentum isn't known

      px = TMath::Sign(1.0, r) * pt * (-y_c / r_c);
      py = TMath::Sign(1.0, r) * pt * (x_c / r_c);
      etap = candidateMomentum.Eta();
      phip = TMath::ATan2(py, px);

      candidateMomentum.SetPtEtaPhiE(pt, etap, phip, candidateMomentum.E());

      // calculate additional track parameters (correct for beamspot position)

      /****************
       *
       * WJF opimization
       *
      d0        = ((x - bsx) * py - (y - bsy) * px) / pt;
      dz        = z - ((x - bsx) * px + (y - bsy) * py) / pt * (pz / pt);
      momentum  = candidateMomentum.P();
      ctgTheta  = 1.0 / TMath::Tan (candidateMomentum.Theta());
      * 
      *
      ***************/


      // 3. time evaluation t = TMath::Min(t_r, t_z)
      //    t_r : time to exit from the sides
      //    t_z : time to exit from the front or the back
      t_r = 0.0; // in [ns]
      int sign_pz = (pz > 0.0) ? 1 : -1;
      if(pz == 0.0) t_z = 1.0E99;
      else t_z = gammam / (pz*1.0E9/c_light) * (-z + HalfLengthMax*sign_pz);

      if(r_c + TMath::Abs(r)  < RADIUS_MAX)
      {
        // helix does not cross the cylinder sides
        t = t_z;
      }
      else
      {
        asinrho = TMath::ASin((RADIUS_MAX*RADIUS_MAX - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c));
        delta = phi_0 - phi;
        if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
        if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
        t1 = (delta + asinrho) / omega;
        t2 = (delta + TMath::Pi() - asinrho) / omega;
        t3 = (delta + TMath::Pi() + asinrho) / omega;
        t4 = (delta - asinrho) / omega;
        t5 = (delta - TMath::Pi() - asinrho) / omega;
        t6 = (delta - TMath::Pi() + asinrho) / omega;

        if(t1 < 0.0) t1 = 1.0E99;
        if(t2 < 0.0) t2 = 1.0E99;
        if(t3 < 0.0) t3 = 1.0E99;
        if(t4 < 0.0) t4 = 1.0E99;
        if(t5 < 0.0) t5 = 1.0E99;
        if(t6 < 0.0) t6 = 1.0E99;

        t_ra = TMath::Min(t1, TMath::Min(t2, t3));
        t_rb = TMath::Min(t4, TMath::Min(t5, t6));
        t_r = TMath::Min(t_ra, t_rb);
        t = TMath::Min(t_r, t_z);
      }

      // 4. position in terms of x(t), y(t), z(t)
      x_t = x_c + r * TMath::Sin(omega * t - phi_0);
      y_t = y_c + r * TMath::Cos(omega * t - phi_0);
      z_t = z + pz*1.0E9 / c_light / gammam * t;
      r_t = TMath::Hypot(x_t, y_t);


      // compute path length for an helix

      /*********
       * WJF optimize 
      alpha = pz*1.0E9 / c_light / gammam;
      l = t * TMath::Sqrt(alpha*alpha + r*r*omega*omega);
      *******/

      if(r_t > 0.0)
      {

        // store these variables before cloning
        /****************
         * WJF: optimization
         *
        candidate->D0 = d0*1.0E3;
        candidate->DZ = dz*1.0E3;
        candidate->P  = momentum;
        candidate->CtgTheta = ctgTheta;
        candidate->Phi = phip;
        *******************/
        candidate->PT = pt;

        mother = candidate;
        candidate = static_cast<Candidate*>(candidate->Clone());

        candidate->InitialPosition = candidatePosition;
        candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, candidatePosition.T() + t*c_light*1.0E3);

        // do this ASAP (performance) 
        if(removeEndcaps){
          if( fabs(candidate->Position.Perp()) < 0.9999*RADIUS_MAX*1E3) continue; // 0.9999 since radii of hit is slightly smaller than radii of surface
        }

        if(removeBarrel){
          if( fabs( candidate->Position.Z() ) < 0.9999*HalfLengthMax*1E3) continue;
        }

        candidate->Momentum = candidateMomentum;
        //candidate->L = pathLength*1.0E3; // WJF optimization
        candidate->Xd = xd*1.0E3;
        candidate->Yd = yd*1.0E3;
        candidate->Zd = zd*1.0E3;
        candidate->AddCandidate(mother);

        // add surface ID to hit
        candidate->SurfaceID = SurfaceID; 

        // Add candidate to output array  
        fHitOutputArray->Add(candidate);

      }
    }
  }
}

void HitFinder::Process()
{

  if(m_debug) std::cout << "Process()" << std::endl;
  int SurfaceID(0);

  // Creat hits for all barrel layers
  clock_t begin = clock();
  int hitNumber(0);
  for(auto barrelRadius : fBarrelLayerRadii){
    std::vector<TLorentzVector> hits;
    bool removeEndcaps(true);
    bool removeBarrel(false);
    ParticlePropagator(barrelRadius, fBarrelLength, SurfaceID, removeEndcaps, removeBarrel); // all barrels have the same length 

    SurfaceID++;
  }
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << " Time for call of ParticlePropagator()x3: " << elapsed_secs << std::endl;

  // Create hits for endcap layers
  for(auto endcapZ : fEndcapZPositions){
    std::vector<TLorentzVector> hits;
    bool removeEndcaps(false);
    bool removeBarrel(true);
    ParticlePropagator(fEndCapRadius, endcapZ, SurfaceID, removeEndcaps, removeBarrel); 
    SurfaceID++;
  }

}

//------------------------------------------------------------------------------
