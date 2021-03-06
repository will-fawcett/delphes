#include "classes/TrackFitter.h" 
#include "classes/HitCollection.h"
#include <cmath>
#include <algorithm>
#include <utility>
#include "TMath.h"

#include <ctime>


inline float quickDeltaPhi(float phi1, float phi2){
  if(phi1<0) phi1+= 2*M_PI;
  if(phi2<0) phi2+= 2*M_PI;
  float dPhi= fabs(phi1 - phi2);
  if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
  return dPhi; 
}

bool sortByHypot( const std::pair< float, float>& a, const std::pair< float, float>& b){
  return (TMath::Hypot(a.first, a.second) < TMath::Hypot(b.first, b.second));
}

inline float quotient(float r, float r2, float param1, float param2){
  return pow(r2,4)*param1*param1 - r*r * r2*r2 *(r2*r2 - 4*param2*param2);
}

// print out map 
void printNewHitMap(std::map<std::string, std::vector<Hit*>> theMap){
  for(const auto& thep : theMap){
    std::cout << "ID: " << thep.first << "\t" << thep.second.size() << std::endl;
  }
}

std::map<std::string, std::vector<Hit*> > TrackFitter::associateHitsSimplePattern(hitContainer& hc, Location& loc) const{

  // Separate collection of hits into "layer-eta"phi" regions

  // Some pre-defined knowledge about the tracker
  const float barrelLength = 2250; // [mm] 

  // WJF: remove set-of-locations functionality
  //std::vector<std::string> setOfLocations; 

  std::map<std::string, std::vector<Hit*> > newMap; 
  for(const auto layer : m_layerIDs){
    for(Hit* hit : hc[layer]){
      std::string locationString = loc.locationFromHit(hit);
      newMap[ locationString ].push_back(hit); 
      //setOfLocations.push_back(locationString);
    }
  }

  // add set of locations to the Location object
  //loc.addSetOfLocationsStrings(setOfLocations);

  // debug 
  //printNewHitMap(newMap); 

  return newMap;
}

std::vector<Hit*> concatenateVector(std::vector<Hit*>& A, std::vector<Hit*>& B){
  // https://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
  std::vector<Hit*> AB;
  AB.reserve( A.size() + B.size() ); // preallocate memory
  AB.insert( AB.end(), A.begin(), A.end() );
  AB.insert( AB.end(), B.begin(), B.end() );
  return AB;
}

std::vector<Hit*> concatenateHitsBasedOnLocations(std::map<std::string, std::vector<Hit*>>& hitMap, std::vector<std::string>& locations){
  // Return a vector of all hits in all of the regions selected by locations 
  // Can probably make this more efficient ? 
  std::vector<Hit*> newVec;
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



std::vector<Hit*> concatenateHitsBasedOnLocations_Jon(const std::map<std::string, std::vector<Hit*>>& hitMap, const std::vector<std::string>& locations) {

  // JB attempt, keeps giving out of range error 
  std::vector<const std::vector<Hit*>*> vectorsToAdd;
  vectorsToAdd.reserve(locations.size() );
  std::size_t output_size(0);
  for (const std::string& location : locations) {
    try{
      vectorsToAdd.push_back(&hitMap.at(location) );
      output_size += vectorsToAdd.back()->size();
    }
    catch(const std::out_of_range& oor){
      std::cout << "out f range error with: " << location << std::endl;
      std::cout << oor.what() << std::endl;
    }
  }
  std::vector<Hit*> outputVec;
  outputVec.reserve(output_size);
  for (const std::vector<Hit*>* vecPtr : vectorsToAdd)
    std::copy(vecPtr->begin(), vecPtr->end(), std::back_inserter(outputVec) );
  return outputVec;
}




bool TrackFitter::associateHitsSimple(hitContainer& hc, float minZ, float maxZ){

    // create a Location object (really just a function ... ) 
    Location loc(0.06, 0.1);
    //loc.printProperties(); 

    // mapping of hits to eta-phi locations 
    std::map<std::string, std::vector<Hit*>> hitMap = this->associateHitsSimplePattern(hc, loc); 

    //////////////////////////////////////////
    // Simplest possible algorithm 
    //////////////////////////////////////////

    // get the inner and outer barrel radii
    const int innerLayerID = m_layerIDs.at(0);
    const int outerLayerID = m_layerIDs.back();
    const int middleLayerID = 1; 
    const float rInner = hc[innerLayerID].at(0)->HitRadius;
    const float rOuter = hc[outerLayerID].at(0)->HitRadius; 

    // reserve some space for the tracks (performance)  
    m_tracks.clear();
    m_tracks.reserve( hc[outerLayerID].size() ); 

    // Calculate the phi window in which the hits in the outer layer must match
    //float trackPtMin = 1.0; // [GeV] (minimum track pT to consider for phiWindow calculation)
    const float trackPtMin = 2.0; // [GeV] (minimum track pT to consider for phiWindow calculation)
    const float bendingRadius = 1000 * trackPtMin/1.199; // [mm]
    float phiWindow = fabs( acos(rInner / (2*bendingRadius)) - acos(rOuter / (2*bendingRadius)) );
    phiWindow *= 2; // multiply by two, to have the deviation travelling in either direction. 


    // Draw a line between the hit in the innermost and outermost layer
    // See if there is a hit on the line in the intermediate layer (within some tolerance)
    for(const auto& innerHit : hc[innerLayerID]){

      const float zInner = innerHit->Z;
      const float phiInner = innerHit->Phi;

      // get locations (areas) for other hits
      std::string innerHitLocation = loc.locationFromHit(innerHit); 
      std::vector<std::string> outerHitLocations  = loc.listOfLocationsInLayer(innerHitLocation, outerLayerID );
      std::vector<std::string> middleHitLocations = loc.listOfLocationsInLayer(innerHitLocation, middleLayerID);
      //std::cout << "There are " << outerHitLocations.size() << " outer hit locations: " << std::endl;
      //for(auto l : outerHitLocations) std::cout << "\t " <<  l << std::endl;

      // get vector of hits defined by list of locations
      std::vector<Hit*> outerHitVector  = concatenateHitsBasedOnLocations(hitMap, outerHitLocations);
      std::vector<Hit*> middleHitVector = concatenateHitsBasedOnLocations(hitMap, middleHitLocations);

      //for(const auto& outerHit : hc[outerLayerID]){
      for(const auto& outerHit : outerHitVector){

        // must be within phi criteria  
        if( quickDeltaPhi(phiInner, outerHit->Phi) > phiWindow) continue; 

        // calculate parameters of line from inner hit to outer hit 
        const float zOuter = outerHit->Z;
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

          float middleZresiduum = intermediateHit->Z - intersect;
          // only select if intermediate hit matches within tolerance along Z  
          if(fabs(middleZresiduum) < m_tolerance){

            // reject the intermediate hit if it is also outside the phi window
            if( quickDeltaPhi(phiInner, intermediateHit->Phi) > phiWindow) continue; 

            // Three hits are matched -> a track 
            std::vector<Hit*> matchedHits;
            matchedHits.reserve(4); // prevent vector from having to grow 
            matchedHits.push_back(innerHit);
            matchedHits.push_back(intermediateHit);
            matchedHits.push_back(outerHit);

            /****************
            std::cout << "Track creation" << std::endl;
            std::cout << "PHI: Inner: " << innerHit->Phi << "\tMiddle: " << intermediateHit->Phi << "\touter: " << outerHit->Phi << std::endl;
            std::cout << "R:   Inner: " << innerHit->HitRadius << "\tMiddle: " << intermediateHit->HitRadius << "\touter: " << outerHit->HitRadius << std::endl;
            std::cout << "Layer:   Inner: " << innerHit->SurfaceID << "\tMiddle: " << intermediateHit->SurfaceID << "\touter: " << outerHit->SurfaceID << std::endl;
            *****************/

            myTrack aTrack(matchedHits, beamlineIntersect, middleZresiduum);
            m_tracks.push_back( aTrack ); 
            //m_tracks.emplace_back(  ); 
          }
        }
      }
    }

    return true;
}

float TrackFitter::calculateRPhiWindowInToOut(const float a, const float b, const float r){
  /***********
   * Calculates the maximum deviation in phi that a particle could have traversed 
   * when travelling from a hit point (a, b) in an inner barrel layer to an outer barrel 
   * layer with radius r 
   * ************************/

  // Radius of trajectory
  float R = r/2.0;

  // centres of the circle (alpha1, beta2) and (alpha2, beta1)
  float quotient = sqrt( (4*R*R - a*a - b*b) * (a*a + b*b) );
  float q = quotient / (a*a + b*b); 

  float alpha1 = 0.5 * ( a + b*q );
  float alpha2 = 0.5 * ( a - b*q );

  float beta1 = 0.5 * ( b + a*q );
  float beta2 = 0.5 * ( b - a*q );

  // angles
  float phi1 = atan2(beta1, alpha2);

  float phi2 = atan2(beta2, alpha1);

  if(m_debug){
    std::cout << "a, b, r : " << a << " " << b << " " << r  << std::endl;
    std::cout << "(a1, b2) (" << alpha1 << ", " << beta2 << ") : phi = " << phi2/M_PI << "pi" << std::endl;
    std::cout << "(a2, b1) (" << alpha2 << ", " << beta1 << ") : phi = " << phi1/M_PI << "pi" <<  std::endl;
    std::cout << "" << std::endl; 
  }
  
  // return the max deviation
  float deltaPhi = acos(cos(phi1 - phi2))/2.0; 
  if( isnan(deltaPhi) ){
    std::cerr << "ERROR: calculateRPhiWindowInToOut() deltaPhi is NAN!" << std::endl;
    return 0;
  }

  return deltaPhi;

}

bool TrackFitter::calculateRPhiWindowOutToIn(const float r2, const float a, const float b){
  /***********
   * Calculates the coordinates of the intersection of:
   * - a circle defined as touching the origin and point (a,b)
   * - with another circle of radius r2 which is centered at the origin
   * Note that: r2^2 < a^2 + b^2
   * ************************/

  // check geometry of coordinates
  if ((a*a + b*b) < r2*r2){
    std::cerr << "ERROR: hit coordinates are inside the barrel layer" << std::endl;
    return false;
  }

  // Calculate center of circle (0, 0) -- (a, b)
  float alpha = a/2;
  float beta  = b/2; 
  // radius of the circle
  float rad = sqrt(alpha*alpha + beta*beta); 

  // coordinates of the intersection
  float x1 = (r2*r2*alpha + sqrt( quotient(rad, r2, alpha, beta) )) / (2*rad*rad);
  float x2 = (r2*r2*alpha - sqrt( quotient(rad, r2, alpha, beta) )) / (2*rad*rad);

  float y1 = (r2*r2*beta + sqrt( quotient(rad, r2, beta, alpha) )) / (2*rad*rad);
  float y2 = (r2*r2*beta - sqrt( quotient(rad, r2, beta, alpha) )) / (2*rad*rad);

  // dont know which comination are on the circle with radius r2, check
  /******************
   * WJF: will probably have an analytic solition! 
  std::cout << "rad: " << rad << std::endl;  
  float h1 = TMath::Hypot(x1, y1);
  float h2 = TMath::Hypot(x1, y2);
  float h3 = TMath::Hypot(x2, y1);
  float h4 = TMath::Hypot(x2, y2);
  std::cout << "(x1, y1) : h1 = " << h1 << std::endl;
  std::cout << "(x2, y1) : h2 = " << h2 << std::endl;

  std::vector< std::pair<float, float> > combinations; // could use map, but with integer index so it's the same as a vector
  combinations.push_back( std::make_pair( x1, y1 ) ); 
  combinations.push_back( std::make_pair( x1, y2 ) ); 
  combinations.push_back( std::make_pair( x2, y1 ) ); 
  combinations.push_back( std::make_pair( x2, y2 ) ); 

  // sort the combinations 
  // TODO: Write correct sorting angle
  std::sort(combinations.begin(), combinations.end(), sortByHypot); 

  //
  std::pair<float, float> c1, c2;
  c1 = combinations.at(1);
  c2 = combinations.at(2); 

  // now find angles 
  // TODO: check this gives the correct angle (!) 
  float phi1 = atan2( c1.second, c1.first ) ;
  float phi2 = atan2( c2.second, c2.first ) ;
  *****************************/

  // return these angles

  return true;  
}



float TrackFitter::calculateZWindowForNextLevel(float y0, float x0, float y2, float x1){
  /***********************************************
   * Calculate the parameters of the straight line passing through the coordinates (x0, y0), (x1, y1) where y1=0
   * Return the x coordinate of the line at y2
   * ********************************************/
  float y1 = 0.0;
  
  // line parameters
  LineParameters params;
  params.calculateLineParameters(x0, y0, x1, y1);

  return (y2 - params.y_intercept())/params.gradient(); 
}

bool TrackFitter::associateHitsLinearOutToIn(hitContainer hitMap, float minZ, float maxZ){
  /***********************************************
   * Hit association algorithm
   * - Assumes concentric barrel layers, no endcaps
   * - Innermost layer labelled with 0, incremented as layers increase
   *
   * The algorithm:
   * - Starts from a hit in the outer layer, calculates a window in the layer benieth in which to search for other hits
   * - Loops over hits in the next innermost layer, if within the search window, then this hit is assigned the hit in the above layer
   *   -- repeat untill the last layer
   * *********************************************/

  // Reverse layer IDs 
  std::reverse(m_layerIDs.begin(), m_layerIDs.end()); // should count from 3 .. 2 .. 1 

  /*****************8
  // loop over all hits in outer layer
  std::vector<Hit*> hitCollectionVec; 
  for(Hit * outerHit : hitMap[ layers.at(0) ]){
    HitCollection collection(outerHit);

    // match to any hits in subsequent layers
    for(int layerID : layers){

      // skip the outer layer (alredy done)
      if(layerID == layers.at(0)) continue; 

      // get the r coordinate of the next-innermost layer
      float rInner = hitMap[layerID-1].at(0)->HitRadius; // should be able to optimize this away? 

      // loop over all hits in the inner layers
      for(Hit * innerHit : hitMap[layerID]){

        // find window for hits in the next layer to be assigned to this one
        float r = hit->HitRadius;
        float z = hit->Z;
        float zLeft  = calculateZWindowForNextLevel(r, z, rInner, minZ); 
        float zRight = calculateZWindowForNextLevel(r, z, rInner, maxZ); 

        if(zLeft < innerHitZ && innerHitZ < zRight){
          collection->addHit(innerHit); 
        }
      }
  }
  **********************/
  

  /******************************8
  // start from outermost barrel layer, and work inwards
  std::map<int, std::vector<HitCollection> > newHitMap; 
  std::vector<HitCollection> hitCollections; 
  for(int layerID : layers){

    if(layerID == layers.back()) continue; // don't execute algorithm for innermost layer

    // get the r coordinate of the next-innermost layer
    float rInner = hitMap[layerID-1].at(0)->HitRadius; // should be able to optimize this away? 

    // loop over all hits in layer
    for(Hit* hit : hitMap[layerID]){

      HitCollection collection(hit);
      newHitMap[layerID].push_back( collection );

      // find window for hits in the next layer to be assigned to this one
      float r = hit->HitRadius;
      float z = hit->Z;
      float zLeft  = calculateZWindowForNextLevel(r, z, rInner, minZ); 
      float zRight = calculateZWindowForNextLevel(r, z, rInner, maxZ); 

      // loop over all hits in next layer
      for(Hit* innerHit : hitMap[layerID - 1]){
        float innerHitZ = innerHit->Z; 

        // if hit within window, assign to the bit above 
        if(zLeft < innerHitZ && innerHitZ < zRight){
          collection.addHit(innerHit); 
        }
      }
    } // loop over hits in layer
  }
  ***********************/

  // extract all hits again
  //std::vector<HitCollection> allHits;
  for(auto layer : m_layerIDs){
    for(Hit* hit : hitMap[layer]){
      //m_associatedHitCollection.push_back(HitCollection(hit));
      m_associatedHitCollection.push_back(HitCollection(hit));
    }
  }

  // surfaces labelled 2 .. 1 ... 0. 0 being innermost 
  for(auto & hit : m_associatedHitCollection){
    int layerID = hit.SurfaceID;
    if(layerID == m_layerIDs.back()) continue; // no more layers inside 

    // calculate r of next-lower level
    float rInner = hitMap[layerID-1].at(0)->HitRadius; 

    // calculate search window in z
    float r = hit.Perp();
    float z = hit.Z();
    float zLeft  = calculateZWindowForNextLevel(r, z, rInner, minZ); 
    float zRight = calculateZWindowForNextLevel(r, z, rInner, maxZ); 

    // calculate search window in phi
    //calculateRPhiWindowOutToIn()

    for(HitCollection jhit : m_associatedHitCollection){
      // only assign if hit is in next-lower level
      if(jhit.SurfaceID != layerID-1) continue;
      if(zLeft < jhit.Z() and jhit.Z() < zRight){
        // add hit 
        hit.addHit(&jhit);
      }
    }
  }


  return true;
} // end associateHitsLinearOutToIn


bool TrackFitter::associateHitsLinearInToOut(hitContainer hitMap, float minZ, float maxZ){
  /****************************************************
   * Fills the m_associatedHitCollection member variable, which is of type std::vector<HitCollection>
   * Each element of this vector is a HitCollection object
   *
   * For each hit in the innermost layer, r--z and r--phi matching is performed. 
   * All hits that lie within this volume are associated to the innermost hit (via the HitCollection object)
   * (the vector will then have length equal to that of the number of hits in the innermost layer)
   * ***********************************************/


  /////////////////////
  // Fill HitCollection
  /////////////////////
  
  // first reserve some space (performance)
  int numHits(0);
  for(auto layer : m_layerIDs){
    numHits += hitMap[layer].size();
  }
  m_associatedHitCollection.reserve(numHits+1);


  // Fill hits 
  for(auto layer : m_layerIDs){
    for(Hit* hit : hitMap[layer]){
      m_associatedHitCollection.push_back(HitCollection(hit));
    }
  }

  // deltaPhi parameter for the search window. Only depends on the separation of the inner and outermost layers
  // radius of outermost barrel layer
  float outermostRadius = hitMap[ m_layerIDs.back() ].at(0)->HitRadius;
  Hit* anInnerHit = hitMap[ 0 ].at(0);
  float maxPhiDeviation = calculateRPhiWindowInToOut(anInnerHit->X, anInnerHit->Y, outermostRadius); // may only be necessary once 
  float innermostRadius = anInnerHit->HitRadius;

  // associate hit algorithm 
  // It's very important to loop over the references to the object
  // if a copy is made, then the pointer storing below wouldn't work
  for(HitCollection& innerHit : m_associatedHitCollection){

    // first loop over inner hits 
    if(innerHit.SurfaceID != m_layerIDs.at(0)) continue; 

    // parameters of for the search window
    float z = innerHit.Z();
    float zRight  = calculateZWindowForNextLevel(innermostRadius, z, outermostRadius, minZ); 
    float zLeft = calculateZWindowForNextLevel(innermostRadius, z, outermostRadius, maxZ); 

    // second loop over outer hits 
    for(HitCollection& outerHit : m_associatedHitCollection){
      if(outerHit.SurfaceID == m_layerIDs.at(0)) continue; // skip first layer (only want seeds from first layer)

      // determine if matched
      if( (zLeft < outerHit.Z() && outerHit.Z() < zRight) ){ // split into two loops to save on DeltaPhi calculation

        // calculate search window
        float phiDeviation = innerHit.DeltaPhi(outerHit);

        if( phiDeviation < maxPhiDeviation ){
          innerHit.addHit(&outerHit); // we're storing pointers to object stored in anoter vector
          // usually it's a bad idea to store pointers to objects defined on the stack
          // but m_associatedHitCollection won't be modified until it's deleted 
        }
      }
    } // loop over outer hits
  } // loop over inner hits



  return true;
}


bool TrackFitter::combineHitsToTracksMatchingInnerAndOutermost(){
  /**********************************
   * This function creates tracks from hit collections. 
   *
   * The algorithm proceeds as follows:
   * For each HitCollection object (stored inside m_associatedHitCollection), all possible combinations of the innermost and outermost hit are made.
   * A straight line is then drawn between the two points. 
   * If in _all_ intermediate layers, there is also a hit that lies on the line (within some tolerance), then the straight line is returned as a track
   *
   * ********************************/


  float tolerance = 0.4; // mm  2*sqrt(0.04) = 0.4 , 0.04 mm is size of pixel 

  // Create pairs of innermost and outermost hits
  int innerLayerID = m_layerIDs.at(0);
  int outerLayerID = m_layerIDs.back();


  if(innerLayerID == outerLayerID){
    std::cerr << "ERROR: Only on layer! Impossible to form tracks" << std::endl;
    return false; 
  }

  // Number of intermediate layers
  int numIntermediateLayers = m_layerIDs.size() - 2;
  
  if(m_debug){
    std::cout << "begin combineHitsToTracksMatchingInnerAndOutermost():" << std::endl;
    std::cout << "inner layer is : " << innerLayerID << std::endl;
    std::cout << "outer layer is : " << outerLayerID << std::endl; 
    std::cout << "There are " << numIntermediateLayers << " intermediate layers" << std::endl;
    std::cout << "There are " << m_associatedHitCollection.size() << " hits in the innermost layer." << std::endl;
  }

  // Loop over all matched hit collections
  for(const HitCollection& hitCollection : m_associatedHitCollection){
    if(hitCollection.SurfaceID != innerLayerID) continue; // just to make sure were really talking about hit collections in the inner layer 
    if(m_debug) hitCollection.printMatchedHits();

    std::map<int, std::vector<Hit*>> layerToHitMap = hitCollection.getLayerToHitMap();
    for(Hit* innerHit : layerToHitMap[innerLayerID]){
      int numTracksCreatedFromThisHit(0);
      for(Hit* outerHit : layerToHitMap[outerLayerID]){

        LineParameters trackLine;
        trackLine.calculateLineParameters( innerHit->Z, innerHit->HitRadius, outerHit->Z, outerHit->HitRadius );

        if(m_debug){
          std::cout << "Inner hit (x, y, z): " << innerHit->X << ", " << innerHit->Y << ", " << innerHit->Z << " . r =" << innerHit->HitRadius << std::endl;
          std::cout << "outer hit (x, y, z): " << outerHit->X << ", " << outerHit->Y << ", " << outerHit->Z << " . r =" << outerHit->HitRadius << std::endl;
          std::cout << "Line parameters, grad: " << trackLine.gradient() << " \tintercept: " << trackLine.x_intercept() << std::endl;
        }

        /****************************************
         * in order to find the tolerance within which hits in the inner layers can be associated
         * vary the coordinates of the hit by 2*sqrt(pixel size) in the inner and outer layers
         * two new lines are constructed, which creates a window in which the hit inside the inner layers could be
         * if each of the inner layers has a hit, then it's a track!
         * 
         *            outer hit +/- 2*sqrt(pixel size)
         * ----------- |---o---|----------------
         *           /   /   /
         * ---------|-------|-------------------  search window (should also be 2*sqrt(pixel size)) 
         *         /   /   /
         * -------|---o---|---------------------
         *        inner hit +/- 2*sqrt(pixel size)
         *
         * Since this tolerance is just 2*sqrt(pixel size) then we don't need to calculate these extra lines
         * but would do if the geometry was different 
         *
         ****************************************/

        // Detect if there is a hit in all of the intermediate layers (within tolerance), return track if true 
        int nMatchedHits(0);
        std::vector<Hit*> trackHits;
        trackHits.push_back(innerHit);

        // loop over intermediate hits (could be made more efficient by calculating radius and zintercept only once?)
        for(auto layer : m_layerIDs){
          if(layer == innerLayerID || layer == outerLayerID) continue; // only intermediate layers

          if(m_debug) std::cout << "There are " << layerToHitMap[layer].size() << " loosely matched hits in the intermediate layer" << std::endl;
          if(layerToHitMap[layer].size() == 0) continue; // must have at least one loosely matched hit in the intermediate layers

          // radius of this intermediate layer
          float radius = layerToHitMap[layer].at(0)->HitRadius; 

          // z coordinate of intersection of line and this barrel layer
          float zCoordinate = (radius - trackLine.x_intercept()) / trackLine.gradient();

          for(Hit* intermediateHit : layerToHitMap[layer]){

            if(fabs(intermediateHit->Z - zCoordinate) < tolerance){
              if(m_debug) std::cout << "\tThere is a matched hit!" << std::endl;
              trackHits.push_back( intermediateHit );
              nMatchedHits++;
            }
          }
        } // end loop over intermediate layers 
        trackHits.push_back( outerHit );

        // if there is a correct number of matches, add the track!
        if(nMatchedHits == numIntermediateLayers){

          // Now create a straight line with all points (re-fit the track with the extra hits) 
          myTrack aTrack(trackHits, 0, 0);
          m_tracks.push_back( aTrack ); 
          numTracksCreatedFromThisHit++;
        }

      } // end loop over hits in outer layer
      if(m_debug){
        std::cout << "Number of tracks created from this hit: " << numTracksCreatedFromThisHit << std::endl;
        GenParticle * particle = (GenParticle*) innerHit->Particle.GetObject();
        std::cout << "particle had pT=" << particle->PT << " GeV" << std::endl;

        
      }
    } // end loop over hits in inner layer 
    if(m_debug) std::cout << "" << std::endl;
  } // end loop over hit collections

  return true;
}

bool TrackFitter::combineHitsToTracksInToOut(){
  /**********************************
  * Create tracks from hit collections.
  * This algorithm works as follows. 
  * For each HitCollection object (stored inside m_associatedHitCollection)
  * all possible combinations of hits are made in each layer. Each set then contains one hit in each layer. 
  * A straight-line track is then created for each set of hits using least squares to extract the line parameters.
  * If the line does not intersect the beamline within 3-sigma of the luminous region, then the track is rejected.
  ************************************/

  // loop over inner hits 
  for(const auto& innerHit : m_associatedHitCollection){
    if(innerHit.SurfaceID != 0) continue;

    // remove cases where there are fewer than two hits
    if(innerHit.countAssignedHits() < 2) continue; 

    // generate sets of three hits
    std::vector< std::vector<Hit*> > hitCombinations = innerHit.makeHitCollection();
    if(m_debug) std::cout << "This hit has " << hitCombinations.size() << " combinations." <<  std::endl;

    // for each combination, do a straight line fit 
    for(auto& combination : hitCombinations){
      if(m_debug) std::cout << "This combination has " << combination.size() << " hits." << std::endl;


      myTrack aTrack(combination, 0, 0);

      // if track has z0 outside of the defined window, reject the track
      float maxZ = m_parameters.at(1);
      if(fabs(aTrack.Z0()) > maxZ){
        //std::cout << "Track found outside luminous region" << std::endl;
        continue;
      }

      m_tracks.push_back(aTrack); // add to collection of tracks 
    }
  }
  return true;
}


// algorithm not yet written
bool TrackFitter::combineHitsToTracksOutToIn(){
  return false;
}


bool TrackFitter::AssociateHits(hitContainer& hc){

  // Associate hits together such that tracks can be formed
  // The hit association algorithm is determined 
  
  // Make sure these are empty, in case the function is called with a different algorithm 
  m_associatedHitCollection.clear();
  m_tracks.clear();


  switch (fitType) { 
    case linearOutToIn:{
      float minZ = m_parameters.at(0);
      float maxZ = m_parameters.at(1);
      if(this->associateHitsLinearOutToIn(hc, minZ, maxZ)){
        return this->combineHitsToTracksOutToIn();
      }
      else{
        return false;
      }
    }
    case linearInToOut:{
      float minZ = m_parameters.at(0);
      float maxZ = m_parameters.at(1);
      if(this->associateHitsLinearInToOut(hc, minZ, maxZ)){
        return this->combineHitsToTracksInToOut();
      }
      else{
        return false;
      }
    }
    case linearInnerAndOuterMatching:{
       float minZ = m_parameters.at(0);
       float maxZ = m_parameters.at(1);
       if(this->associateHitsLinearInToOut(hc, minZ, maxZ)){
         return this->combineHitsToTracksMatchingInnerAndOutermost();
       }
       else return false;
    }
    case simpleLinear:{
       float minZ = m_parameters.at(0);
       float maxZ = m_parameters.at(1);
       return this->associateHitsSimple(hc, minZ, maxZ);
    }
    case MAX:{
      return false;
    }

  }
  return true; 
}

std::vector <myTrack> TrackFitter::GetTracks(){
  return m_tracks;
}


void TrackFitter::ApplyCurvatureCut(float cut){
  // loop over all tracks, remove those with curvature difference deemed too large
  std::vector< myTrack > newVec;
  for(const myTrack& track: m_tracks){
    if( fabs(track.kappa_bc() - track.kappa_nbc()) < cut){
      newVec.push_back(track);
    }
  }
  m_tracks = newVec; 
}

void TrackFitter::ApplyPtDependantCurvatureCut(std::vector<float> pT_thresholds, std::vector<float> kappaThresholds){

  // loop over all tracks, remove those with curvature less than cut for applied pT threshold
  std::vector< myTrack > newVec;
  for(const myTrack& track : m_tracks){
    float deltaKappa = fabs(track.kappa_bc() - track.kappa_nbc()); 
    bool trackPasses(true);

    // test if the track passes the cuts  
    for(int i=0; i<pT_thresholds.size(); ++i){
      if( track.Pt() > pT_thresholds.at(i) ){
        trackPasses = track.testKappaThreshold(kappaThresholds.at(i));
      }
    }
    if(trackPasses) newVec.push_back(track);
  }
  m_tracks = newVec; 
}
