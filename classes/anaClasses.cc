#include "classes/anaClasses.h" 
#include "classes/HitCollection.h"
#include <cmath>
#include <utility>
#include "TMath.h"
#include "TGraphErrors.h"

#include <ctime>


bool sortByHypot( const std::pair< float, float>& a, const std::pair< float, float>& b){
  return (TMath::Hypot(a.first, a.second) < TMath::Hypot(b.first, b.second));
}

inline float quotient(float r, float r2, float param1, float param2){
  return pow(r2,4)*param1*param1 - r*r * r2*r2 *(r2*r2 - 4*param2*param2);
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

inline lineParameters calculateLineParameters(float x0, float y0, float x1, float y1){
  // Calculate the parameters of the straight line passing through the coordinates (x0, y0), (x1, y1)

  float m = (y0 - y1) / (x0 - x1); 
  float c = (x1*y0 - x0*y1) / (x1 - x0);

  lineParameters params;
  params.gradient = m;
  params.intercept = c; 
  return params;  

}


float TrackFitter::calculateZWindowForNextLevel(float y0, float x0, float y2, float x1){
  /***********************************************
   * Calculate the parameters of the straight line passing through the coordinates (x0, y0), (x1, y1) where y1=0
   * Return the x coordinate of the line at y2
   * ********************************************/
  float y1 = 0.0;
  
  // line parameters
  lineParameters params = calculateLineParameters(x0, y0, x1, y1);

  return (y2 - params.intercept)/params.gradient; 
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
      float rInner = hitMap[layerID-1].at(0)->Perp(); // should be able to optimize this away? 

      // loop over all hits in the inner layers
      for(Hit * innerHit : hitMap[layerID]){

        // find window for hits in the next layer to be assigned to this one
        float r = hit->Perp();
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
    float rInner = hitMap[layerID-1].at(0)->Perp(); // should be able to optimize this away? 

    // loop over all hits in layer
    for(Hit* hit : hitMap[layerID]){

      HitCollection collection(hit);
      newHitMap[layerID].push_back( collection );

      // find window for hits in the next layer to be assigned to this one
      float r = hit->Perp();
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
    float rInner = hitMap[layerID-1].at(0)->Perp(); 

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
  float outermostRadius = hitMap[ m_layerIDs.back() ].at(0)->Perp();
  Hit* anInnerHit = hitMap[ 0 ].at(0);
  float maxPhiDeviation = calculateRPhiWindowInToOut(anInnerHit->X, anInnerHit->Y, outermostRadius); // may only be necessary once 
  float innermostRadius = anInnerHit->Perp();

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
    std::cout << "inner layer is : " << innerLayerID << std::endl;
    std::cout << "outer layer is : " << outerLayerID << std::endl; 
    std::cout << "There are " << numIntermediateLayers << " intermediate layers" << std::endl;
    std::cout << "There are " << m_associatedHitCollection.size() << " hits in the innermost layer." << std::endl;
    std::cout << "" << std::endl;
  }

  // Loop over all matched hit collections
  for(const HitCollection& hitCollection : m_associatedHitCollection){
    if(hitCollection.SurfaceID != innerLayerID) continue; // just to make sure were really talking about hit collections in the inner layer 
    if(m_debug) hitCollection.printMatchedHits();

    std::map<int, std::vector<Hit*>> layerToHitMap = hitCollection.getLayerToHitMap();
    for(Hit* innerHit : layerToHitMap[innerLayerID]){
      for(Hit* outerHit : layerToHitMap[outerLayerID]){

        lineParameters trackLine = calculateLineParameters( innerHit->Z, innerHit->Perp(), outerHit->Z, outerHit->Perp() );

        if(m_debug){
          std::cout << "Inner hit: " << innerHit->X << ", " << innerHit->Y << ", " << innerHit->Z << std::endl;
          std::cout << "outer hit: " << outerHit->X << ", " << outerHit->Y << ", " << outerHit->Z << std::endl;
          std::cout << "Line parameters, grad: " << trackLine.gradient << " \tintercept: " << trackLine.intercept << std::endl;
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
          float radius = layerToHitMap[layer].at(0)->Perp(); 

          // z coordinate of intersection of line and this barrel layer
          float zCoordinate = (radius - trackLine.intercept) / trackLine.gradient; 

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
          myTrack aTrack = simpleLinearLeastSquaresFit(trackHits);
          m_tracks.push_back( aTrack ); 
        }

      } // end loop over hits in outer layer
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

      myTrack aTrack = simpleLinearLeastSquaresFit(combination);

      // if track has z0 outside of the defined window, reject the track
      float maxZ = m_parameters.at(1);
      if(fabs(aTrack.z0) > maxZ){
        //std::cout << "Track found outside luminous region" << std::endl;
        continue;
      }
      m_tracks.push_back(aTrack); // add to collection of tracks 
    }
  }
}

myTrack TrackFitter::simpleLinearLeastSquaresFit(std::vector< Hit* > hits) const {
  std::vector< std::pair<float, float> > coordinates;
  for(Hit* hit : hits){
    coordinates.push_back( std::make_pair(hit->Z, hit->Perp()) ); // coordinates in (r, z)
  }
  lineParameters parameters = simpleLinearLeastSquaresFit(coordinates); 
  return myTrack(parameters, hits); 
}


//TrackFitter::simpleLinearLeastSquaresFit(std::vector<float> xvals, std::vector<float> yvals) const{
lineParameters TrackFitter::simpleLinearLeastSquaresFit(std::vector<std::pair<float, float> > coordinates) const{
  // Function to do simple linear least squares fitting (for a straight line)
  // with the parameters y = mx + c. 
  // Extracts the best fit for m and c. 
  // Takes a vector of the coordinates {xi, yi} 

  float X(0), Y(0), XX(0), XY(0);
  float n = coordinates.size();

  if(m_debug) std::cout << "coordinates (x, y)" << std::endl;
  for(auto coord : coordinates){
    float xi = coord.first;
    float yi = coord.second;
    if(m_debug) std::cout << xi << ", " << yi << std::endl;
    X += xi;
    Y += yi;
    XX += xi*xi;
    XY += xi*yi;
  }

  // gradient
  float m = (XY*n - X*Y) / (XX*n -X*X);

  // intercept
  float c = (Y - m*X) / n; 

  if(m_debug){
    std::cout << "X: " << X << " Y: " << Y << " XX: " << XX << " XY:" << XY << std::endl;
    std::cout << "(m, c) : " << m << ", " << c << ")" << std::endl; 
  }

  lineParameters params;
  params.gradient = m;
  params.intercept = c; 
  return params;  
}

// algorithm not yet written
bool TrackFitter::combineHitsToTracksOutToIn(){
  return false;
}


bool TrackFitter::AssociateHits(hitContainer hc){
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
    case MAX:{
      return false;
    }

  }
  return true; 
}

std::vector <myTrack> TrackFitter::GetTracks(){
  return m_tracks;
}




void myTrack::calculateTrackParameters( cartesianCoordinate coord ){
  // function to calculate d0, z0 (and in principle theta, phi, qOverP)
  // relative to the specified coordinate!!!!!
  // TODO: add calculation if with respect to NOT (0, 0, 0)
  
  // parameters relative to (0, 0, 0)
  phi = atan2( (m_gradient + m_intercept), 1);
  z0 = -1 * (m_intercept / m_gradient);

  d0=0; // can't calculate this with straight line
  theta=0; // not calculable with 2D line
  qOverP=0; // straight line : infinite radius of curvature (therefore infinite momentum)  

}
