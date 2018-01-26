#include "classes/anaClasses.h" 
#include <cmath>
#include "TMath.h"

void sort(std::vector< std::pair< float, float> >){
  x = p.first
  y = p.second
  h = TMath::Hypot(x, y);
  h < something return; 
}

inline float quotient(float r, float r2, float param1, float param2){
  return pow(r2,4)*param1*param1 - r*r * r2*r2 *(r2*r2 - 4*param2*param2);
}

bool TrackFitter::calculateRPhiWindow(const float r2, const float a, const float b){
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
  sort(combinations); 

  //
  std::pair<float, float> c1, c2;
  c1 = combinations.at(1);
  c2 = combinations.at(2); 

  // now find angles 
  // TODO: check this gives the correct angle (!) 
  float phi1 = atan( c1.second, c1.first ) ;
  float phi2 = atan( c2.second, c2.first ) ;

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
  float m = (y0 - y1) / (x0 - x1); 
  float c = (x1*y0 - x0*y1) / (x1 - x0);

  // x = (y-c)/m
  return (y2 - c)/m; 
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

  // Get layer IDs from hitContainer
  std::vector<int> layers;
  for(auto const& key : hitMap){
    layers.push_back(key.first);
  }
  std::reverse(layers.begin(), layers.end()); // should count from 3 .. 2 .. 1 

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
  for(auto layer : layers){
    for(Hit* hit : hitMap[layer]){
      //m_associatedHitCollection.push_back(HitCollection(hit));
      m_associatedHitCollection.push_back(HitCollection(hit));
    }
  }

  // surfaces labelled 2 .. 1 ... 0. 0 being innermost 
  for(auto & hit : m_associatedHitCollection){
    int layerID = hit.SurfaceID;
    if(layerID == layers.back()) continue; // no more layers inside 

    // calculate r of next-lower level
    float rInner = hitMap[layerID-1].at(0)->Perp(); 

    // calculate search window in z
    float r = hit.Perp();
    float z = hit.Z();
    float zLeft  = calculateZWindowForNextLevel(r, z, rInner, minZ); 
    float zRight = calculateZWindowForNextLevel(r, z, rInner, maxZ); 

    // calculate search window in phi
    calculateRPhiWindow()
bool TrackFitter::calculateRPhiWindow(float r2, float a, float b){



    for(HitCollection jhit : m_associatedHitCollection){
      // only assign if hit is in next-lower level
      if(jhit.SurfaceID != layerID-1) continue;
      if(zLeft < jhit.Z() and jhit.Z() < zRight){
        // add hit 
        hit.addHit(&jhit);
      }
    }
  }

  for(auto hit : m_associatedHitCollection){
    hit.printHit();
  }

  return true;
} // end associateHitsLinearOutToIn


bool TrackFitter::associateHitsLinearInToOut(hitContainer hitMap, float minZ, float maxZ){

  // get all barrel layers
  std::vector<int> layers;
  for(auto const& key : hitMap){
    layers.push_back(key.first);
  }

  // Fill HitCollection
  for(auto layer : layers){
    for(Hit* hit : hitMap[layer]){
      m_associatedHitCollection.push_back(HitCollection(hit));
    }
  }

  // associate hit algorithm 
  for(HitCollection& innerHit : m_associatedHitCollection){

    if(innerHit.SurfaceID == layers.at(0)){

      // parameters of for the search window
      float r = innerHit.Perp();
      float z = innerHit.Z();

      for(auto outerHit : m_associatedHitCollection){

        if(outerHit.SurfaceID == layers.at(0)) continue; // skip first layer (only want seeds from first layer)

        // calculate search window
        float rOuter = outerHit.Perp();
        float zRight  = calculateZWindowForNextLevel(r, z, rOuter, minZ); 
        float zLeft = calculateZWindowForNextLevel(r, z, rOuter, maxZ); 

        //TODO add r-phi matching (!)
        std::pair<float, float> phiCoords = calculateRPhiWindow()


        // determine if matched
        if(zLeft < outerHit.Z() && outerHit.Z() < zRight){
          innerHit.addHit(&outerHit);
        }
      }
    }
  }

  return true;
}

// create tracks from hit collections 
bool TrackFitter::combineHitsToTracksInToOut(){

  for(auto hitCollection : m_associatedHitCollection){

    // consider innermost hits 
    if(hitCollection.SurfaceID == 0){
      hitCollection.printHit();

      // now make some tracks ... TODO

    }
  }
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
    case MAX:{
      return false;
    }

  }
  return true; 
}

std::vector <myTrack> TrackFitter::GetTracks(){

  // perform linear fit on the hit collections 

  std::vector<myTrack> temp; 
  return temp;
}


void HitCollection::printAssignedHits(){
  /******
  std::cout << "Hit has: " << std::endl;
  for(const auto &surface : m_assignedHits){ // returns a pair
      std::cout << "\t" << surface.second.size() << " hits in layer " << surface.first << std::endl;
  }
  ********/
}

std::string HitCollection::hitInfo(){
  std::ostringstream s;
  s << "position: (" << m_hit->X << ", " << m_hit->Y << ", " << m_hit->Z << ")" 
    << "\t surface: " << m_hit->SurfaceID
    << "\t has " << countAssignedHits() << " hits. "; 
  return s.str();
}

void HitCollection::printHit(){
  std::cout << hitInfo() << std::endl;
}




