#include "classes/anaClasses.h" 
#include <cmath>
#include <utility>
#include "TMath.h"

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
  //if(beta1 <= 0.0) phi1 += 2*M_PI;

  float phi2 = atan2(beta2, alpha1);
  //if(beta2 <= 0.0) phi2 += 2*M_PI; 

  //std::cout << "a, b, r : " << a << " " << b << " " << r  << std::endl;
  //std::cout << "(a1, b2) (" << alpha1 << ", " << beta2 << ") : phi = " << phi2/M_PI << "pi" << std::endl;
  //std::cout << "(a2, b1) (" << alpha2 << ", " << beta1 << ") : phi = " << phi1/M_PI << "pi" <<  std::endl;
  //std::cout << "" << std::endl; 
  
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

  // get all barrel layers
  std::vector<int> layers;
  for(auto const& key : hitMap){
    layers.push_back(key.first);
  }

  /////////////////////
  // Fill HitCollection
  /////////////////////
  
  // first reserve some space (performance)
  int numHits(0);
  for(auto layer : layers){
    numHits += hitMap[layer].size();
  }
  m_associatedHitCollection.reserve(numHits+1);


  // Fill hits 
  for(auto layer : layers){
    for(Hit* hit : hitMap[layer]){
      m_associatedHitCollection.push_back(HitCollection(hit));
    }
  }

  // associate hit algorithm 
  // It's very important to loop over the references to the object
  // if a copy is made, then the pointer storing below wouldn't work
  for(HitCollection& innerHit : m_associatedHitCollection){

    // first loop over inner hits 
    if(innerHit.SurfaceID != layers.at(0)) continue; 

    // parameters of for the search window
    float r = innerHit.Perp();
    float z = innerHit.Z();

    // second loop over outer hits 
    for(HitCollection& outerHit : m_associatedHitCollection){
      if(outerHit.SurfaceID == layers.at(0)) continue; // skip first layer (only want seeds from first layer)

      // calculate search window
      float rOuter = outerHit.Perp();
      float zRight  = calculateZWindowForNextLevel(r, z, rOuter, minZ); 
      float zLeft = calculateZWindowForNextLevel(r, z, rOuter, maxZ); 

      // doesn't really need to be calculated for each hit ... 
      float maxPhiDeviation = calculateRPhiWindowInToOut(innerHit.X(), innerHit.Y(), rOuter);

      std::cout << "1->2 " << calculateRPhiWindowInToOut(532, 0, 582) << std::endl;
      std::cout << "1->3 " << calculateRPhiWindowInToOut(532, 0, 632) << std::endl;
      std::cout << "2->3 " << calculateRPhiWindowInToOut(582, 0, 632) << std::endl;

      float phiDeviation = innerHit.DeltaPhi(outerHit);

      // determine if matched
      if( (zLeft < outerHit.Z() && outerHit.Z() < zRight) && (phiDeviation < maxPhiDeviation) ){
        std::cout << "add hit " << phiDeviation << std::endl;
        innerHit.addHit(&outerHit); // we're storing pointers to object stored in anoter vector
        // usually it's a bad idea to store pointers to objects defined on the stack
        // but m_associatedHitCollection won't be modified until it's deleted 
      }
    } // loop over outer hits
  } // loop over inner hits


  for(auto hitCollection : m_associatedHitCollection){
    if(hitCollection.SurfaceID == 0){
      hitCollection.printMatchedHits();
      //hitCollection.printAssignedHitPointers();
    }
  }

  return true;
}

// create tracks from hit collections 
bool TrackFitter::combineHitsToTracksInToOut(){

  for(auto hitCollection : m_associatedHitCollection){

    // consider innermost hits 
    if(hitCollection.SurfaceID == 0){
      //hitCollection.printHit();
      hitCollection.printMatchedHits();

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


float HitCollection::DeltaPhi(HitCollection& h) const  {
  
  return acos( cos ( this->Phi() - h.Phi() ));
      //return m_position.DeltaPhi(h.GetPosition()); // used TLorentzVector DeltaPhi function
}

void HitCollection::printMatchedHits(int level) const{
  std::cout << this->hitInfo() << std::endl; 
  int iHit(1);

  // loop over pointers of subhits  
  for(HitCollection* subHit : m_assignedHits){
    for(int i=0; i<level+1; ++i) std::cout << "\t";
    std::cout << "h" << iHit << " " << subHit->hitInfo() <<  "deltaPhi " << this->DeltaPhi(*subHit) << std::endl;
    iHit++;

    // recursively print out hits belonging to any subhits
    if(subHit->countAssignedHits() > 0){
      subHit->printMatchedHits(level+1);
    }

  }
}

void HitCollection::printMatchedHits() const{
  this->printMatchedHits(0);
}


std::string HitCollection::hitInfo() const{
  std::ostringstream s;
  s << "position: (" << m_hit->X << ", " << m_hit->Y << ", " << m_hit->Z << ")"  // cartesian
    << " (" << m_hit->Perp() << "," << m_hit->Phi() << ")" // (r, phi, z) cylindrical polars (z already printed, so not shown here)
    << "\t surface: " << m_hit->SurfaceID
    << "\t has " << this->countAssignedHits() << " hits. "; 
  return s.str();
}

void HitCollection::printHit() const{
  std::cout << this->hitInfo() << std::endl;
}
