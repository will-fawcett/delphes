#include "classes/myTrack.h"

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

bool myTrack::isNotFake() const{

  // If all of the hits associated to this track have the same particle ID, the the track is not a fake track
  //
  // If this hit is from pileup, the pileup particle is unlikely to be stored in the ROOT file, and so the reference to it will break
  // In this case we cannot definitely tell if the track is fake or not, however we can match by pT. If all three particles have (exactly) the same pT
  // then they should have originated from the same particle 

  bool puFlag(false);
  std::vector<int> uniqueIDs;
  std::vector<int> ptIDs;
  for(Hit* hit : m_associatedHits){
    if(hit->IsPU){
      puFlag=true;
      int intPtKeVID = hit->intPtKeVID; 
      ptIDs.push_back(intPtKeVID);
    }
    else{
      // only fill if the UID exists (i.e. is not a PU particle) 
      auto uniqueID = dynamic_cast<GenParticle*>(hit->Particle.GetObject())->GetUniqueID();
      uniqueIDs.push_back(uniqueID);
    }
  }
  
  if(puFlag){
    // check if the pT IDs are the same
    for(int i=0; i<ptIDs.size()-1; ++i){
      if(ptIDs.at(i) != ptIDs.at(i+1)) return false;
    }
  }
  else{
    // check if the unique IDs are the same
    for(int i=0; i<uniqueIDs.size()-1; ++i){
      if(uniqueIDs.at(i) != uniqueIDs.at(i+1) ) return false;
    }
  }
  return true;
}
