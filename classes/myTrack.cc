#include "classes/myTrack.h"

bool myTrack::calculateTrackParameters( cartesianCoordinate coord ){
  // function to calculate d0, z0 (and in principle theta, phi, qOverP)
  // relative to the specified coordinate!!!!!
  // TODO: add calculation if with respect to NOT (0, 0, 0)
  
  switch(trackParamAlgo){
    case beamlineConstraint:{
      return this->trackParametersBeamlineConstraint()
    }
    case noBeamlineConstraint:{
      return false;
    }
  }

  // parameters relative to (0, 0, 0)
  phi = atan2( (m_gradient + m_intercept), 1);
  z0 = -1 * (m_intercept / m_gradient);

  d0=0; // can't calculate this with straight line
  theta=0; // not calculable with 2D line
  qOverP=0; // straight line : infinite radius of curvature (therefore infinite momentum)  

}




bool myTrack::trackParametersBeamlineConstraint(){
  // assumes the track originates from (0, 0, z0)
  // Calculates the track parameters using the origin and two other points, corresponding to the innermost and outermost hit
  
  if(m_associatedHits.size() != 3){
    std::cerr << "ERROR: more than three hits associate to this track. Algorith is not compatible." << std::endl;
    return false;
  }

  Hit * hit1 = m_associatedHits.at(0);
  Hit * hit2 = m_associatedHits.at(1);
  Hit * hit3 = m_associatedHits.at(2);

  // Check hits are properly orderd radially. Shouldn't really happen
  if(hit1.Perp() < hit2.Perp()){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    return false;
  }
  if(hit2.Perp() < hit3.Perp(){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    return false;
  }

  
  // copy of hit coordinates
  float x1 = hit1->X;
  float y1 = hit1->Y; 
  float z1 = hit1->Z; 

  float x3 = hit3->X;
  float y3 = hit3->Y; 
  float z3 = hit3->Z; 

  // Calculate the track parameters given the two points and the beamline constraint
  // This can be solved analytically ... 
   
  // radius of curvature 
  float radius = ( (x1*x1 - x3*x3) + (y1*y1 - y3*y3) ) / ( 2*(x1 - x3) );

  // kappa (1/radius) 
  float kappa = 1/radius; 

  // transverse momentum
  float pT = 1.199 * fabs(radius); // for a 4T magnetic field 

  // Now calculate the parameters in the longitudinal plane
  // Follows calculations by A. Schoning

  float r01 = hypotf(x1, y1);
  float r03 = hypotf(x3, y3); 
  float r13 = hypotf( fabs(x3-x1) , fabs(y3-y1) ); // not sure if fabs actuall needed ...  
  float cord13 = x1*y3 - y1*x3; 
  float PHI1 = 2 * arcsin( chord13 / (r13*r03) );
  float PHI3 = 2 * arcsin( chord13 / (r13*r01) );
  float s1 = radius * PHI1;
  float s2 = radius * PHI2; 

  float z0 = z1 - s1* ( z3 - z1 ) / (s3 - s1);

  float theta = atan2( hit1.Y, hit1.Perp() );

  // set the track parameters for this track
  parameters.pT = pT; 
  parameters.d0 = 0.0; // by definition 
  parameters.z0 = z0;
  parameters.theta = theta; 



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
