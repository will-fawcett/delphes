#include "classes/myTrack.h"
#include <iostream>
#include "classes/UtilityFunctions.h"

bool myTrack::calculateTrackParameters( cartesianCoordinate coord, trackParamAlgo algorithm ){
  // function to calculate d0, z0 (and in principle theta, phi, qOverP)
  // relative to the specified coordinate!!!!!
  // TODO: add calculation if with respect to NOT (0, 0, 0)
  
  switch (algorithm) {
    case beamlineConstraint:{
      return this->trackParametersBeamlineConstraint();
    }
    case noBeamlineConstraint:{
      return this->trackParametersNoBeamlineConstraint();
    }
    case MAXIMUM:{
      return false;
    }
  }

  return true; 

}


bool myTrack::trackParametersNoBeamlineConstraint(){
  // calculates the track parameters only using the three points in the triplet

  std::cout << "trackParametersNoBeamlineConstraint()" << std::endl;

  if(m_associatedHits.size() != 3){
    std::cerr << "ERROR: more than three hits associate to this track. Algorith is not compatible." << std::endl;
    return false;
  }
  Hit * hit1 = m_associatedHits.at(0);
  Hit * hit2 = m_associatedHits.at(1);
  Hit * hit3 = m_associatedHits.at(2);
  // Check hits are properly orderd radially. Shouldn't really happen
  if(hit1->Perp() > hit2->Perp()){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    return false;
  }
  if(hit2->Perp() > hit3->Perp()){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    return false;
  }
  
  // calculate the center of the circle described by the points (xi, yi) for i = 1, 2, 3

  // copy of hit coordinates
  float x1 = hit1->X;
  float y1 = hit1->Y; 
  float z1 = hit1->Z; 

  float x2 = hit2->X;
  float y2 = hit2->Y; 
  float z2 = hit2->Z; 

  float x3 = hit3->X;
  float y3 = hit3->Y; 
  float z3 = hit3->Z; 

  float r1 = hit1->Perp();
  float r2 = hit2->Perp();
  float r3 = hit3->Perp();

  // circle described by (x-a)^2 + (y-b)^2 = r^2 
  float b = ( (x3 - x1)*(r2*r2 - r1*r1) - (x2 - x1)*(r3*r3 - r1*r1) ) / ( 2*( (y2-y1)*(x3-x1) - (y3-y1)*(x2-x1)  ) );

  float a = ( r2*r2 - r1*r1 - 2*b*(y2 - y1) ) / ( 2*(x2 - x1));

  float radius = hypotf( (x1-a), (y1-b) );

  std::cout << "(xi, yi) : "
    << x1 << ", " << y1 << "\t"
    << x2 << ", " << y2 << "\t"
    << x3 << ", " << y3 << "\t"
    << std::endl;

  std::cout << "(a, b) and radius : " << a << ", " << b << "\t" << radius << std::endl;

  // transverse momentum
  m_pT = 1.199 * fabs(radius/1000); // for a 4T magnetic field and radius in [m], divide by 1000 as length units in [mm]

  // don't know yet how to calculate d0 (and phi) 
  m_phi = 0;
  m_d0 = 0;
  
  // longitudinal parameters calculated from least-squares fit of the three hit points
  lineParameters params = simpleLinearLeastSquaresFit(m_associatedHits);
  m_z0 = params.x_intercept; 

  m_theta = 0; // careful how this is calculated, may only want this to be defined from [0; pi] 


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
  if(hit1->Perp() > hit2->Perp()){
    std::cerr << "ERROR: hits not in the correct order!" << std::endl;
    return false;
  }
  if(hit2->Perp() > hit3->Perp()){
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

  // Calculate the parameters in the longitudinal plane
  // Follows calculations by A. Schoning

  float r01 = hypotf(x1, y1);
  float r03 = hypotf(x3, y3); 
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
  //float radiusAndre = ( (x1*x1 - x3*x3) + (y1*y1 - y3*y3) ) / ( 2*(x1 - x3) ); // not sure if formulea in Andres paper is quite correct ... ? 
  float radius = sqrt(a*a + b*b); 

  // kappa (1/radius) 
  float kappa = 1/radius; 

  // transverse momentum
  float pT = 1.199 * fabs(radius/1000); // for a 4T magnetic field and radius in [m], divide by 1000 as length units in [mm]

  // arc lengths 
  float s1 = radius * PHI1;
  float s3 = radius * PHI3; 

  // phi angle given by line tangent to the circle at (0,0) 
  float phi = atan2(-a, b);

  // Assign the track parameters
  m_z0 = z1 - s1* ( z3 - z1 ) / (s3 - s1);
  m_theta = atan2( (s3-s1), (z3-z1) );  // CHECK this is correct, may only want from [0:pi]
  m_d0 = 0.0; // by definition (beamline constraint) 
  m_pT = pT;
  m_phi = phi;
  return true;

}


bool myTrack::isFake() const{
  return (!this->isNotFake()); 
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
