/** \class FakeTrackRemover
 *
 *  Associate tracks to a vertex  
 *
 *  \authors W. Fawcett
 *
 */


#include "modules/FakeTrackRemover.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"
#include "classes/LineParameters.h"

#include "TVector3.h"

#include <stdexcept>
#include <iostream>
#include <string>

void printHit(Candidate* hit){
  TLorentzVector pos = hit->Position;
  std::cout << "candidate Position: r: " << pos.Perp() << " x: " << pos.X() << " y: " << pos.Y() << " z: " << pos.Z() << std::endl;
  std::cout << "Xd: " << hit->Xd << " Yd: " << hit->Yd << " Zd: " << hit->Zd << std::endl;
  std::cout << "SurfaceID " << hit->SurfaceID << std::endl;
}


inline float calcDeltaPhiLimit(float r1, float r2){
  // Input arguments:
  // r1: inner radius [mm]
  // r2: outer radius [mm]
  //
  // Calculates bending radius [mm] for track of 2 GeV in 4 T magnetic field 
  //  pT[GeV] = 0.29976 * B[T] * R[m] 
  //  ==> R[m] = 3.3366 * pT[GeV]/B[T]
  float twiceRadius = 1000 * 3.3366 * 2 / 4.0; // bending radius  
  twiceRadius *= 2; // twice the radius  
  return fabs(acos(r1 / twiceRadius) - acos(r2 / twiceRadius));
}

inline TVector3 fillTV3(float r, float phi, float z){

  // calculate cartesian coordinates
  //float theta = 2*atan( exp(-eta) );
  //float z = r*sin(theta); 
  float x = r*cos(phi);
  float y = r*sin(phi);

  TVector3 tempVec(x, y, z);
  return tempVec; 
}

//------------------------------------------------------------------------------

FakeTrackRemover::FakeTrackRemover() 
{
}

//------------------------------------------------------------------------------

FakeTrackRemover::~FakeTrackRemover()
{
}

//------------------------------------------------------------------------------

void FakeTrackRemover::Init()
{

  // debug 
  m_debug = GetBool("debug", false); 

  // calculate deltaphi limits (note, input units are [mm])
  m_deltaPhi12Limit = calcDeltaPhiLimit(552.0, 582.0);
  m_deltaPhi23Limit = calcDeltaPhiLimit(582.0, 612.0);
  m_deltaPhi13Limit = calcDeltaPhiLimit(552.0, 612.0); 

  // "should be" 
  //12: 0.00912392
  //23: 0.0091387
  //13: 0.0182626

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "HitToTrack/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // Output arrays 
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void FakeTrackRemover::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
}

//------------------------------------------------------------------------------

void FakeTrackRemover::Process()
{
  ////////////////////////////////////////////////
  // loops over all input tracks
  // Applies selection cuts to attempt to remove fake tracks
  ////////////////////////////////////////////////
  
  // debug 
  int nTracksOriginal(0);
  int nTracksOriginalPt2(0);
  int nTracksSurviving(0);
  int nTracksSurvivingPt2(0);
  int nFakeOriginal(0);
  int nFakeSurviving(0);
  
  // Match tracks to vertex 
  Candidate *track;
  fItTrackInputArray->Reset();
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {

    float deltaKappa = track->kappa_013 - track->kappa_123;
    float pT = track->PT;

    // debug
    if(m_debug){
      nTracksOriginal++;
      if(track->PT > 2) nTracksOriginalPt2++;
      if(track->IsFake) nFakeOriginal++;
    }
    
    // some deltaKappa cuts
    if(pT > 0 && fabs(deltaKappa)    > 0.003) continue;
    if(pT > 2.0 && fabs(deltaKappa)  > 0.002) continue; // 0.002 is super efficient!
    if(pT > 2.0 && fabs(deltaKappa)  > 0.0015) continue; // new, not harsh  
    if(pT > 3.5 && fabs(deltaKappa)  > 0.001) continue; // new 
    if(pT > 50.0 && fabs(deltaKappa) > 0.0005) continue;

    // get hits
    auto hits = track->GetCandidates();
    //std::cout << "This track has " << hits->GetEntries() << " candidates" << std::endl; 
    Candidate* hit1 = static_cast<Candidate*>(hits->At(1));
    Candidate* hit2 = static_cast<Candidate*>(hits->At(2));
    Candidate* hit3 = static_cast<Candidate*>(hits->At(3));

    // Delta-phi constraints
    float deltaPhi12 = fabs(hit1->Phi - hit2->Phi);
    float deltaPhi23 = fabs(hit2->Phi - hit3->Phi);
    float deltaPhi13 = fabs(hit1->Phi - hit3->Phi);

    // WJF: remove delta phi ... somehow failing? Maybe problem with Phi calculation?? 
    if(deltaPhi12 > m_deltaPhi12Limit) continue;
    if(deltaPhi23 > m_deltaPhi23Limit) continue;
    if(deltaPhi13 > m_deltaPhi13Limit) continue;


    // recalcualte zresiduum (could add this calculation to HitToTrack::CalculateTrackParametersTriplet ), and add it to the track class?  
    LineParameters params;
    params.calculateLineParameters(hit1->Position.Z(), hit1->Position.Perp(), hit3->Position.Z(), hit3->Position.Perp());
    float intersect = (582 - params.y_intercept())/params.gradient();
    float zresiduum = hit2->Position.Z() - intersect; // should not be larger than 0.5
    float beamlineIntersect = params.x_intercept();

    // Calculate sign of angular differences between points (in 3D)
    TVector3 p1 = fillTV3(hit1->Position.Perp(), hit1->Phi, 0); 
    TVector3 p2 = fillTV3(hit2->Position.Perp(), hit2->Phi, 0); 
    TVector3 p3 = fillTV3(hit3->Position.Perp(), hit3->Phi, 0); 
    TVector3 zhat(0, 0, 1); // z direction vector

    // z component of the cross product of two angles
    float z_phi12 = p1.Cross(p2)*zhat / (p1.Mag() * p2.Mag()); 
    float z_phi13 = p1.Cross(p3)*zhat / (p1.Mag() * p3.Mag());
    float z_phi23 = p2.Cross(p3)*zhat / (p2.Mag() * p3.Mag());
    float z_12m23 = z_phi12 * z_phi23;
    float z_12p23 = z_phi12 + z_phi13; 

    // z residuum 
    const float zresiduumCut = 0.1; 
    if(fabs(zresiduum) > zresiduumCut) continue; 

    // z residuum over eta (extra information)
    if(fabs( zresiduum/fabs(track->Eta) ) > 0.4) continue; 

    // Phi constraint, such that the differences in phi angles are in the same sense when moving from one tracking layer to the next
    if(z_12m23 < 0) continue;

    // beamline intersect
    if(fabs(beamlineIntersect) > 200) continue; 

    if(m_debug){
      nTracksSurviving++;
      if(track->PT > 2) nTracksSurvivingPt2++;
      if(track->IsFake) nFakeSurviving++;
    }

    fOutputArray->Add(track);

  }

  if(m_debug){
    std::cout << "FakeTrackRemover::Process(): nTracksOriginal: " << nTracksOriginal << std::endl;
    std::cout << "FakeTrackRemover::Process(): nTracksOriginalPt2: " << nTracksOriginalPt2 << std::endl;
    std::cout << "FakeTrackRemover::Process(): nTracksSurviving: " << nTracksSurviving << std::endl;
    std::cout << "FakeTrackRemover::Process(): nTracksSurvivingPt2: " << nTracksSurvivingPt2 << std::endl;
    std::cout << "FakeTrackRemover::Process(): nFakeOriginal: " << nFakeOriginal << std::endl;
    std::cout << "FakeTrackRemover::Process(): nFakeSurviving: " << nFakeSurviving << std::endl;
  }



}

//------------------------------------------------------------------------------
