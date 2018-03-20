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
#include <string>

// testing
#include <ctime>

using namespace std;




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

  m_debug=false;

  if(m_debug) std::cout << "HitToTrack::Init()" << std::endl;

  ////////////////////////////////////////////////////
  // Input parameters (copied from ParticlePropagator) 
  ////////////////////////////////////////////////////
  fInputArray = ImportArray(GetString("InputArray", "Delphes/hit"));
  //fInputArray   = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();



  // Output arrays 
  fTrackOutputArray = ExportArray(GetString("OutputArray", "tracks"));

  ///////////////////////////////
  // parameters for barrel layers
  ///////////////////////////////
  //fBarrelLength = GetDouble("BarrelLength", 1.0); // barrel length [m]
  //ExRootConfParam barrelLayersParam = GetParam("BarrelLayerRadii"); // barrel layer radii [m]
  //Long_t size = barrelLayersParam.GetSize();
  //for(int i = 0; i < size; ++i){
  //fBarrelLayerRadii.push_back( barrelLayersParam[i].GetDouble() ); 
  //}

  ///////////////////////////////
  // parameters for endcap layers
  ///////////////////////////////
  //ExRootConfParam endcapZPositions = GetParam("EndCapZ"); // z positions of the endcaps. Only needed in the positive sense (symmetry assumed) 
  //fEndCapRadius = GetDouble("EndCapRadius", 1.0); // end cap radius [m] 
  //for(int i=0; i<endcapZPositions.GetSize(); ++i){
  //if( endcapZPositions[i].GetDouble() < 0.0){
  //std::cerr << "WARNIGN: HitToTrack: Negative endcap z positions detected, which will be ignored. The user need only input positive values, and a symmetrical detector design is assumed." << std::endl;
  //continue;
  //}
  //fEndcapZPositions.push_back( endcapZPositions[i].GetDouble() ); 
  //}

  if(m_debug) std::cout << "Init(): Defined input and output arrays" << std::endl;

}

//------------------------------------------------------------------------------

void HitToTrack::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------


void HitToTrack::Process()
{

  if(m_debug) std::cout << "Process()" << std::endl;

  Candidate *candidate;

  // loop over the hits in the event
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    std::cout << candidate->HitRadius << std::endl;
    //for(int i=0; i<fNHistograms; ++i) windowHists.at(i)->Fill(candidate->DZ, candidate->PT);
  }
  

  /********************
  int SurfaceID(0);

  // Creat hits for all barrel layers
  //clock_t begin = clock();
  int hitNumber(0);
  for(auto barrelRadius : fBarrelLayerRadii){
    std::vector<TLorentzVector> hits;
    bool removeEndcaps(true);
    bool removeBarrel(false);
    ParticlePropagator(barrelRadius, fBarrelLength, SurfaceID, removeEndcaps, removeBarrel); // all barrels have the same length 
    SurfaceID++;
  }

  //clock_t end = clock();
  //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cout << " Time for call of ParticlePropagator()x3: " << elapsed_secs << std::endl;

  // Create hits for endcap layers
  for(auto endcapZ : fEndcapZPositions){
    std::vector<TLorentzVector> hits;
    bool removeEndcaps(false);
    bool removeBarrel(true);
    ParticlePropagator(fEndCapRadius, endcapZ, SurfaceID, removeEndcaps, removeBarrel); 
    SurfaceID++;
  }
  *************/

}

//------------------------------------------------------------------------------
