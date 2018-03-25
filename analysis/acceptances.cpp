// Delphes classes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

// plotting stuff (shared)
//#include "plotting.h"

// c++ libs
#include <iostream>
#include <sstream>
#include <exception>
#include <map>
#include <ctime>

// stuff for ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TNamed.h"





// Rule for sorting vector of Track* objects 
bool reverse(const Track* i, const Track* j){
  return i->PT > j->PT; 
}

//------------------------------------------------------------------------------

void printTime(clock_t begin, clock_t end, std::string message){
  std::cout << message << ": time elapsed: " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
}

//------------------------------------------------------------------------------

struct Plots {
  std::vector<TH1*> track1Pt; 
  std::vector<TH1*> track2Pt; 
  std::vector<TH1*> track3Pt; 
  std::vector<TH1*> track4Pt; 

};


void BookHistograms(ExRootResult *result, Plots *plots, std::vector<std::string> branchNames)
{

  // one of each type of histogram for each branchname 
  for(int i=0; i<branchNames.size(); ++i){

    std::string branch = branchNames.at(i);

    plots->track1Pt.push_back(
        result->AddHist1D( branch+"track1Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track2Pt.push_back(
        result->AddHist1D( branch+"track2Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track3Pt.push_back(
        result->AddHist1D( branch+"track3Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track4Pt.push_back(
        result->AddHist1D( branch+"track4Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );

  }
}

//------------------------------------------------------------------------------

void AnalyseEvents(const int nEvents, ExRootTreeReader *treeReader, Plots *plots, std::vector<std::string> branchNames)
{


  std::vector<TClonesArray*> branches;
  for(const auto& bName : branchNames){
    branches.push_back( treeReader->UseBranch(bName) );
  }

  // Define branches
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  //TClonesArray *branchTruthTrack = treeReader->UseBranch("TruthTrack");
  //TClonesArray *branchTrack      = treeReader->UseBranch("Track");

  //TClonesArray* branchTracksFromHit30 = treeReader->UseBranch("TracksFromHit30");
  //TClonesArray* branchPBMatchedHitTracks = treeReader->UseBranch("PBMatchedHitTracks");
  //TClonesArray* branchSmearedTracksFromHits = treeReader->UseBranch("SmearedTracksFromHits");

  



  // Loop over all events
  Long64_t allEntries = treeReader->GetEntries();
  std::cout << "** Chain contains " << allEntries << " events" << std::endl;
  float eventWeight = 1.0/allEntries; // an event weight for normalising histograms
  for(Long64_t entry = 0; entry < allEntries; ++entry)
  {

    // limit number of events looped over
    if(nEvents != -1){
      if(entry > nEvents-1) break; // -1 because humans will count from 1 event
    }

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    // print every 10% complete
    if( entry % int(allEntries/10)==0 ) std::cout << "Event " << entry << " out of " << allEntries << std::endl;
    if( entry % 100==0 ) std::cout << "Event " << entry << " out of " << allEntries << std::endl;

    /**************
    int pCounter(0);
    std::vector<float> particlePts; 
    particlePts.reserve(1000);
    for(auto itParticle=branchParticle->begin(); itParticle != branchParticle->end(); ++itParticle){
      GenParticle* particle = dynamic_cast<GenParticle*>(*itParticle);
      if(particle->Status == 1){
        std::cout << "p" << pCounter << " " << particle->PT << std::endl;
        particlePts.push_back(particle->PT);
        pCounter++;
      }
    }
    std::cout << "" << std::endl;

    std::sort(particlePts.rbegin(), particlePts.rend());
    for(auto pt : particlePts){
      std::cout << pt << std::endl;
    }
    std::cout << "" << std::endl;
    continue;
    ***********/


    for(int iBranch=0; iBranch<branchNames.size(); ++iBranch){
      TClonesArray* branch = branches.at(iBranch);

      //std::cout << "branch: " << branchNames.at(iBranch) << std::endl;

      // Store the tracks, and sort into descending pT order
      int nTracks = branch->GetEntriesFast(); 
      std::vector<Track*> sortedTracks; 
      sortedTracks.reserve(nTracks);
      for(int iTrack=0; iTrack<nTracks; ++iTrack){
        Track* track = (Track*) branch->At(iTrack);
        sortedTracks.push_back(track);
      }
      std::sort( sortedTracks.begin(), sortedTracks.end(), reverse);

      // Loop over the sorted tracks 
      for(int iTrack=0; iTrack<nTracks; ++iTrack){
        Track* track = sortedTracks.at(iTrack);
        if(iTrack == 0) plots->track1Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 1) plots->track2Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 2) plots->track3Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 3) plots->track4Pt.at(iBranch)->Fill(track->PT, eventWeight);
      }


        //GenParticle* particle = (GenParticle*) track->Particle.GetObject();
        //std::cout << "track pt: " << track->PT << "\tparticle pt: " << particle->PT << "\t resolution: " << track->PT - particle->PT << std::endl;
        /********
        // pt resolution
        auto can = (Candidate*) branch->At(iTrack);
        auto collection = can->GetCandidates();
        std::cout << "branch: " << branchNames.at(iBranch) << " nCandidates: " << collection->GetEntriesFast() << std::endl; 
        auto hit0 = collection->At(0); // particle that made it (?)
        auto hit1 = collection->At(1); // hit1
        auto hit2 = collection->At(2); // hit2
        auto hit3 = collection->At(3); // hit3 
        **********/

      /********
      int iTrack(0);
      for(auto itTrack=branch->begin(); itTrack != branch->end(); ++itTrack){
        Track* track = dynamic_cast<Track*>(*itTrack);

        plots->track1Pt.at(i)->Fill();
      }
      **********/

    }
    
        //for(auto itTrack=branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
        //}


    



  } // end loop over entries

} // end AnalyseEvents



//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, Plots *plots)
{
  result->Print("pdf");
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  gROOT->SetBatch(1);
  //gROOT->ProcessLine("#include <vector>");

  std::string appName = "acceptances";
  std::cout << "Will execute " << appName << std::endl;


  // Get input agruments
  std::string inputFile = argv[1]; 
  int nEvents = atoi(argv[2]);
  nEvents = atoi(argv[2]);


  // info relating to input arguments
  if(nEvents < 0){
    std::cout << "INFO: Will run over all entries." << std::endl;
  }
  else{
    std::cout << "INFO: Will run over " << nEvents << " entries. " <<  std::endl;
  }


  // branches to loop over
  std::vector<std::string> branchNames = {
        // default delphes tracks
        "TruthTrack",
        "Track",
        "PBMatchedTracks",

        // tracks from hits 
        "TracksFromHit30", 
        "SmearedTracksFromHits",
        "PBMatchedHitTracks", 
  };

  // control analysis
  TChain *chain = new TChain("Delphes");
  std::cout << "Adding " << inputFile << " to chain" << std::endl;
  chain->Add(inputFile.c_str());


  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  Plots *plots = new Plots;
  BookHistograms(result, plots, branchNames);
  AnalyseEvents(nEvents, treeReader, plots, branchNames);

  // Turn off, never really want this ... 
  int doPrintHistograms = 0; 
  if(doPrintHistograms) PrintHistograms(result, plots);

  std::string outputFile = "test.root";
  std::cout << "Writing to file: " << outputFile << std::endl;
  result->Write(outputFile.c_str());

  std::cout << "** Exiting..." << std::endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

  return 0;
}
