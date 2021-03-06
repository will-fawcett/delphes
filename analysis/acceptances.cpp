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
#include <sys/stat.h>
#include <getopt.h> // argument parsing! 
#include <glob.h> 

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

// global N_JETS
int N_JETS = 6; 
bool m_debug = false;

//------------------------------------------------------------------------------

// Fast way to test if a file exists
inline bool fileExists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

//------------------------------------------------------------------------------

// Rule for sorting vector of Track* objects 
bool reverseTrack(const Track* i, const Track* j){
  return i->PT > j->PT; 
}

// Rule for sorting vector of Jet* objects 
bool reverseJet(const Jet* i, const Jet* j){
  return i->PT > j->PT; 
}

// Rule for sorting vector of Candidate* objects
bool reverseCandidate(const Candidate* i, const Candidate* j){
  return i->PT > j->PT; 
}

//------------------------------------------------------------------------------
// a simple glob https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
inline std::vector<std::string> glob(const std::string& pat){
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    std::vector<std::string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(std::string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

//------------------------------------------------------------------------------

void printTime(clock_t begin, clock_t end, std::string message){
  std::cout << message << ": time elapsed: " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
}

//------------------------------------------------------------------------------

struct Plots {

  TH1* numEntries;
  std::vector<TH1*> fakeTrackPt;
  std::vector<TH1*> trueTrackPt;

  std::vector<TH1*> allTrackPt;
  std::vector<TH1*> allTrackEta;
  std::vector<TH1*> track1Eta;
  std::vector<TH2*> allTrackPt_Eta;
  std::vector<TH2*> track1Pt_Eta;

  std::vector<TH1*> electron1Pt;
  std::vector<TH1*> electron2Pt;
  std::vector<TH1*> lepton1Pt;
  std::vector<TH1*> lepton2Pt;
  std::vector<TH1*> muon1Pt;
  std::vector<TH1*> muon2Pt;

  std::vector<TH1*> track1Pt; 
  std::vector<TH1*> track2Pt; 
  std::vector<TH1*> track3Pt; 
  std::vector<TH1*> track4Pt; 


  std::vector<TH1*> ptRes;
  std::vector<TH2*> ptRes_pt;
  std::vector<TH2*> ptRes_eta;
  std::vector<TH3*> ptRes_pt_eta;

  std::vector<TH1*> ptResRaw;
  std::vector<TH2*> ptResRaw_pt;
  std::vector<TH2*> ptResRaw_eta;
  std::vector<TH3*> ptResRaw_pt_eta;

  std::vector<TH1*> z0Res;
  std::vector<TH2*> z0Res_pt;
  std::vector<TH2*> z0Res_eta;
  std::vector<TH3*> z0Res_pt_eta;

  std::vector<TH1*> d0Res;
  std::vector<TH2*> d0Res_pt;
  std::vector<TH2*> d0Res_eta;
  std::vector<TH3*> d0Res_pt_eta;

  std::vector<TH1*> phiRes;
  std::vector<TH2*> phiRes_pt;
  std::vector<TH2*> phiRes_eta;
  std::vector<TH3*> phiRes_pt_eta;

  //std::vector<TH1*> etaRes;
  //std::vector<TH2*> etaRes_pt;
  //std::vector<TH2*> etaRes_eta;

  std::vector<std::vector<TH1*>> jetiPt;

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, Plots *plots, std::vector<std::string> trackBranchNames, std::vector<std::string> jetBranchNames, bool doResolutionPlots)
{

  // book keeping
  plots->numEntries = result->AddHist1D("numEntries", "", "", "", 2, 0, 2, 0, 0);

  // one of each histogram for each branchname 
  for(int i=0; i<trackBranchNames.size(); ++i){

    std::string branch = trackBranchNames.at(i);

    // Fake and true tracks
    plots->fakeTrackPt.push_back(
        result->AddHist1D(branch+"_fakeTrackPt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->trueTrackPt.push_back(
        result->AddHist1D(branch+"_trueTrackPt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->allTrackPt.push_back(
        result->AddHist1D(branch+"_allTrackPt", "", "", "", 1000, 0, 1000, 0, 0)
        );

    // Track eta / pt 
    plots->allTrackEta.push_back(
        result->AddHist1D(branch+"_allTrackEta", "", "", "", 200, -10, 10, 0, 0)
        );
    plots->track1Eta.push_back(
        result->AddHist1D(branch+"_track1Eta", "", "", "", 200, -10, 10, 0, 0)
        );

    plots->allTrackPt_Eta.push_back(
        result->AddHist2D((branch+"_allTrackPt_Eta").c_str(), "", "", "", 200, 0, 200, 200, -10, 10)
        );
    plots->track1Pt_Eta.push_back(
        result->AddHist2D((branch+"_track1Pt_eta").c_str(), "" ,"", "", 200, 0, 200, 200, -10, 10)
        );

    // Nth track pT
    plots->track1Pt.push_back(
        result->AddHist1D( branch+"_track1Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track2Pt.push_back(
        result->AddHist1D( branch+"_track2Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track3Pt.push_back(
        result->AddHist1D( branch+"_track3Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );
    plots->track4Pt.push_back(
        result->AddHist1D( branch+"_track4Pt", "", "", "", 1000, 0, 1000, 0, 0)  
        );

    // Lepton pT
    plots->electron1Pt.push_back(
        result->AddHist1D(branch+"_electron1Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->electron2Pt.push_back(
        result->AddHist1D(branch+"_electron2Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->muon1Pt.push_back(
        result->AddHist1D(branch+"_muon1Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->muon2Pt.push_back(
        result->AddHist1D(branch+"_muon2Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->lepton1Pt.push_back(
        result->AddHist1D(branch+"_lepton1Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );
    plots->lepton2Pt.push_back(
        result->AddHist1D(branch+"_lepton2Pt", "", "", "", 1000, 0, 1000, 0, 0)
        );

    if(doResolutionPlots){
      // resolution histograms
      float ptMin(-10), ptMax(10);
      plots->ptRes.push_back(
          result->AddHist1D(branch+"_ptRes", "pt resolution", "#deltap_{T}/p_{T}", "Number of Tracks", 1000, ptMin, ptMax, 0, 0)
          );
      plots->ptRes_pt.push_back(
          result->AddHist2D(branch+"_ptRes_pt", "pt resolution", "#deltap_{T}/p_{T}", "Truth p_{T} [GeV]", 1000, ptMin, ptMax, 200, 0, 1000)
          );
      plots->ptRes_eta.push_back(
          result->AddHist2D(branch+"_ptRes_eta", "eta resolution", "#deltap_{T}/p_{T}", "Truth #eta", 1000, ptMin, ptMax, 200, -2, 2)
          );
      plots->ptRes_pt_eta.push_back(
          result->AddHist3D(branch+"_ptRes_pt_eta", "pt resolution", "#deltap_{T}/p_{T}", "Truth p_{T} [GeV]", "Truth #eta", 1000, ptMin, ptMax, 200, 0, 200, 100, -2, 2)
          );

      // wont be able to calculate pT resolutions for tracks with > 200 GeV
      float ptMinRaw(-200), ptMaxRaw(200); 
      plots->ptResRaw.push_back(
          result->AddHist1D(branch+"_ptResRaw", "pt resolution", "#deltap_{T} [GeV]", "Number of Tracks", 1000, ptMinRaw, ptMaxRaw, 0, 0)
          );
      plots->ptResRaw_pt.push_back(
          result->AddHist2D(branch+"_ptResRaw_pt", "pt resolution", "#deltap_{T} [GeV]", "Truth p_{T} [GeV]", 1000, ptMinRaw, ptMaxRaw, 200, 0, 2000)
          );
      plots->ptResRaw_eta.push_back(
          result->AddHist2D(branch+"_ptResRaw_eta", "pt resolution", "#deltap_{T} [GeV]", "Truth #eta", 1000, ptMinRaw, ptMaxRaw, 200, -2, 2)
          );
      plots->ptResRaw_pt_eta.push_back(
          result->AddHist3D(branch+"_ptResRaw_pt_eta", "pt resolution", "#deltap_{T} [GeV]", "Truth p_{T} [GeV]", "Truth #eta", 4000, ptMinRaw, ptMaxRaw, 200, 0, 200, 40, -2, 2)
          ); // lots of bins needed since pT resolution degreades over many decades

      // z0
      float z0Min(-20), z0Max(20); // large range for low pT tracks needed
      plots->z0Res.push_back(
          result->AddHist1D(branch+"_z0Res", "z0 resolution", "#deltaz_{0} [mm]", "Number of Tracks", 1000, z0Min, z0Max)
          );
      plots->z0Res_pt.push_back(
          result->AddHist2D(branch+"_z0Res_pt", "z0 resolution", "#deltaz_{0} [mm]", "Truth p_{T} [GeV]", 1000, z0Min, z0Max, 200, 0, 2000)
          );
      plots->z0Res_eta.push_back(
          result->AddHist2D(branch+"_z0Res_eta", "z0 resolution", "#deltaz_{0} [mm]", "Truth #eta", 1000, z0Min, z0Max, 200, -2, 2)
          );
      plots->z0Res_pt_eta.push_back(
          result->AddHist3D(branch+"_z0Res_pt_eta", "z0 resolution", "#deltaz_{0} [mm]", "Truth p_{T} [GeV]", "Truth #eta", 2000, z0Min, z0Max, 200, 0, 200, 100, -2, 2)
          ); // particularly tricky, lots of bins

      // d0
      float d0Min(-50), d0Max(50); // large range for low pT tracks needed
      plots->d0Res.push_back(
          result->AddHist1D(branch+"_d0Res", "d0 resolution", "#deltad_{0} [mm]", "Number of Tracks", 1000, d0Min, d0Max)
          );
      plots->d0Res_pt.push_back(
          result->AddHist2D(branch+"_d0Res_pt", "d0 resolution", "#deltad_{0} [mm]", "Truth p_{T} [GeV]", 1000, d0Min, d0Max, 200, 0, 2000)
          );
      plots->d0Res_eta.push_back(
          result->AddHist2D(branch+"_d0Res_eta", "d0 resolution", "#deltad_{0} [mm]", "Truth #eta", 1000, d0Min, d0Max, 200, -2, 2)
          );
      plots->d0Res_pt_eta.push_back(
          result->AddHist3D(branch+"_d0Res_pt_eta", "d0 resolution", "#deltad_{0} [mm]", "Truth p_{T} [GeV]", "Truth #eta", 2000, d0Min, d0Max, 200, 0, 200, 100, -2, 2)
          ); // particularly tricky, lots of bins

    float phiMin(-0.1), phiMax(0.1);
    plots->phiRes.push_back(
          result->AddHist1D(branch+"_phiRes", "phi resolution", "#delta#phi", "Number of Tracks", 1000, phiMin, phiMax)
        );
    plots->phiRes_pt.push_back(
          result->AddHist2D(branch+"_phiRes_pt", "phi resolution", "#delta#phi", "Truth p_{T} [GeV]", 1000, phiMin, phiMax, 200, 0, 2000)
        );
    plots->phiRes_eta.push_back(
          result->AddHist2D(branch+"_phiRes_eta", "phi resolution", "#delta#phi", "Truth #eta", 1000, phiMin, phiMax, 200, -2, 2)
        );
    plots->phiRes_pt_eta.push_back(
          result->AddHist3D(branch+"_phiRes_pt_eta", "phi resolution", "#delta#phi", "Truth p_{T} [GeV]", "Truth #eta", 1000, phiMin, phiMax, 200, 0, 200, 100, -2, 2)
        );

    } // end of doResolutionPlots


  } // end loop over track branch names
  

  for(int i=0; i<jetBranchNames.size();++i){
    // one histogram for jet multiplicity
    std::string branch = jetBranchNames.at(i);
    std::vector<TH1*> tempVec;
    for(int j=0; j<N_JETS; ++j){
      std::string jetMulti = std::to_string(j+1);
      tempVec.push_back(
          result->AddHist1D( branch+"_jet"+jetMulti+"Pt", "", "", "", 1000, 0, 1000, 0, 0)
          );
    }
    plots->jetiPt.push_back(tempVec);
  } // end loop over jet branch names

  std::cout << "End BookHistograms" << std::endl;

}

//------------------------------------------------------------------------------

void AnalyseEvents(const int nEvents, bool hasPileup, ExRootTreeReader *treeReader, Plots *plots, std::vector<std::string> trackBranchNames, std::vector<std::string> jetBranchNames, bool doResolutionPlots)
{

  ///////////////////////////
  // Define branches
  ///////////////////////////

  if(m_debug) std::cout << "AnalyseEvents()" << std::endl;
  if(m_debug) std::cout << "Extracting branches" << std::endl;

  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");

  // Initialise track branches
  std::vector<TClonesArray*> trackBranches;
  for(const auto& bName : trackBranchNames){
    trackBranches.push_back( treeReader->UseBranch(bName) );
  }

  // Initialise jet branches
  std::vector<TClonesArray*> jetBranches;
  for(const auto& bName : jetBranchNames){
    jetBranches.push_back( treeReader->UseBranch(bName) );
  }



  ///////////////////////
  // Loop over all events
  ///////////////////////
  
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
    if( entry % 100==0 ) std::cout << "Event " << entry << " out of " << allEntries << std::endl; // print every 100 events


    ////////////////////////////////////////
    // Loop over the different track branches
    ////////////////////////////////////////
    
    if(m_debug) std::cout << "Being track loop" << std::endl;
    plots->numEntries->Fill(1);
    for(int iBranch=0; iBranch<trackBranchNames.size(); ++iBranch){
      TClonesArray* branch = trackBranches.at(iBranch);

      // Store the tracks, and sort into descending pT order
      int nTracks = branch->GetEntriesFast(); 
      std::vector<Track*> sortedTracks; 
      sortedTracks.reserve(nTracks);
      for(int iTrack=0; iTrack<nTracks; ++iTrack){
        Track* track = static_cast<Track*>( branch->At(iTrack) ); 
        sortedTracks.push_back(track);
      }
      std::sort( sortedTracks.begin(), sortedTracks.end(), reverseTrack);
      if(m_debug) std::cout << "\tExtracted tracks and sorted" << std::endl;


      // Loop over the sorted tracks 
      std::vector<Track*> leptons;
      std::vector<Track*> electrons;
      std::vector<Track*> muons;
      for(int iTrack=0; iTrack<nTracks; ++iTrack){
        Track* track = sortedTracks.at(iTrack); // ok to leave as Track*  
        if(iTrack == 0) plots->track1Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 1) plots->track2Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 2) plots->track3Pt.at(iBranch)->Fill(track->PT, eventWeight);
        if(iTrack == 3) plots->track4Pt.at(iBranch)->Fill(track->PT, eventWeight);

        plots->allTrackPt.at(iBranch)->Fill(track->PT, eventWeight);
        if(track->IsFake == 1){
          plots->fakeTrackPt.at(iBranch)->Fill(track->PT, eventWeight);
        }
        else{
          plots->trueTrackPt.at(iBranch)->Fill(track->PT, eventWeight);
        }

        // Calculate resolutions
        // Some of the might not make sense physically
        // Note: particle = static_cast<Candidate*>(candidate->GetCandidates()->At(0)); // this is "particle" 
        if(!hasPileup){
          GenParticle* particle = static_cast<GenParticle*>(track->Particle.GetObject());
          if(iTrack==0){
            plots->track1Eta.at(iBranch)->Fill(particle->Eta, eventWeight);
            plots->track1Pt_Eta.at(iBranch)->Fill(track->PT, particle->Eta, eventWeight);
          }
          plots->allTrackPt_Eta.at(iBranch)->Fill(track->PT, particle->Eta, eventWeight);
          plots->allTrackEta.at(iBranch)->Fill(particle->Eta, eventWeight);

          // Lepton plots
          int pid = particle->PID;
          if(abs(pid) == 11 || abs(pid) == 13){ // lepton
            leptons.push_back(track);
            //plots->lepton1Pt
            if( abs(pid) == 11 ){ // electron
              electrons.push_back(track);
            }
            if( abs(pid) == 13 ) { // muon
              muons.push_back(track);
            }
          }

          if(doResolutionPlots){
            double ptResolution = track->PT - particle->PT;  
            double z0Resolution = track->DZ - particle->DZ;
            double d0Resolution = track->D0 - particle->D0;
            double etaResolution = track->Eta - particle->Eta; 
            double cotThetaResolution = track->CtgTheta - particle->CtgTheta; 
            double phiResolution = track->Phi - particle->Phi;

            // Resolution plots
            plots->ptResRaw.at(iBranch)->Fill(ptResolution);
            plots->ptResRaw_pt.at(iBranch)->Fill(ptResolution, particle->PT);
            plots->ptResRaw_eta.at(iBranch)->Fill(ptResolution, particle->Eta);
            plots->ptResRaw_pt_eta.at(iBranch)->Fill(ptResolution, particle->PT, particle->Eta);

            plots->ptRes.at(iBranch)->Fill(ptResolution/particle->PT);
            plots->ptRes_pt.at(iBranch)->Fill(ptResolution/particle->PT, particle->PT);
            plots->ptRes_eta.at(iBranch)->Fill(ptResolution/particle->PT, particle->Eta);
            plots->ptRes_pt_eta.at(iBranch)->Fill(ptResolution/particle->PT, particle->PT, particle->Eta);

            plots->z0Res.at(iBranch)->Fill(z0Resolution);
            plots->z0Res_pt.at(iBranch)->Fill(z0Resolution, particle->PT);
            plots->z0Res_eta.at(iBranch)->Fill(z0Resolution, particle->Eta);
            plots->z0Res_pt_eta.at(iBranch)->Fill(z0Resolution, particle->PT, particle->Eta);

            plots->d0Res.at(iBranch)->Fill(d0Resolution);
            plots->d0Res_pt.at(iBranch)->Fill(d0Resolution, particle->PT);
            plots->d0Res_eta.at(iBranch)->Fill(d0Resolution, particle->Eta);
            plots->d0Res_pt_eta.at(iBranch)->Fill(d0Resolution, particle->PT, particle->Eta);

            plots->phiRes.at(iBranch)->Fill(phiResolution);
            plots->phiRes_pt.at(iBranch)->Fill(phiResolution, particle->PT);
            plots->phiRes_eta.at(iBranch)->Fill(phiResolution, particle->Eta);
            plots->phiRes_pt_eta.at(iBranch)->Fill(phiResolution, particle->PT, particle->Eta);
          }
          



        }
      } // end loop over tracks

      // fill plots for leptons 
      for(int iLep=0; iLep < leptons.size(); ++iLep){
        Track* lepton = leptons.at(iLep);
        if(iLep == 0) plots->lepton1Pt.at(iBranch)->Fill( lepton->PT, eventWeight );
        if(iLep == 1) plots->lepton2Pt.at(iBranch)->Fill( lepton->PT, eventWeight );
      }
      for(int iEle=0; iEle < electrons.size(); ++iEle){
        Track* electron = electrons.at(iEle);
        if(iEle == 0) plots->electron1Pt.at(iBranch)->Fill( electron->PT, eventWeight );
        if(iEle == 1) plots->electron2Pt.at(iBranch)->Fill( electron->PT, eventWeight );
      }
      for(int iMuon=0; iMuon < muons.size(); ++iMuon){
        Track* muon = muons.at(iMuon);
        if(iMuon == 0) plots->muon1Pt.at(iBranch)->Fill( muon->PT, eventWeight );
        if(iMuon == 1) plots->muon2Pt.at(iBranch)->Fill( muon->PT, eventWeight );
      }


    } // end loop over track branches
    

    ////////////////////////////////////////
    // Loop over the different jet branches
    ////////////////////////////////////////

    for(int iBranch=0; iBranch<jetBranches.size(); ++iBranch){
      TClonesArray* branch = jetBranches.at(iBranch);

      // Make sure jets are ordered in pT
      int nJets = branch->GetEntriesFast();
      std::vector<Jet*> sortedJets;
      sortedJets.reserve(nJets);
      for(int iJet = 0; iJet<nJets; ++iJet){
        Jet* jet = static_cast<Jet*>(branch->At(iJet));
        sortedJets.push_back(jet);
      }
      std::sort(sortedJets.begin(), sortedJets.end(), reverseJet);

      // Loop over sorted jets
      for(int iJet=0; iJet<nJets; ++iJet){
        Jet* jet = sortedJets.at(iJet);
        if(iJet < N_JETS){
          plots->jetiPt.at(iBranch).at(iJet)->Fill(jet->PT, eventWeight);
        }
      }

    } // end loop over jet branches


  } // end loop over entries

} // end AnalyseEvents



//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, Plots *plots)
{
  result->Print("pdf");
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{

  gROOT->SetBatch(1);
  //gROOT->ProcessLine("#include <vector>");

  std::string appName = "acceptances";
  std::cout << "Will execute " << appName << std::endl;

  std::string outputFile;
  std::string inputGlob; bool inputGlobSet(false);
  std::string inputFile; bool inputFileSet(false);

  // Get input agruments
  //std::string inputFile = argv[1]; 
  //std::string outputFile = argv[2];
  //int nEvents = atoi(argv[3]);
  //int hasPileup = atoi(argv[4]);
  int nEvents(-1);
  bool hasPileup(false);
  bool overwrite(false);
  bool doResolutionPlots(false);


  struct option longopts[] = {
    // These options set a flag
    {"overwrite",   no_argument,        0, 'w'},
    {"hasPileup",   no_argument,        0, 'p'},
    {"doResolution", no_argument,       0, 'r'},
    {"debug",   no_argument,            0, 'd'},
    // These options don't set a flag
    {"input",       required_argument,  0, 'i'},
    {"output",      required_argument,  0, 'o'},
    {"nEvents",     required_argument,  0, 'n'},
    {"inputGlob",   required_argument,  0, 'g'},
    {0,0,0,0}
  };

  int option_index = 0;
  int oc; 
  // if option letter is followed by a colon, then the option requires an argument
  while((oc = getopt_long(argc, argv, "wprdi:o:n:g:", longopts, &option_index)) != -1){
    switch(oc){
      case 'o':
        outputFile = optarg;  break;
      case 'w':
        overwrite = true;     break;
      case 'd':
        m_debug = true;     break;
      case 'i':
        inputFileSet = true;
        inputFile = optarg;   break;
      case 'n':
        nEvents = std::atoi(optarg);    break;
      case 'g':
        inputGlobSet = true;
        inputGlob = optarg;   break;
      case 'p':
        hasPileup = true;     break;
      case 'r':
        doResolutionPlots = true; break;
      case ':':
        std::cout << "WARNING: missing argument" << std::endl; 
        break;
      case '?':
        std::cout << "WARNING: unknown argument" << std::endl; 
        break;
    } // end of switch
  } // end of while(oc) 

  std::cout << "Input arguments:" << std::endl; 
  std::cout << "\tinputFile: " << inputFile << std::endl; 
  std::cout << "\toutputFile: " << outputFile << std::endl; 
  std::cout << "\tinputGlob: " << inputGlob << std::endl; 
  std::cout << "\tnEvents: " << nEvents << std::endl; 
  std::cout << "\toverwrite: " << overwrite << std::endl; 
  std::cout << "\thasPileup: " << hasPileup << std::endl; 

  if(!(inputGlobSet || inputFileSet)){
    std::cerr << "ERROR: need to choose either a single inputFile, or an input glob pattern" << std::endl; // not neither
    return 1;
  }
  if(inputGlobSet && inputFileSet){
    std::cerr << "ERROR: need to choose either a single inputFile, or an input glob pattern" << std::endl; // not both
    return 1;
  }


  // test if outputFile exists (exit if it does and overwrite not set)
  if(fileExists(outputFile) && overwrite){
    std::cout << "WARNING: " <<  outputFile << " file exists, will overwrite" << std::endl;
  }
  else if(fileExists(outputFile)){
    std::cerr << "ERROR: " <<  outputFile << " file exists, will NOT overwrite" << std::endl;
    return 0;
  }

  // info relating to input arguments
  if(nEvents < 0){
    std::cout << "INFO: Will run over all entries." << std::endl;
  }
  else{
    std::cout << "INFO: Will run over " << nEvents << " entries. " <<  std::endl;
  }


  // track branches to loop over
  std::vector<std::string> trackBranchNames = {
        // default delphes tracks
        "TruthTrack",
        "Track",
        "PBMatchedTracks",

        // tracks from hits 
        "TracksFromHit30", 
        "SmearedTracksFromHits",
        "PBMatchedHitTracks", 
  };

  // Jet branches to loop over
  std::vector<std::string> jetBranchNames = {
        // default delphes jets
        "TruthTrackJets",
        "SmearedTrackJets",
        "PBMatchedTrackJets",

        // Jets from tracks from hits 
        "TruthHitTrackJets",
        "SmearedHitTrackJets",
        "PBMatchedHitTrackJets",
  };

  // control analysis
  TChain *chain = new TChain("Delphes");
  if(inputFileSet){
    std::cout << "Adding " << inputFile << " to chain" << std::endl;
    chain->Add(inputFile.c_str());
  }
  else if(inputGlobSet){
    std::vector<std::string> fileList = glob(inputGlob);
    for(const auto& fName : fileList){
      std::cout << "Adding " << fName << " to chain" << std::endl;
      chain->Add(fName.c_str());
    }
  }


  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  Plots *plots = new Plots;
  BookHistograms(result, plots, trackBranchNames, jetBranchNames, doResolutionPlots);
  AnalyseEvents(nEvents, hasPileup, treeReader, plots, trackBranchNames, jetBranchNames, doResolutionPlots);

  // Turn off, never really want this ... 
  int doPrintHistograms = 0; 
  if(doPrintHistograms) PrintHistograms(result, plots);

  std::cout << "Writing to file: " << outputFile << std::endl;
  result->Write(outputFile.c_str(), overwrite);

  std::cout << "** Exiting..." << std::endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

  return 0;
}
