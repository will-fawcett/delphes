//#include "modules/TestClass.h"
//#include "modules/ClusteredJet.h"

// Delphes classes
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "modules/FastJetFinder.h"

// my analysis classes
#include "classes/TrackFitter.h"

// plotting stuff (shared)
#include "plotting.h"

// c++ libs
#include <iostream>
#include <sstream>
#include <exception>
#include <map>

// stuff for ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "TNamed.h"

// fastjet (presumably the version inside Delphes)
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

// attempt with eigen
//#include <Eigen/Core>
//#include "tricktrack/RiemannFit.h"
//

#include <ctime>


//------------------------------------------------------------------------------

void printTime(clock_t begin, clock_t end, std::string message){
  std::cout << message << ": time elapsed: " << double(end - begin) / CLOCKS_PER_SEC << std::endl;
}

//------------------------------------------------------------------------------

struct TrackStruct
{
  Float_t pT;
  Float_t d0;
  Float_t z0;
  Float_t eta;
  Float_t phi;
  Float_t kappa_013;
  Float_t kappa_123;
  Int_t isFake; 

  Float_t zresiduum;
  Float_t beamlineIntersect;


  // hit information
  Float_t hit3pT;
  Float_t hit3rho; 
  Float_t hit3eta;
  Float_t hit3phi; 

  Float_t hit2pT;
  Float_t hit2rho; 
  Float_t hit2eta;
  Float_t hit2phi; 

  Float_t hit1pT;
  Float_t hit1rho; 
  Float_t hit1eta;
  Float_t hit1phi; 

  Float_t z_phi12; 
  Float_t z_phi23; 
  Float_t z_phi13; 

  // don't forget to modify structDecoder string if this struct is edited

};



//------------------------------------------------------------------------------

struct TrackVectorStruct{
  std::vector<Float_t> pT;
  std::vector<Int_t> isFake;
  Int_t nTracks;
};

//------------------------------------------------------------------------------

TrackVectorStruct fillTrackVectorStruct(std::vector<myTrack>& theTracks){
  
  TrackVectorStruct trackHolder; 
  for(const myTrack& track : theTracks){

    trackHolder.pT.push_back( track.Pt() );
    trackHolder.isFake.push_back( track.isFake() );

  }
  trackHolder.nTracks = theTracks.size();
  return trackHolder;
}

//------------------------------------------------------------------------------

TrackStruct fillTrackStruct(const myTrack& track){ 
  // Function to fill a TrackStruct with the relevant info from the myTrack class 

  TrackStruct trackHolder;
  trackHolder.pT = track.Pt();
  trackHolder.d0 = track.D0();
  trackHolder.z0 = track.Z0();
  trackHolder.eta = track.Eta();
  trackHolder.phi = track.Phi();
  trackHolder.kappa_013 = track.kappa_bc();
  trackHolder.kappa_123 = track.kappa_nbc();
  trackHolder.isFake = track.isFake(); 

  trackHolder.zresiduum = track.GetZresiduum();
  trackHolder.beamlineIntersect = track.GetBeamlineIntersect();

  // hit information
  Hit * tempHit3 = track.GetHitAtLayer(2);
  Hit * tempHit2 = track.GetHitAtLayer(1);
  Hit * tempHit1 = track.GetHitAtLayer(0);

  trackHolder.hit3pT  = tempHit3->PT;
  trackHolder.hit3rho = tempHit3->HitRadius;
  trackHolder.hit3eta = tempHit3->Eta;
  trackHolder.hit3phi = tempHit3->Phi;

  trackHolder.hit2pT  = tempHit2->PT;
  trackHolder.hit2rho = tempHit2->HitRadius;
  trackHolder.hit2eta = tempHit2->Eta;
  trackHolder.hit2phi = tempHit2->Phi;

  trackHolder.hit1pT  = tempHit1->PT;
  trackHolder.hit1rho = tempHit1->HitRadius;
  trackHolder.hit1eta = tempHit1->Eta;
  trackHolder.hit1phi = tempHit1->Phi;

  // signed aangular differences in x--y plane 
  TVector3 p1(tempHit1->X, tempHit1->Y, 0); 
  TVector3 p2(tempHit2->X, tempHit2->Y, 0); 
  TVector3 p3(tempHit3->X, tempHit3->Y, 0); 
  TVector3 zhat(0, 0, 1); // z direction vector
  
  trackHolder.z_phi12 = p1.Cross(p2)*zhat / (p1.Mag() * p2.Mag());
  trackHolder.z_phi13 = p1.Cross(p3)*zhat / (p1.Mag() * p3.Mag());
  trackHolder.z_phi23 = p2.Cross(p3)*zhat / (p2.Mag() * p3.Mag());

  return trackHolder; 
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

// using hitContainer = std::map<int, std::vector<Hit*> >; 
hitContainer fillHitContainer(TClonesArray* hitBranch){
  // Fill a  hitContainer object with the hits form a particular branch

    hitContainer hc; 
    for(auto itHit = hitBranch->begin(); itHit != hitBranch->end(); ++itHit){
      Hit * hit = dynamic_cast<Hit*>(*itHit);
      int SurfaceID = hit->SurfaceID; 
      hc[SurfaceID].push_back(hit); 
    }
    return hc;
}

//------------------------------------------------------------------------------

int countFakes(std::vector<myTrack>& theTracks){

    // count number of fake tracks 
    int nFakes(0);
    for(const auto& track : theTracks){
      if(track.isFake()) nFakes++;
    }
    return nFakes;
}

//------------------------------------------------------------------------------

void AnalyseEvents(const int nEvents, ExRootTreeReader *treeReader, Plots *plots, float vertexDistSigma, int nVertexSigma, std::string oFileName, float zresiduumTolerance)
{

  // Define branches
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  //TClonesArray *branchTruthTrack = treeReader->UseBranch("TruthTrack");
  TClonesArray *branchTrack      = treeReader->UseBranch("Track");

  // Different hit branches, corresponding to the different geometry 
  TClonesArray *branchHit50        = treeReader->UseBranch("Hits50");
  TClonesArray *branchHit40        = treeReader->UseBranch("Hits40");
  TClonesArray *branchHit30        = treeReader->UseBranch("Hits30");
  TClonesArray *branchHit20        = treeReader->UseBranch("Hits20");
  TClonesArray *branchHit10        = treeReader->UseBranch("Hits10");

  // Declare track structs 
  //TrackVectorStruct tracks50;
  //TrackVectorStruct tracks40;
  //TrackVectorStruct tracks30;
  //TrackVectorStruct tracks20;
  //TrackVectorStruct tracks10;
  
  TrackStruct tracks50;
  TrackStruct tracks40;
  TrackStruct tracks30;
  TrackStruct tracks20;
  TrackStruct tracks10;

  // Try to add a new ROOT file 
  std::string newOFileName = oFileName.erase( oFileName.size()-5, 5); // remove last 5 charachters ".root"
  newOFileName += "_tracks.root";
  std::cout << "New string name: " << newOFileName << std::endl;

  TFile* oFile = new TFile(newOFileName.c_str(), "RECREATE");
  TTree* tree10 = new TTree("Tracks10", "A tree with tracks 10");
  TTree* tree20 = new TTree("Tracks20", "A tree with tracks 20");
  TTree* tree30 = new TTree("Tracks30", "A tree with tracks 30");
  TTree* tree40 = new TTree("Tracks40", "A tree with tracks 40");
  TTree* tree50 = new TTree("Tracks50", "A tree with tracks 50");

  TString structDecoder = "pT/F:d0/F:z0/F:eta/F:phi/F:kappa_013/F:kappa_123/F:isFake/I:zresiduum/F:beamlineIntersect/F";
  structDecoder += ":hit3pT/F:hit3rho/F:hit3eta/F:hit3phi/F";
  structDecoder += ":hit2pT/F:hit2rho/F:hit2eta/F:hit2phi/F";
  structDecoder += ":hit1pT/F:hit1rho/F:hit1eta/F:hit1phi/F";
  structDecoder += ":z_phi12/F:z_phi23/F:z_phi13/F";

  TBranch* branchTrack10 = tree10->Branch("tracks10", &tracks10, structDecoder);
  TBranch* branchTrack20 = tree20->Branch("tracks20", &tracks20, structDecoder);
  TBranch* branchTrack30 = tree30->Branch("tracks30", &tracks30, structDecoder);
  TBranch* branchTrack40 = tree40->Branch("tracks40", &tracks40, structDecoder);
  TBranch* branchTrack50 = tree50->Branch("tracks50", &tracks50, structDecoder);

  //TTree testTree = new Tree("something","TestTree");
  //Holder holder;
  //testTree->Branch("event.",&holder);
  //testTree->Fill();
  //testTree->ResetBranchAddresses();

  // testing
  /**********
  TrackHolder * th = new TrackHolder();
  TrackHolder holder;
  //TBranch* testBranch = tree->Branch("testBranch", "TrackHolder", &holder, 64000, 0);
  TBranch* testBranch = tree->Branch("event.", &holder);

  std::vector<Float_t> testpT = {0.1, 0.2, 0.3};
  std::vector<Int_t> testIsFake = {0, 1, 1};
  Int_t test_nTracks = 3;
  holder.pT = testpT;
  holder.isFake = testIsFake;
  holder.nTracks = test_nTracks; 

  tree->Fill();
  **********/




  // based on the spread of vertices, calculate the maximum tolerance along the beam line to which a track can point
  float maxZ = nVertexSigma*vertexDistSigma;
  float minZ = -1*maxZ; 
  std::vector<float> parameters;
  parameters.push_back(minZ);
  parameters.push_back(maxZ);

  // Loop over all events
  Long64_t allEntries = treeReader->GetEntries();
  int nEventsCorrectlyIdentifiedVertex(0);
  std::cout << "** Chain contains " << allEntries << " events" << std::endl;
  // an event weight for normalising histograms
  std::vector<int> layerIDs;
  float eventWeight = 1.0/allEntries; 
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



    //////////////////////////////////////////
    // Count the number of tracks from Delphes
    //////////////////////////////////////////
    
    int nTracks(0);
    int nTracks1GeV(0);
    int nTracks10GeV(0);
    for(auto itTrack=branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      Track* track = dynamic_cast<Track*>(*itTrack);
      nTracks++;
      if(track->PT < 1.0) continue;
      nTracks1GeV++;
      if(track->PT < 10) continue;
      nTracks10GeV++;
    }
    plots->nDelphesTracks->Fill(nTracks, eventWeight);
    plots->nDelphesTracks1GeV->Fill(nTracks1GeV, eventWeight);
    plots->nDelphesTracks10GeV->Fill(nTracks10GeV, eventWeight);

    
    // Fill the hitContainer with all the hits in the event (geometry defined by the surfaces) 
    //clock_t begin_fillHitContainer = clock();
    hitContainer hc50 = fillHitContainer(branchHit50); 
    //clock_t end_fillHitContainer = clock();
    //printTime(begin_fillHitContainer, end_fillHitContainer, "Fill hit container");
    hitContainer hc40 = fillHitContainer(branchHit40); 
    hitContainer hc30 = fillHitContainer(branchHit30); 
    hitContainer hc20 = fillHitContainer(branchHit20); 
    hitContainer hc10 = fillHitContainer(branchHit10); 


    // get list of unique layer IDs (only need filling once)
    if(entry == 0){
      for(auto const& key : hc50){
        layerIDs.push_back(key.first);
      }
      std::sort(layerIDs.begin(), layerIDs.end()); // layers should be ascengin 
      std::cout << "Sorted Layer IDs: " << std::endl;
      for(auto l : layerIDs) std::cout << "\t" << l << std::endl;
    }


    // contains the rules to associate hits together, and then the collection of hits into tracks 
    //tf.debug();
    
    std::vector<hitContainer> theHitContainers = {hc10, hc20, hc30, hc40, hc50};
    int counter(0);
    for(auto& hc : theHitContainers){

      // fill histogram with the pT of hits in the outermost layer
      int nHitsPt2(0);
      for(auto& hit : hc[ layerIDs.back() ]){
        plots->nHitsPt.at(counter)->Fill(hit->PT);
        plots->nHitsEta.at(counter)->Fill(hit->Eta);
        if(hit->PT > 2.0){
          nHitsPt2++;
          plots->nHitsEtaPt2.at(counter)->Fill(hit->Eta);
        }
      }

      // number of hits with pT > 2 GeV
      plots->nDelphesHitsPt2.at(counter)->Fill(nHitsPt2); 

      //float tolerance = 0.4; // mm  2*sqrt(0.04) = 0.4 , 0.04 mm is size of pixel 
      //float tolerance = 0.5; // [mm]  2*sqrt(0.04) = 0.4 , 0.04 mm is size of pixel, got to 0.5 to be conservative (this is a preselection)
      TrackFitter tf(fitTypes::simpleLinear, parameters, layerIDs, zresiduumTolerance); 
      bool associated = tf.AssociateHits(hc);

      if(associated){
        std::vector< myTrack > baseTracks = tf.GetTracks();
        /*********
        std::cout << "Number of base tracks: " << baseTracks.size() << std::endl;
        int nBaseFakes = countFakes(baseTracks);
        std::cout << "\tOf which " << nBaseFakes << " are fake." << std::endl;
        *******/
        
        for(const myTrack& track : baseTracks){
          
          //if(track.isFake()) continue; // only true tracks

          if(counter==0){
            tracks10 = fillTrackStruct(track); 
            tree10->Fill();
          }
          else if(counter==1){
            tracks20 = fillTrackStruct(track);
            tree20->Fill();
          }
          else if(counter==2){
            tracks30 = fillTrackStruct(track);
            tree30->Fill();
          }
          else if(counter==3){
            tracks40 = fillTrackStruct(track);
            tree40->Fill();
          }
          else if(counter==4){
            tracks50 = fillTrackStruct(track);
            tree50->Fill();
          }
        }


        //tf.ApplyCurvatureCut( 0.01 ); // reject tracks with curvature difference greater than this number [mm^-1]
        //tf.ApplyCurvatureCut( 0.005 ); // reject tracks with curvature difference greater than this number [mm^-1]

        // first argument: vector of pT thresholds. Second argument: deltaKappa threshoolds
        tf.ApplyPtDependantCurvatureCut( {0, 2, 10, 50}, {0.005, 0.004, 0.0025, 0.001} );

        std::vector< myTrack > theTracks = tf.GetTracks(); 
        int nFakes = countFakes(theTracks);

        /***************
        std::cout << "Number of tacks after cut: " << theTracks.size() << std::endl;
        std::cout << "\tOf which " << nFakes << " are fake." << std::endl;
        std::cout << "Removed " << baseTracks.size()-theTracks.size() << " tracks" << std::endl;
        std::cout << "Removed " << nBaseFakes - nFakes << " fakes" << std::endl;
        std::cout << "" << std::endl;
        ********************/

        std::vector<Hit*> outerHits = hc[2];
        plots->nDelphesHits.at(counter)->Fill( outerHits.size(), eventWeight); // number of hits in outer layer
        plots->nRecoTracks.at(counter)->Fill(theTracks.size(), eventWeight);
        plots->nRecoTracksMatched.at(counter)->Fill(theTracks.size() - nFakes, eventWeight);

        int nTracksPt2(0);
        int nMatchedTracksPt2(0);

        // loop over tracks to fill fake rate plots
        for(const auto& track : theTracks){

          //track.printHitInfo();

          plots->recoTrackPt.at(counter)->Fill(track.Pt());
          plots->recoTrackEta.at(counter)->Fill(track.Eta()); // might want to check theta calculation
          //if(track.Pt()>2) nTracksPt2++;
          float outerHitPt = track.GetHitPtAtLayer(2);
          float outerHitEta = track.GetHitAtLayer(2)->Eta; 
          if(outerHitPt > 2) nTracksPt2++;

          float curvatureDifference =  track.kappa_bc() - track.kappa_nbc();
          if(track.isFake()){
            // FAKE tracks 
            plots->recoTrackPt_fake.at(counter)->Fill(track.Pt());
            plots->recoTrackEta_fake.at(counter)->Fill(track.Eta()); 
            if(counter == 1){
              plots->curvatureDifference_fake->Fill(curvatureDifference); 
              static int fakePointCounter = 0;
              plots->curvatureDifference_curvature_fake->SetPoint(fakePointCounter,curvatureDifference , track.kappa_bc());
              //std::cout << "fakePointCounter: " << fakePointCounter << std::endl;
              fakePointCounter++;
            }
          }
          else{
            // TRUE tracks 
            plots->recoTrackPtResolution.at(counter)->Fill(track.Pt() - outerHitPt); 
            plots->recoTrackPt_true.at(counter)->Fill(track.Pt());
            plots->recoTrackHitPt_true.at(counter)->Fill(outerHitPt);
            plots->recoTrackHitEta_true.at(counter)->Fill(outerHitEta);
            plots->recoTrackEta_true.at(counter)->Fill(track.Eta()); 
            //if(track.Pt()>2) nMatchedTracksPt2++;
            if(track.GetHitPtAtLayer(2) > 2){
              nMatchedTracksPt2++; 
              plots->recoTrackHitEtaPt2_true.at(counter)->Fill(outerHitEta);
            }

            if(counter == 1){
              plots->curvatureDifference_true->Fill(curvatureDifference); 
              static int truePointCounter = 0;
              plots->curvatureDifference_curvature_true->SetPoint(truePointCounter, curvatureDifference, track.kappa_bc());
              truePointCounter++;
            }
          }
        } // end of loop over tracks 
        plots->nRecoTracksPt2.at(counter)->Fill(nTracksPt2);
        plots->nRecoTracksMatchedPt2.at(counter)->Fill(nMatchedTracksPt2);

      } // end of if(associated)
      //std::cout << "i: " << counter << std::endl;
      ++counter;
    } // end loop over hitContainers 

  } // end loop over entries
  tree10->Write();
  tree20->Write();
  tree30->Write();
  tree40->Write();
  tree50->Write();
} // end AnalyseEvents



//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, Plots *plots)
{
  result->Print("pdf");
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  //Matrix3xNd hits = Matrix3xNd::Random(3,3);
  //Matrix3Nd hits_cov = Matrix3Nd::Random(9,9);

  gROOT->SetBatch(1);
  gROOT->ProcessLine("#include <vector>");

  // configurable(>) parameters
  float vertexDistSigma = 53.0; // standard deviation (1 sigma) of the distribution of vertices along z
  int nVertexSigma = 4;

  std::string appName = "basicTrackReconstructon";
  std::cout << "Will execute " << appName << std::endl;
  std::string inputFile = argv[1]; // doesn't complain about cast? Maybe compiler can deal with it :p
  std::string outputFile = argv[2];
  //int doPrintHistograms = atoi(argv[3]);
  int doPrintHistograms = 0; 
  float zresiduumTolerance = atof(argv[3]);
  nVertexSigma = atoi(argv[4]);
  int nEvents(-1);

  if(argc > 5){
    nEvents = atoi(argv[5]);
    std::cout << "INFO: Will run over " << nEvents << std::endl;
  }
  // 
  if(argc < 5)
  {
    std::cout << " Usage: " << appName << " input_file output_file [optional: nEvents]" << std::endl;
    std::cout << " 1) input_file: ROOT file containing delphes output,"       << std::endl;
    std::cout << " 2) output_file: ROOT file name that will contain output histograms" << std::endl;
    //std::cout << " 3) bool: write files to .pdf" << std::endl;
    std::cout << "3) float: tolerance for the beamline constraint (recommended 0.5 [mm])" << std::endl; 
    std::cout << " 4) number of sigma of for the beamline constraint" << std::endl;
    std::cout << " 5) nEvents (optional): number of events to be processed (for testing, default is all events)." << std::endl;
    return 1;
  }

  // control analysis
  TChain *chain = new TChain("Delphes");
  std::cout << "Adding " << inputFile << " to chain" << std::endl;
  chain->Add(inputFile.c_str());



  /*************
  // test with a struct
  TestStruct track;
  track.pT = 5.0;
  track.z0 = 10;
  track.kappa = 0.001; 

  tree->Branch("track", &track, "pT/F:z0/F:kappa/F");
  tree->Fill();
  track.pT = 10.0;
  track.z0 = 11.0;
  track.kappa = 0.02;
  tree->Fill(); 
  tree->Write();
  oFile->Close();
  delete oFile;
  *************/


  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  Plots *plots = new Plots;
  BookHistograms(result, plots);
  AnalyseEvents(nEvents, treeReader, plots, vertexDistSigma, nVertexSigma, outputFile, zresiduumTolerance);

  if(doPrintHistograms) PrintHistograms(result, plots);

  std::cout << "Writing to file: " << outputFile << std::endl;
  result->Write(outputFile.c_str());

  std::cout << "** Exiting..." << std::endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

  return 0;
}
