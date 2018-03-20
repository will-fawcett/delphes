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

// fastjet (presumably the version inside Delphes)
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;

//------------------------------------------------------------------------------

std::string PrintTLorentz(TLorentzVector &v){
  std::ostringstream s;
  s <<  " (" << v.Pt() << ", " << v.Eta() << ", " << v.Phi() << ", " << v.M() << ")";
  return s.str();
}

std::string PrintPseudoJet(PseudoJet &j){
  std::ostringstream s;
  s <<  " (" << j.pt() << ", " << j.eta() << ", " << j.phi_std() <<  ")";
  return s.str();
}

//------------------------------------------------------------------------------
//


using MyClusteredJet = std::pair<PseudoJet, std::vector<float> >;
using JetWithTracks = std::pair<PseudoJet, std::vector<Track*> >;

// Rule for sorting vector of tracks 
bool reverse(const Track* i, const Track* j){
  return i->PT > j->PT; 
}

// Add tracks to MyClusteredJet
JetWithTracks addTracksToJet(const MyClusteredJet& clusteredJet, TClonesArray* branchTrack){

  // use unique ID to match the track
  std::vector<float> trackIDs = clusteredJet.second;
  std::vector<Track*> matchedTracks;


  //std::cout << "addTracksToJet: number of input tracks: " << trackIDs.size() << std::endl;
  for(auto trackID : trackIDs){
    for(auto itTrack=branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      Track* track = dynamic_cast<Track*>(*itTrack);
      if(track->GetUniqueID() == trackID){
        matchedTracks.push_back(track);
        break; // once the track is matched, break
        // TODO: There is more than 1 track with the SAME uniqueID !!!!!!!!!!!!!!!!!!!!!!!!!
        // Might not be assigning the "right" track!!!!!!!!!!!!!!!
      }
    }
  }
  //std::cout << "addTracksToJet: number of matched tracks: " << matchedTracks.size() << std::endl;

  // Sort matched tracks
  std::sort( matchedTracks.begin(), matchedTracks.end(), reverse ); // greater ensures sorting is descending, for pT

  JetWithTracks output;
  output.first = clusteredJet.first;
  output.second = matchedTracks;
  return output;
}


std::vector<JetWithTracks> addTracksToJet(const std::vector<MyClusteredJet>& clusteredJets, TClonesArray* branchTrack){
    
  std::vector<JetWithTracks> outputJetsWithTracks;
  for(auto aClusteredJet : clusteredJets){
    JetWithTracks temp = addTracksToJet(aClusteredJet, branchTrack);
    outputJetsWithTracks.push_back(temp);
  }
  return outputJetsWithTracks;
}

// Use FastJetFinder inside delphes to cluster jets
std::vector<MyClusteredJet> jetCluster(std::vector<Track*> inputTracks, JetDefinition* definition, const float minJetPt, const float minTrackPt, bool debug){

    // Convert tracks to PseudoJets, dress with unique ID
    std::vector<PseudoJet> inputList;
    for(auto track : inputTracks){
      TLorentzVector momentum = track->P4();
      if(momentum.Pt() < minTrackPt) continue; // only add tracks above a certain jet pT thresold
      PseudoJet jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
      jet.set_user_index(track->GetUniqueID()); // add info to the jet
      inputList.push_back(jet);
    }

    if(debug) std::cout << "Number of input tracks to clustering (after min track pT selection): " << inputList.size() << std::endl;

    // run clustering
    ClusterSequence sequence(inputList, *definition);

    // sort pT
    std::vector<PseudoJet> outputList = sorted_by_pt(sequence.inclusive_jets(minJetPt));
    if(debug) std::cout << "There are " << outputList.size() << " output jets" << std::endl;


    //////////////////////
    // Extract tracks uIDs that were associated to a jet
    //////////////////////
    std::vector<MyClusteredJet> outputJets;
    for(auto fastJet : outputList){
      if(debug) std::cout << PrintPseudoJet(fastJet) << std::endl;

      // extract tracks that were actually associated to jets
      std::vector<PseudoJet> clusteredTracks = sequence.constituents(fastJet);
      if(debug) std::cout << "\tThis jet had: " << clusteredTracks.size() << " tracks." << std::endl;
      std::vector<float> uIDs;
      for(auto track : clusteredTracks){
        //std::cout << "ID: " << track.user_index() << std::endl;
        uIDs.push_back(track.user_index());
      }

      //MyClusteredJet outputJet(fastJet, uIDs);
      MyClusteredJet outputJet;
      outputJet.first  = fastJet;
      outputJet.second = uIDs;
      outputJets.push_back(outputJet);
    }
    
    return outputJets;
}

std::vector<MyClusteredJet> jetCluster(TClonesArray* branchTrack, JetDefinition* definition, const float minJetPt, const float minTrackPt, bool debug){
  std::vector<Track*> outputTracks;
  int nIteratedTracks(0);
  for(auto itTrack = branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
    outputTracks.push_back( dynamic_cast<Track*>(*itTrack) );
    nIteratedTracks++;
  }
  //std::cout << "nIteratedTracks: " << nIteratedTracks << std::endl;
  return jetCluster(outputTracks, definition, minJetPt, minTrackPt, debug);
}


//------------------------------------------------------------------------------

void PrintTrack(Track *track)
{
  std::cout << "pT "       << track->PT       << std::endl;
  std::cout << "eta "      << track->Eta      << std::endl;
  std::cout << "phi "      << track->Phi      << std::endl;
  std::cout << "CotTheta " << track->CtgTheta << std::endl;
  std::cout << "Charge "   << track->Charge   << std::endl;
  std::cout << "DZ "       << track->DZ       << std::endl;
  std::cout << "D0 "       << track->D0       << std::endl;
}

inline float calculatDeltaPhi(float phi1, float phi2){
  return acos( cos( phi1 - phi2 ));
}

//------------------------------------------------------------------------------

std::map<std::string, float> findPrimaryBin(TClonesArray* branchTrack, float binWidth, float slideStep, int beamMinZ, int beamMaxZ){

    /*************************
     *
     * Simple fixed bin algorithm. May split PV
     *
    // let bin width be 1 mm ?
    TH1F * eventBinnedZpT = new TH1F("eventBinnedZpT", "", 600, -300, 300); // for 5mm use 120 bins
    for(auto itTrack = branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      float trackPt = static_cast<Track*>(*itTrack)->PT;
      float zPosition = static_cast<Track*>(*itTrack)->DZ;
      plots->binnedZpT->Fill(zPosition, trackPt);
      eventBinnedZpT->Fill(zPosition, trackPt);
    } // end iteration over tracks

    // get bin with largest sum(pT)
    Int_t binWithPtMax = eventBinnedZpT->GetMaximumBin();

    // get z range of bin with largest sum(pT)
    TAxis * xaxis = static_cast<TAxis*>(eventBinnedZpT->GetXaxis());
    float zMin = xaxis->GetBinLowEdge(binWithPtMax);
    float zMax = xaxis->GetBinUpEdge(binWithPtMax);
    float zWidth = xaxis->GetBinWidth(binWithPtMax);
    delete eventBinnedZpT;

    *
    *************************/

    int nBins = (beamMaxZ-beamMinZ)/binWidth;

    // define histograms for sliding window algorithm
    //std::vector<TH1F*> windowHists;
    //std::vector<TH1F*> nTrackHists;
    std::vector< std::pair< TH1F*, TH1F*> > hists;
    for(int i=0; i< binWidth/slideStep; ++i){

      //windowHists.push_back( new TH1F( (std::to_string(i)+"window").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i) );
      //nTrackHists.push_back( new TH1F( (std::to_string(i)+"nTrack").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i) );
      TH1F * windowHist = new TH1F( (std::to_string(i)+"window").c_str(),  "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i );
      TH1F * nTrackHist = new TH1F( (std::to_string(i)+"nTrack").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i );
      hists.push_back( std::make_pair(windowHist, nTrackHist) );
    }

    //std::cout << "nHistos: " << binWidth/slideStep << " nbins: " << nBins << " maxZ: " << beamMaxZ << " minZ: " << beamMinZ << std::endl;

    // loop over all tracks, fill histograms
    for(auto itTrack = branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      float trackPt   = dynamic_cast<Track*>(*itTrack)->PT;
      float zPosition = dynamic_cast<Track*>(*itTrack)->DZ;
      //for(auto hist : windowHists)  hist->Fill(zPosition, trackPt);
      for(auto histPair : hists){
        histPair.first->Fill(zPosition, trackPt);
        histPair.second->Fill(zPosition, 1); // number of tracks, maybe could have used a TH2D?
      }
    }

    // Extract the largest sum(pT)
    float previousMaxPt(0), zBinLow(0), zBinHigh(0), zBinWidth(0), zCentroid(0);
    int nTracks(0);
    //for(auto hist : windowHists){
    for(auto histPair : hists){
      Int_t binWithPtMax = histPair.first->GetMaximumBin();
      TAxis * xaxis = static_cast<TAxis*>(histPair.first->GetXaxis());
      float maxPt = histPair.first->GetBinContent(binWithPtMax);
      if(maxPt > previousMaxPt){
        previousMaxPt = maxPt;
        zBinLow   = xaxis->GetBinLowEdge(binWithPtMax);
        zBinHigh  = xaxis->GetBinUpEdge(binWithPtMax);
        zBinWidth = xaxis->GetBinWidth(binWithPtMax);
        zCentroid = xaxis->GetBinCenter(binWithPtMax);
        //std::cout << "low " << zBinLow << " \thigh " << zBinHigh << "\t center " << zCentroid << std::endl;
        nTracks = histPair.second->GetBinContent(binWithPtMax);
      }
    }

    // store results
    std::map<std::string, float> results;
    results["zMin"]   = zBinLow;
    results["zMax"]   = zBinHigh;
    results["zWidth"] = zBinWidth;
    results["zPositon"] = zCentroid;
    results["nTracks"] = nTracks;

    // delete histograms
    /********
    for(auto hist : windowHists){
      delete hist;
    }
    windowHists.clear();
    ********/
    for(auto histPair : hists){
      delete histPair.first;
      delete histPair.second;
    }
    hists.clear();

    return results;

} // end of find primary bin


//------------------------------------------------------------------------------

struct TestPlots
{

  TH1 *z0Res;
  TH2 *z0Res_pt;
  TH2 *z0Res_eta;
  TH3 *z0Res_pt_eta;

  TH1 *d0Res;
  TH2 *d0Res_pt;
  TH2 *d0Res_eta;
  TH3 *d0Res_pt_eta;

  TH1 *CtgThetaRes;
  TH2 *CtgThetaRes_pt;
  TH2 *CtgThetaRes_eta;
  TH3 *CtgThetaRes_pt_eta;

  TH1 *phiRes;
  TH2 *phiRes_pt;
  TH2 *phiRes_eta;
  TH3 *phiRes_pt_eta;

  TH1 *pt;
  TH1 *logpt;
  TH1 *truthPt;
  TH1 *ptRes;
  TH2 *ptRes_pt;
  TH2 *ptRes_eta;
  TH3 *ptRes_pt_eta;

  TH1 *ptResRaw;
  TH2 *ptResRaw_pt;
  TH2 *ptResRaw_eta;
  TH3 *ptResRaw_pt_eta;

  TH1 *eta;
  TH1 *etaRes;
  TH2 *etaRes_eta;
  TH3 *etaRes_pt_eta;

  TH2 *track_eta_phi_pt;
  TH2 *jet_eta_phi_pt;

  // For triggers? 
  std::vector<TH1*> truthTrackNpT;
  TH1 *truthTrackPt100;

  // Track jets from Delphes
  TH1 *allJetPt;
  TH1 *nJets;
  std::vector<TH1*> jetNPt;
  std::vector<TH1*> jetNEta;
  std::vector<TH1*> jetNPhi;
  std::vector<TH1*> jetNTrackPt;
  std::vector<TH1*> jetNTrackMulti;

  // Jets clustered from all Delphes tracks
  TH1 *nominalJets;
  std::vector<TH1*> nominalJetNPt;
  std::vector<TH1*> nominalJetNEta;
  std::vector<TH1*> nominalJetNPhi;
  std::vector<TH1*> nominalJetNTrackPt;
  std::vector<TH1*> nominalJetNTrackMulti;

  // Jets clustered from tracks only associated to the PB
  TH1 *nAssociatedJets;
  std::vector<TH1*> associatedJetNPt;
  std::vector<TH1*> associatedJetNEta;
  std::vector<TH1*> associatedJetNPhi;
  std::vector<TH1*> associatedJetNTrackPt;
  std::vector<TH1*> associatedJetNTrackMulti;

  // Primary vertices
  TH1 *binnedZpT;
  TH1 *binnedZnVertices;

  TH1 *misidPV;
  TH1 *misidPVLogx;

  //TH1 *trackMultiplicityInPB;
  //TH1 *vertexMultiplicityInPB;
  //TH2 *nVertexVsZVertexPosition;
  //TH2* nVertexVsZPBPosition;
  //TH2* PBvZVertexPosition;
  
  // Occupancy plots
  TH2 *trackOccupancy;
  TProfile *trackOccupancyProf;

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots, bool calculateTrackParameters)
{

  // Track occupancy
  plots->trackOccupancy = result->AddHist2D("trackOccupancy", "", "Track multi per tower", "Track pT [GeV]", 100, 0, 100, 100, 0, 300);
  plots->trackOccupancyProf = result->AddProfile("trackOccupancyProf", "", "Track multi per tower", "Number of towers", 100, 0, 100);



  // trigger
  for(int i=0; i<19; ++i){
    int pTcut = 5*(i+1);
    plots->truthTrackNpT.push_back(
        result->AddHist1D("truthTrackPt"+std::to_string(pTcut), "", "Truth Track pT after offline cut of "+std::to_string(pTcut), "", 100, 0, 100)
        );
  }
  plots->truthTrackPt100 =result->AddHist1D("truthTrack100", "", "", "",  100, 0, 100, 0, 0);

  plots->misidPVLogx = result->AddHist1D("misidPVLogx", "Misidentified primary vertices", "Distance between PB and true PV [mm]", "Number of events", 200, 0, 400, 1, 0); // probably will be some overflow
  plots->misidPV = result->AddHist1D("misidPV", "Misidentified primary vertices", "Distance between PB and true PV [mm]", "Number of events", 200, 0, 200, 0, 0); // probably will be some overflow

  // primary bin 
  //plots->trackMultiplicityInPB    = result->AddHist1D("trackMultiplicityInPB", "Track multiplicity in PB", "", "", 100, 0, 100);
  //plots->vertexMultiplicityInPB   = result->AddHist1D("vertexMultiplicityInPB", "Track multiplicity in PB", "", "", 100, 0, 100);
  //plots->nVertexVsZVertexPosition = result->AddHist2D("nVertexVsZVertexPosition", "Number of vertices in PB", "z position of true PV", 100, 0, 100, 100, -300, 300);
  //plots->nVertexVsZPBPosition     = result->AddHist2D("nVertexVsZPBPosition", "Number of vertices in PB", "z position of PB", 100, 0, 100, 100, -300, 300);
  //plots->PBvZVertexPosition       = result->AddHist2D("PBvZVertexPosition", "z position of PB", "z position of PV", 100, -300, 300, 100, -300, 300);


  // jet histograms (from Delphes track jets)
  for(int i=0; i<7; ++i){
    std::string i_str = std::to_string(i+1);

    plots->jetNPt.push_back(
        result->AddHist1D( "jet"+i_str+"Pt", "", i_str+" Jet p_{T} [GeV]", "",  1000, 0, 1000, 0, 0)
        );
    plots->jetNEta.push_back(
        result->AddHist1D( "jet"+i_str+"Eta", "", i_str+" Jet #eta", "",  100, -2.5, 2.5)
        );
    plots->jetNPhi.push_back(
        result->AddHist1D( "jet"+i_str+"Phi", "", i_str+" Jet #phi", "",  100, -M_PI, M_PI)
        );
    // track-inside-jet properties
    plots->jetNTrackPt.push_back(
      result->AddHist1D( "jet"+i_str+"TrackPt", "", "p_{T} of tracks inside jet "+i_str+" [GeV]", "", 100, 0, 20)
        );
    plots->jetNTrackMulti.push_back(
      result->AddHist1D( "jet"+i_str+"TrackMulti", "", "Multiplicity of tracks inside jet "+i_str, "", 100, 0, 100)
        );

    // Jets clustered from all tracks (outside of Delphes)
    plots->nominalJetNPt.push_back(
        result->AddHist1D( "nominalJet"+i_str+"Pt", "", i_str+" Jet p_{T} [GeV]", "",  1000, 0, 1000, 0, 0)
        );
    plots->nominalJetNEta.push_back(
        result->AddHist1D( "nominalJet"+i_str+"Eta", "", i_str+" Jet #eta", "",  100, -2.5, 2.5)
        );
    plots->nominalJetNPhi.push_back(
        result->AddHist1D( "nominalJet"+i_str+"Phi", "", i_str+" Jet #phi", "",  100, -M_PI, M_PI)
        );
    // track-inside-jet properties for nominal jets
    plots->nominalJetNTrackPt.push_back(
      result->AddHist1D( "nominalJet"+i_str+"TrackPt", "", "p_{T} of tracks inside jet "+i_str+" [GeV]", "", 100, 0, 20)
        );
    plots->nominalJetNTrackMulti.push_back(
      result->AddHist1D( "nominalJet"+i_str+"TrackMulti", "", "Multiplicity of tracks inside jet "+i_str, "", 100, 0, 100)
        );

    // Jets associated to PV
    plots->associatedJetNPt.push_back(
        result->AddHist1D( "associatedJet"+i_str+"Pt", "", i_str+" Jet p_{T} [GeV]", "",  1000, 0, 1000, 0, 0)
        );
    plots->associatedJetNEta.push_back(
        result->AddHist1D( "associatedJet"+i_str+"Eta", "", i_str+" Jet #eta", "",  100, -2.5, 2.5)
        );
    plots->associatedJetNPhi.push_back(
        result->AddHist1D( "associatedJet"+i_str+"Phi", "", i_str+" Jet #phi", "",  100, -M_PI, M_PI)
        );
    // track-inside-jet properties for associated jets
    plots->associatedJetNTrackPt.push_back(
      result->AddHist1D( "associatedJet"+i_str+"TrackPt", "", "p_{T} of tracks inside jet "+i_str+" [GeV]", "", 100, 0, 20)
        );
    plots->associatedJetNTrackMulti.push_back(
      result->AddHist1D( "associatedJet"+i_str+"TrackMulti", "", "Multiplicity of tracks inside jet "+i_str, "", 100, 0, 100)
        );

  }

  plots->nJets           = result->AddHist1D( "nJets", "nJets", "Number of Jets", "", 100, 0, 100, 0, 0 );
  plots->nAssociatedJets = result->AddHist1D( "nAssociatedJets", "nAssociatedJets", "Number of Jets", "", 100, 0, 100, 0, 0 );
  plots->allJetPt        = result->AddHist1D("allJetPt", "", "Jet p_{T} (all jets) [GeV]", "", 100, 0, 1000);

  // z pt
  plots->binnedZpT = result->AddHist1D( "binnedZpT", "", "z position [mm]", "Sum(p_{T}) [GeV]", 600, -300, 300);

  // Plots only used if track parameters are calculated
  if(calculateTrackParameters){

    // Basic kinematics
    // pt
    plots->pt = result->AddHist1D(
        "pt", "pt", "Track p_{T} [GeV]", "Number of Tracks",
        1000, 0, 500, 0, 1);
    plots->logpt = result->AddHist1D(
        "logpt", "logpt", "Track p_{T} [GeV]", "Number of Tracks",
        1000, 0, 500, 1, 1);
    plots->truthPt = result->AddHist1D(
        "truthPt", "pt", "Truth Track p_{T} [GeV]", "Number of Tracks",
        1000, 0, 1000, 0, 1);
    float ptMinRaw(-200), ptMaxRaw(200); // wont be able to calculate pT resolutions for tracks with > 200 GeV
    plots->ptResRaw = result->AddHist1D(
        "ptResRaw", "pt resolution", "#deltap_{T} [GeV]", "Number of Tracks",
        1000, ptMinRaw, ptMaxRaw, 0, 0);
    plots->ptResRaw_pt = result->AddHist2D(
        "ptResRaw_pt", "pt resolution", "#deltap_{T} [GeV]", "Truth p_{T} [GeV]",
        1000, ptMinRaw, ptMaxRaw, 200, 0, 2000);
    plots->ptResRaw_eta = result->AddHist2D(
        "ptResRaw_eta", "pt resolution", "#deltap_{T} [GeV]", "Truth #eta",
        1000, ptMinRaw, ptMaxRaw, 200, -2, 2);
    plots->ptResRaw_pt_eta = result->AddHist3D(
        "ptResRaw_pt_eta", "pt resolution", "#deltap_{T} [GeV]", "Truth p_{T} [GeV]", "Truth #eta",
        4000, ptMinRaw, ptMaxRaw, 200, 0, 200, 40, -2, 2); // lots of bins needed since pT resolution degreades over many decades

    float ptMin(-10), ptMax(10);
    plots->ptRes = result->AddHist1D(
        "ptRes", "pt resolution", "#deltap_{T}/p_{T}", "Number of Tracks",
        1000, ptMin, ptMax, 0, 0);
    plots->ptRes_pt = result->AddHist2D(
        "ptRes_pt", "pt resolution", "#deltap_{T}/p_{T}", "Truth p_{T} [GeV]",
        1000, ptMin, ptMax, 200, 0, 1000);
    plots->ptRes_eta = result->AddHist2D(
        "ptRes_eta", "eta resolution", "#deltap_{T}/p_{T}", "Truth #eta",
        1000, ptMin, ptMax, 200, -2, 2);
    plots->ptRes_pt_eta = result->AddHist3D(
        "ptRes_pt_eta", "pt resolution", "#deltap_{T}/p_{T}", "Truth p_{T} [GeV]", "Truth #eta",
        1000, ptMin, ptMax, 200, 0, 200, 100, -2, 2);


    //eta
    plots->eta = result->AddHist1D(
        "eta", "eta", "Track #eta", "Number of Tracks",
        100, -4, 4);
    float etaMin(-0.001), etaMax(0.001);
    plots->etaRes = result->AddHist1D(
        "etaRes", "eta resolution", "#delta#eta", "Number of Tracks",
        1000, etaMin, etaMax);
    plots->etaRes_eta = result->AddHist2D(
        "etaRes_eta", "eta resolution", "#delta#eta", "Truth #eta",
        1000, etaMin, etaMax, 100, -2, 2);
    plots->etaRes_pt_eta = result->AddHist3D(
        "etaRes_pt_eta", "pt resolution", "#deltap_{T}/p_{T}", "Truth p_{T} [GeV]", "Truth #eta",
        1000, etaMin, etaMax, 200, 0, 200, 100, -2, 2);

    // z0
    float z0Min(-20), z0Max(20); // large range for low pT tracks needed
    plots->z0Res = result->AddHist1D(
        "z0Res", "z0 resolution", "#deltaz_{0} [mm]", "Number of Tracks",
        1000, z0Min, z0Max);
    plots->z0Res_pt = result->AddHist2D(
        "z0Res_pt", "z0 resolution", "#deltaz_{0} [mm]", "Truth p_{T} [GeV]",
        1000, z0Min, z0Max, 200, 0, 2000);
    plots->z0Res_eta = result->AddHist2D(
        "z0Res_eta", "z0 resolution", "#deltaz_{0} [mm]", "Truth #eta",
        1000, z0Min, z0Max, 200, -2, 2);
    plots->z0Res_pt_eta = result->AddHist3D(
        "z0Res_pt_eta", "z0 resolution", "#deltaz_{0} [mm]", "Truth p_{T} [GeV]", "Truth #eta",
        2000, z0Min, z0Max, 200, 0, 200, 100, -2, 2); // particularly tricky, lots of bins


    // d0
    float d0Min(-50), d0Max(50);
    plots->d0Res = result->AddHist1D(
        "d0Res", "d0 resolution", "#deltad_{0} [mm]", "Number of Tracks",
        1000, d0Min, d0Max);
    plots->d0Res_pt = result->AddHist2D(
        "d0Res_pt", "d0 resolution", "#deltad_{0} [mm]", "Truth p_{T} [GeV]",
        1000, d0Min, d0Max, 200, 0, 2000);
    plots->d0Res_eta = result->AddHist2D(
        "d0Res_eta", "d0 resolution", "#deltad_{0} [mm]", "Truth #eta",
        1000, d0Min, d0Max, 200, -2, 2);
    plots->d0Res_pt_eta = result->AddHist3D(
        "d0Res_pt_eta", "d0 resolution", "#deltad_{0} [mm]", "Truth p_{T} [GeV]", "Truth #eta",
        1000, d0Min, d0Max, 200, 0, 200, 100, -2, 2);

    // cot(theta)
    float CtgThMin(-0.1), CtgThMax(0.1);
    plots->CtgThetaRes = result->AddHist1D(
        "CtgThetaRes", "CtgThetaRes", "#deltacot(#theta)", "Number of Tracks",
        1000, CtgThMin, CtgThMax);
    plots->CtgThetaRes_pt = result->AddHist2D(
        "CtgThetaRes_pt", "CtgThetaRes_pt", "#deltacot(#theta)", "Truth p_{T} [GeV]",
        1000, CtgThMin, CtgThMax, 200, 0, 2000);
    plots->CtgThetaRes_eta = result->AddHist2D(
        "CtgThetaRes_eta", "CtgThetaRes_eta", "#deltacot(#theta)", "Truth #eta",
        1000, CtgThMin, CtgThMax, 200, -2, 2);
    plots->CtgThetaRes_pt_eta = result->AddHist3D(
        "CtgThetaRes_pt_eta", "CtgTheta resolution", "#deltacot(#theta)", "Truth p_{T} [GeV]", "Truth #eta",
        1000, CtgThMin, CtgThMax, 200, 0, 200, 100, -2, 2);

    // phi
    float phiMin(-0.1), phiMax(0.1);
    plots->phiRes = result->AddHist1D(
        "phiRes", "phi resolution", "#delta#phi", "Number of Tracks",
        1000, phiMin, phiMax);
    plots->phiRes_pt = result->AddHist2D(
        "phiRes_pt", "phi resolution", "#delta#phi", "Truth p_{T} [GeV]",
        1000, phiMin, phiMax, 200, 0, 2000);
    plots->phiRes_eta = result->AddHist2D(
        "phiRes_eta", "phi resolution", "#delta#phi", "Truth #eta",
        1000, phiMin, phiMax, 200, -2, 2);
    plots->phiRes_pt_eta = result->AddHist3D(
        "phiRes_pt_eta", "phi resolution", "#delta#phi", "Truth p_{T} [GeV]", "Truth #eta",
        1000, phiMin, phiMax, 200, 0, 200, 100, -2, 2);

    // track
    plots->track_eta_phi_pt = result->AddHist2D(
        "track_eta_phi_pt", "", "#eta", "#phi",
        400, -2, 2, 100, -M_PI, M_PI);
    plots->jet_eta_phi_pt = result->AddHist2D(
        "jet_eta_phi_pt", "", "#eta", "#phi",
        400, -2, 2, 100, -M_PI, M_PI);
  } // end if calculateTrackParameters


}

void calculateResolutions(TClonesArray *branchTruthTrack, TClonesArray *branchTrack, TestPlots *plots, bool DEBUG){

  GenParticle *particle, *truthParticle;
  Track *track, *truthTrack;

  // loop over all truth tracks
  for(Int_t j=0; j<branchTruthTrack->GetEntriesFast(); ++j){
    if(DEBUG) std::cout << "track j " << j << std::endl;
    truthTrack      = (Track*) branchTruthTrack->At(j);
    truthParticle   = (GenParticle*) truthTrack->Particle.GetObject();
    if(DEBUG) std::cout << "Extract particle info" << std::endl;
    if( fabs(truthTrack->Eta) > 2.0 ) continue; // remove tracks with |eta| > 2.0
    if( truthTrack->PT < 1.0) continue;

    // Check is the truth track is from pileup
    if(truthParticle == NULL) continue;

    // Extract the UID of the particle
    UInt_t truthParticleUID = truthParticle->GetUniqueID();

    // Loop over all tracks in event
    for(Int_t i=0; i<branchTrack->GetEntriesFast(); ++i)
    {
      if(DEBUG) std::cout << "track i: " << i << " out of " << branchTrack->GetEntriesFast() << std::endl;
      track = (Track*) branchTrack->At(i);
      particle = (GenParticle*) track->Particle.GetObject();
      if(DEBUG) std::cout << "extracted track and particle" << std::endl;

      // Check is the reco track is from pileup
      if(particle == NULL) continue;

      UInt_t particleUID = particle->GetUniqueID();

      if(DEBUG) std::cout << "truthParticleID: "  << truthParticle->GetUniqueID() << std::endl;
      if(DEBUG) std::cout << "particleID: "  << particle->GetUniqueID() << std::endl;

      // Match the track with the corresponding truth track
      if(truthParticleUID == particleUID){
        if(DEBUG) std::cout << "Matched track to particle" << std::endl;
        double z0Resolution = track->DZ - truthTrack->DZ;
        double ptResolution = track->PT - truthTrack->PT;
        double d0Resolution = track->D0 - truthTrack->D0;
        double etaResolution = track->Eta - truthTrack->Eta;
        double CtgThetaResolution = track->CtgTheta - truthTrack->CtgTheta;
        double phiResolution = track->Phi - truthTrack->Phi;

        float truthPt = truthTrack->PT;
        float truthEta = truthTrack->Eta;

        plots->pt->Fill(track->PT);
        plots->logpt->Fill(track->PT);
        plots->truthPt->Fill(truthTrack->PT);
        plots->ptResRaw->Fill(ptResolution);
        plots->ptResRaw_pt->Fill(ptResolution, truthPt);
        plots->ptResRaw_eta->Fill(ptResolution, truthEta);
        plots->ptResRaw_pt_eta->Fill(ptResolution, truthPt, truthEta);

        plots->ptRes->Fill(ptResolution/truthPt);
        plots->ptRes_pt->Fill(ptResolution/truthPt, truthPt);
        plots->ptRes_eta->Fill(ptResolution/truthPt, truthEta);
        plots->ptRes_pt_eta->Fill(ptResolution/truthPt, truthPt, truthEta);

        plots->eta->Fill(track->Eta);
        plots->etaRes->Fill(etaResolution);
        plots->etaRes_eta->Fill(etaResolution, truthEta);
        plots->etaRes_pt_eta->Fill(etaResolution, truthPt, truthEta);

        plots->z0Res->Fill(z0Resolution);
        plots->z0Res_pt->Fill(z0Resolution, truthPt);
        plots->z0Res_eta->Fill(z0Resolution, truthEta);
        plots->z0Res_pt_eta->Fill(z0Resolution, truthPt, truthEta);

        plots->d0Res->Fill(d0Resolution);
        plots->d0Res_pt->Fill(d0Resolution, truthPt);
        plots->d0Res_eta->Fill(d0Resolution, truthEta);
        plots->d0Res_pt_eta->Fill(d0Resolution, truthPt, truthEta);

        plots->CtgThetaRes->Fill(CtgThetaResolution);
        plots->CtgThetaRes_pt->Fill(CtgThetaResolution, truthPt);
        plots->CtgThetaRes_eta->Fill(CtgThetaResolution, truthEta);
        plots->CtgThetaRes_pt_eta->Fill(CtgThetaResolution, truthPt, truthEta);

        plots->phiRes->Fill(phiResolution);
        plots->phiRes_pt->Fill(phiResolution, truthPt);
        plots->phiRes_eta->Fill(phiResolution, truthEta);
        plots->phiRes_pt_eta->Fill(phiResolution, truthPt, truthEta);
        if(DEBUG) std::cout << "Filled resolution histograms" << std::endl;


        ///////////////////////////////////
        // track trigger?
        ///////////////////////////////////
        // apply a pT cut to the offline tracks, and fill the truth pT histogram (emulate a trigger?)
        for(int i=0; i<19; ++i){
          int pTcut = 5*(i+1);
          if(track->PT > pTcut){
            plots->truthTrackNpT.at(i)->Fill(truthPt);
          }
        }
        plots->truthTrackPt100->Fill(truthPt);
        // end of fill 

      } //end if match track (/particle)
    } // loop over tracks 
  } // loop over tracks
} // end of calculate resolutions

//------------------------------------------------------------------------------

void AnalyseEvents(const int nEvents, ExRootTreeReader *treeReader, TestPlots *plots, bool DEBUG, bool calculateTrackParameters, float minJetPt, float minTrackPt)
{


  // Define branches
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  TClonesArray *branchTruthTrack = treeReader->UseBranch("TruthTrack");
  TClonesArray *branchTrack      = treeReader->UseBranch("Track");
  TClonesArray *branchTrackJet   = treeReader->UseBranch("TrackJet");
  TClonesArray *branchVertex     = treeReader->UseBranch("Vertex");
  TClonesArray *branchPrimaryBin = treeReader->UseBranch("PrimaryBin");
  TClonesArray *branchTower      = treeReader->UseBranch("Tower");
  //TClonesArray *branchPileupParticle = treeReader->UseBranch("PileupParticle");

  // Setup FastJet
  JetDefinition *definition = 0;
  definition = new JetDefinition(antikt_algorithm, 0.4);


  // Loop over all events
  Long64_t allEntries = treeReader->GetEntries();
  int nEventsCorrectlyIdentifiedVertex(0);
  std::cout << "** Chain contains " << allEntries << " events" << std::endl;
  for(Long64_t entry = 0; entry < allEntries; ++entry)
  {

    // limit number of events looped over
    if(nEvents != -1){
      if(entry > nEvents-1) break; // -1 because humans will count from 1 event
    }

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // print every 10% complete
    if( entry % 100==0 ) std::cout << "Event " << entry << " out of " << allEntries << std::endl;




    // Occupancy estimator using Towers
    std::cout << branchTrackJet->GetEntriesFast() << " towers" << std::endl;
    for(int iTower=0; iTower<branchTrackJet->GetEntriesFast(); ++iTower){
      Tower * tower = static_cast<Tower*>(branchTower->At(iTower)); 

      // coordinates of tower center 
      float towerPhi = tower->Phi;
      float towerEta = tower->Eta; 

      // tower edges
      float towerLowEta = tower->Edges[0];
      float towerUpEta = tower->Edges[1];
      float towerLowPhi = tower->Edges[2];
      float towerUpPhi = tower->Edges[3];

      // loop over all particles (tracks) assocated with the tower 
      if(tower->Particles.GetEntriesFast() < 1) continue; 
      std::cout << tower->Particles.GetEntriesFast() << " particles in this tower" << std::endl;

      TRefArray towerParticles = tower->Particles; 
      for(int iParticle=0; iParticle<towerParticles.GetEntriesFast(); ++iParticle){

        std::cout << "before cast" << std::endl;
        GenParticle * towerParticle = dynamic_cast<GenParticle*>(towerParticles.At(iParticle));
        std::cout << "After cast" << std::endl;

        std::cout << "E" << towerParticle->E << std::endl;
        std::cout << "P" << towerParticle->P << std::endl;
        std::cout << towerParticle->GetName() << std::endl;
        std::cout << "After name" << std::endl;
      }
    }







    /////////////////////////////////////////
    // Calculate track parameter resolutions
    /////////////////////////////////////////

    if(calculateTrackParameters){ // turn off if pileup is ON (too many tracks!)
      calculateResolutions(branchTruthTrack, branchTrack, plots, DEBUG);
    }



    ///////////////////////////////////
    // loop over all reco tracks, and find the "primary bin"
    //////////////////////////////////

    // scan parameters
    int beamMinZ    = -300; // mm
    int beamMaxZ    = 300; // mm
    float binWidth  = 1.0; // mm
    float slideStep = 0.1; // mm
    std::map<std::string, float> PBInfo = findPrimaryBin(branchTrack, binWidth, slideStep, beamMinZ, beamMaxZ);
    float zMin   = PBInfo["zMin"];
    float zMax   = PBInfo["zMax"];
    int nTracksInPB = PBInfo["nTracks"];
    float zWidth = PBInfo["zWidth"];
    float PBCentroid = PBInfo["zPosition"]; 


    ///////////////////////////////////
    // Find the location of the real PV
    ///////////////////////////////////
    float vertexZ(0.0);
    int nVerticesInPB(0);
    for(int i=0; i<branchVertex->GetEntriesFast(); ++i){
      Vertex * vertex = (Vertex*) branchVertex->At(i);
      if(vertex->Z < zMax && vertex->Z > zMin) nVerticesInPB++;
      if(vertex->IsPU) continue;
      else{
        vertexZ = vertex->Z;
        //std::cout << "Event: " << entry << " vertex z: " << vertexZ << std::endl;
      }
    }

    //trackMultiplicityInPB->Fill(nTracksInPB);
    //vertexMultiplicityInPB->Fill(nVerticesInPB);
    //nVertexVsZVertexPosition->(nVerticesInPB, vertexZ);
    //nVertexVsZPBPosition->(nVertexVsZPBPosition, PBCentroid);
    //PBvZVertexPosition->(PBCentroid, vertexZ)
    


    // distance between PB (centre) and PV
    //std::cout << "True vertex position: " << vertexZ << " PB  << std::endl;
    //std::cout << "true vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]\t" << PBCentroid << std::endl;
    for(int i=0; i<branchPrimaryBin->GetEntriesFast(); ++i){
      Vertex * primaryBin = dynamic_cast<Vertex*>(branchPrimaryBin->At(i));
      //std::cout << "primary bin center: " << primaryBin->Z << std::endl; 
    }


    // Check to see that real PV is in the PB
    if(vertexZ > zMin && vertexZ < zMax){
      nEventsCorrectlyIdentifiedVertex++;
      //std::cout << "Vertex is inside PB" << std::endl;
      //std::cout << "vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]\t nVertices: " << nVerticesInPB << std::endl;
    }
    else{
      /*******
      std::cout << "Vertex isn't inside PB" << std::endl;
      std::cout << "vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]\t" << fabs(vertexZ - (zMax+zMin)/2) << "\t nVertices: " << nVerticesInPB << std::endl;
      std::cout << "nTracks in PB " << nTracksInPB << std::endl;
      plots->misidPV->Fill( fabs(vertexZ - (zMax+zMin)/2) );
      plots->misidPVLogx->Fill( fabs(vertexZ - (zMax+zMin)/2) );
      //std::cout << "vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]" << std::endl;
      ***********/
    }

    //std::cout << "This event has " << branchTrack->GetEntriesFast() << std::endl;
    int nTracksIterator(0);
    ///////////////////////////////////
    // Select tracks which belong to the PB
    ///////////////////////////////////
    //std::cout << "There are " << branchTrack->GetEntriesFast() << " input tracks" << std::endl;
    std::vector<Track*> tracksAssociatedToPB;
    std::vector<TLorentzVector> vectorsAssociatedToPB;
    for(auto itTrack = branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      float zPosition = dynamic_cast<Track*>(*itTrack)->DZ;
      //std::cout << "\t" << nTracksIterator << "\t" << zPosition << std::endl; 
      if(zPosition > zMin && zPosition < zMax){
        tracksAssociatedToPB.push_back( dynamic_cast<Track*>(*itTrack) );
        vectorsAssociatedToPB.push_back( dynamic_cast<Track*>(*itTrack)->P4());
        //std::cout << "\tmatch" << std::endl;
      }
      nTracksIterator++;
    }
    //std::cout << "and " << nTracksIterator << " iterated tracks" << std::endl;
    //std::cout << "Tracks associated to PB: " << tracksAssociatedToPB.size() << std::endl;

    ///////////////////////////
    // cluster Tracks into Jet
    ///////////////////////////

    // Jets made from tracks associated to the PB
    std::vector<MyClusteredJet> associatedJets = jetCluster(tracksAssociatedToPB, definition, minJetPt, minTrackPt, false);

    // Cluster of truth tracks (test: should be the same as input collection)
    std::vector<MyClusteredJet> nominalTruthJets = jetCluster(branchTruthTrack, definition, minJetPt, minTrackPt, false);

    // Recluster of reco tracks (test: should be the same as the input collection)
    std::vector<MyClusteredJet> nominalJets = jetCluster(branchTrack, definition, minJetPt, minTrackPt, false);

    // Add tracks to the clusted jets
    std::vector<JetWithTracks> assocJetsWithTracks   = addTracksToJet(associatedJets, branchTrack);
    std::vector<JetWithTracks> nominalJetsWithTracks = addTracksToJet(nominalJets, branchTrack);


    // Fill plots with track jet properties (jets directly from Delphes)
    int nJets = branchTrackJet->GetEntriesFast();
    plots->nJets->Fill( nJets );
    for(int i=0; i<nJets && i<7; ++i){
      Jet * jet = (Jet*) branchTrackJet->At(i);
      plots->jetNPt.at(i)->Fill( jet->PT );
      plots->jetNEta.at(i)->Fill( jet->Eta );
      plots->jetNPhi.at(i)->Fill( jet->Phi );

      // Plot properties of tracks inside jet
      int nTracks(0);
      for(int j=0; j<jet->Constituents.GetEntriesFast(); ++j){
        TObject *object = jet->Constituents.At(j);
        if(object->IsA() == Track::Class()){
          Track * track = (Track*) object;
          plots->jetNTrackPt.at(i)->Fill(track->PT);
          nTracks++;
        }
      }
      plots->jetNTrackMulti.at(i)->Fill(nTracks);
    }

    // Fill plots with track jets (clustered from all Delphes tracks, but outside delphes)
    //std::cout << "Filling plots with nominal jets, there are: " << nominalJetsWithTracks.size() << std::endl; 
    for(unsigned int i=0; i<nominalJetsWithTracks.size() && i<7; ++i){
      PseudoJet nominalJet = nominalJetsWithTracks.at(i).first;
      plots->nominalJetNPt.at(i)->Fill( nominalJet.pt() );
      plots->nominalJetNEta.at(i)->Fill( nominalJet.eta() );
      plots->nominalJetNPhi.at(i)->Fill( nominalJet.phi_std() ); // returns phi in the range -pi : pi
      // plots of the tracks inside the jet
      std::vector<Track*> tracks = nominalJetsWithTracks.at(i).second;
      plots->nominalJetNTrackMulti.at(i)->Fill(tracks.size());
      //std::cout << "nominal jet " << i << " has " << tracks.size() << " tracks." << std::endl;
      for(auto track : tracks){
        plots->nominalJetNTrackPt.at(i)->Fill(track->PT);
      }

    }

    // Fill plots with associated jet properties (and properties of the tracks associated to the jets) 
    for(unsigned int i=0; i<assocJetsWithTracks.size() && i<7; ++i){
      PseudoJet associatedJet = assocJetsWithTracks.at(i).first;
      plots->associatedJetNPt.at(i)->Fill( associatedJet.pt() );
      plots->associatedJetNEta.at(i)->Fill( associatedJet.eta() );
      plots->associatedJetNPhi.at(i)->Fill( associatedJet.phi_std() ); // returns phi in the range -pi : pi
      // plots of the tracks inside the jet
      std::vector<Track*> tracks = assocJetsWithTracks.at(i).second;
      plots->associatedJetNTrackMulti.at(i)->Fill(tracks.size());
      for(auto track : tracks){
        plots->associatedJetNTrackPt.at(i)->Fill(track->PT);
      }
    }


    /*********************
    ///////////////////////////////////
    // Loop over all jets and find which tracks are associated to the PB
    // WJF: idea was to find the jets with tracks that came from the PB
    // This was almost all jets, as there is no "directional" information
    ///////////////////////////////////

    std::vector<Jet*> jetsFromPV;
    std::vector<std::pair<int, int>> jetsFromPVProperties ;
    // Loop over jets (or tracks), select tracks with dz inside the window with the largets pT
    for(int i = 0; i < branchTrackJet->GetEntriesFast(); ++i){
      jet = (Jet*) branchTrackJet->At(i);
      int nConstituents = jet->Constituents.GetEntriesFast();
      int nConstituentsInPV(0);

      for(int j = 0; j < nConstituents; ++j){

        TLorentzVector momentum;
        momentum.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        TObject *object = jet->Constituents.At(j);

        if(object == 0) continue; // Check if the constituent is accessible
        if(object->IsA() == Track::Class()){
          track = (Track*) object;
          if(track->DZ > zMin && track->DZ < zMax){
            nConstituentsInPV++;
          }
          //std::cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << std::endl;
        }
        else{
          // add exception handling if this ever happens?
          std::cerr << "Unidentified object in jet collection (not a track): " << object->GetName() << std::endl;
        }
      }
      if(nConstituentsInPV > 0){
        jetsFromPV.push_back( jet );
        jetsFromPVProperties.push_back( std::make_pair(nConstituents, nConstituentsInPV));
      }
    } // end of loop over jets
    // Loop over jets that are associated to the PV
    plots->nAssociatedJets->Fill(jetsFromPV.size() );
    for(int i=0; i<jetsFromPV.size() && i<7; ++i){
      //std::cout << "jet with " << jetsFromPVProperties.at(i).first << " constituents, and " << jetsFromPVProperties.at(i).second << " associated to PV. Fraction: " << static_cast<float>(jetsFromPVProperties.at(i).second)/static_cast<float>(jetsFromPVProperties.at(i).first) << std::endl;
      plots->associatedJetNPt.at(i)->Fill(jetsFromPV.at(i)->PT);

    }
    *******************************/



  } // end loop over entries

  std::cout << "Of " << allEntries << " events, " << nEventsCorrectlyIdentifiedVertex << " had the vertex correctly identified, i.e. " << static_cast<float>(nEventsCorrectlyIdentifiedVertex)/static_cast<float>(allEntries) << " of events."  << std::endl;
} // end AnalyseEvents



//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("pdf");
}


//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  gROOT->SetBatch(1);
  bool DEBUG = false;
  bool calculateTrackParameters = true;
  //bool calculateTrackParameters = false;

  // configurable(>) parameters
  float minJetPt = 5.0; // GeV
  float minTrackPt = 1.0; // GeV

  //TestClass aTest;
  //aTest.hello();

  std::string appName = "trackParameters";
  std::string inputFile = argv[1]; // doesn't complain about cast? Maybe compiler can deal with it :p
  std::string outputFile = argv[2];
  int doPrintHistograms = atoi(argv[3]);
  int nEvents(-1);
  if(argc > 4){
    nEvents = atoi(argv[4]);
    std::cout << "INFO: Will run over " << nEvents << std::endl;
  }
  // 
  if(argc < 4)
  {
    std::cout << " Usage: " << appName << " input_file output_file [optional: nEvents]" << std::endl;
    std::cout << " input_file: ROOT file containing delphes output,"       << std::endl;
    std::cout << " output_file: ROOT file name that will contain output histograms" << std::endl;
    std::cout << " bool: write files to .pdf" << std::endl;
    std::cout << " nEvents: number of events to be processed (for testing, default is all events)." << std::endl;
    return 1;
  }

  // control analysis
  TChain *chain = new TChain("Delphes");
  std::cout << "Adding " << inputFile << " to chain" << std::endl;
  chain->Add(inputFile.c_str());

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;
  BookHistograms(result, plots, calculateTrackParameters);
  AnalyseEvents(nEvents, treeReader, plots, DEBUG, calculateTrackParameters, minJetPt, minTrackPt);

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

    /***********************
    // plots for associated jets (manualy assigned tracks to jets)
    std::cout << "\nEvent "<< entry << std::endl;
    std::cout << "Manual track assignemnt:" << std::endl;
    std::cout << "There are " << associatedJets.size() << " associated jets." << std::endl;
    for(unsigned int i=0; i<associatedJets.size() && i<7; ++i){
      PseudoJet associatedJet = associatedJets.at(i).first;
      plots->associatedJetNPt.at(i)->Fill( associatedJet.pt() );
      plots->associatedJetNEta.at(i)->Fill( associatedJet.eta() );
      plots->associatedJetNPhi.at(i)->Fill( associatedJet.phi_std() ); // returns phi in the range -pi : pi
      // use unique ID to match the track
      std::vector<float> associatedTrackIDs = associatedJets.at(i).second;
      int nTracks(0);
      for(auto itTrack=branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
        Track* track = dynamic_cast<Track*>(*itTrack);
        for(auto trackID : associatedTrackIDs){
          if(track->GetUniqueID() == trackID){
            nTracks++;
            plots->associatedJetNTrackPt.at(i)->Fill(track->PT);
            std::cout << "\t\t" << track->PT << std::endl;
          }
        }
      }
      std::cout << "jet " << i << " had " << nTracks << " tracks." << std::endl;
      plots->associatedJetNTrackMulti.at(i)->Fill(nTracks);
    }
    plots->nAssociatedJets->Fill( associatedJets.size() );
    *******************/
