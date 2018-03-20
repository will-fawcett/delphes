#include "ExRootAnalysis/ExRootResult.h"
#include "TH1.h"

struct Plots
{


  // Number of true tracks
  TH1* nDelphesTracks;
  TH1* nDelphesTracks1GeV;
  TH1* nDelphesTracks10GeV;
  std::vector<TH1*> nDelphesHits; 
  std::vector<TH1*> nDelphesHitsPt2; 
  
  // number of reconstructed tracks
  std::vector<TH1*> nRecoTracks; 
  std::vector<TH1*> nRecoTracksPt2; 
  std::vector<TH1*> nRecoTracksMatched;
  std::vector<TH1*> nRecoTracksMatchedPt2;
  TH1* fractionOfFakeTracks;

  TH1* trueParticlePt_numParticles;
  TH1* trueParticlePt_numRecoTracks;


  // Fake rate as a function of pT, eta
  std::vector<TH1*> recoTrackPt;
  std::vector<TH1*> recoTrackEta;

  std::vector<TH1*> recoTrackPt_fake;
  std::vector<TH1*> recoTrackEta_fake; 

  std::vector<TH1*> recoTrackPt_true;
  std::vector<TH1*> recoTrackEta_true; 

  std::vector<TH1*> recoTrackHitPt_true;
  std::vector<TH1*> recoTrackHitEta_true;
  std::vector<TH1*> recoTrackHitEtaPt2_true;
  
  std::vector<TH1*> nHitsPt; 
  std::vector<TH1*> nHitsEta;  
  std::vector<TH1*> nHitsEtaPt2;  
  std::vector<TH1*> recoTrackPtResolution;

  // Simple histogram of curvature difference 
  TH1* curvatureDifference_true;
  TH1* curvatureDifference_fake;

  TGraph* curvatureDifference_curvature_true;
  TGraph* curvatureDifference_curvature_fake;

  // plots to show the magnitude of the sign difference between angles
  std::vector<TH1*> z_phi12_true;
  std::vector<TH1*> z_phi13_true;
  std::vector<TH1*> z_phi23_true;
  std::vector<TH1*> z_phi12m23_true;
  std::vector<TH1*> z_phi12p23_true;

  std::vector<TH1*> z_phi12_fake;
  std::vector<TH1*> z_phi13_fake;
  std::vector<TH1*> z_phi23_fake;
  std::vector<TH1*> z_phi12m23_fake;
  std::vector<TH1*> z_phi12p23_fake;


};

void BookHistograms(ExRootResult *result, Plots *plots)
{

  int hitMultiplicity = 12000;
  int hitMultiplicityBins = hitMultiplicity/10; 

  for(int i=0; i<5; ++i ){
    std::string trackerID = std::to_string((i+1)*10);

    plots->nDelphesHits.push_back(
        result->AddHist1D("nDelphesHits_"+trackerID, "Number of hits in outermost layer", "", "", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );
    plots->nDelphesHitsPt2.push_back(
        result->AddHist1D("nDelphesHitsPt2_"+trackerID, "Number of hits in outermost layer with pT>2", "", "", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );
    plots->nRecoTracks.push_back(
        result->AddHist1D("nRecoTracks_"+trackerID, "Tracks+fakes", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );
    plots->nRecoTracksMatched.push_back(
        result->AddHist1D("nRecoTracksMatched_"+trackerID, "Matched tracks", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );
    plots->nRecoTracksPt2.push_back(
        result->AddHist1D("nRecoTracksPt2_"+trackerID, "Tracks+fakes", "Number of tracks pT>2", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );
    plots->nRecoTracksMatchedPt2.push_back(
        result->AddHist1D("nRecoTracksMatchedPt2_"+trackerID, "Matched tracks pT>2GeV", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0)
        );

    int nBinsPt = 1000;
    int pTMax = 500;
    plots->recoTrackPt.push_back(
        result->AddHist1D("recoTrackPt_"+trackerID, "Reco track pT", "", "", nBinsPt, 0, pTMax, 0, 0)
        );
    plots->recoTrackPt_fake.push_back(
        result->AddHist1D("recoTrackPt_fake_"+trackerID, "Fake Reco track pT", "", "", nBinsPt, 0, pTMax, 0, 0)
        );
    plots->recoTrackPt_true.push_back(
        result->AddHist1D("recoTrackPt_true_"+trackerID, "True Reco track pT", "", "", nBinsPt, 0, pTMax, 0, 0)
        );
    plots->recoTrackHitPt_true.push_back(
        result->AddHist1D("recoTrackHitPt_true_"+trackerID, "True Reco track, pT of hit", "", "", nBinsPt, 0, pTMax, 0, 0)
        );
    plots->recoTrackHitEta_true.push_back(
        result->AddHist1D("recoTrackHitEta_true_"+trackerID, "True Reco track, Eta of hit", "", "", 100, -5, 5, 0, 0)
        );
    plots->recoTrackHitEtaPt2_true.push_back(
        result->AddHist1D("recoTrackHitEtaPt2_true_"+trackerID, "True Reco track, Eta of hit, hit pT > 2 GeV", "", "", 100, -5, 5, 0, 0)
        );
    plots->nHitsPt.push_back(
        result->AddHist1D("nHitsPt_"+trackerID, "Hit Pt", "", "", nBinsPt, 0, pTMax, 0, 0)
        );
    plots->nHitsEta.push_back(
        result->AddHist1D("nHitsEta_"+trackerID, "Hit Eta", "", "", 100, -5, 5, 0, 0)
        );
    plots->nHitsEtaPt2.push_back(
        result->AddHist1D("nHitsEtaPt2_"+trackerID, "Hit Eta, hit pT > 2 GeV", "", "", 100, -5, 5, 0, 0)
        );

    plots->recoTrackEta.push_back(
        result->AddHist1D("recoTrackEta_"+trackerID, "Reco track eta", "", "", 100, -5, 5, 0, 0)
        );
    plots->recoTrackEta_fake.push_back(
        result->AddHist1D("recoTrackEta_fake_"+trackerID, "Fake Reco track eta", "", "", 100, -5, 5, 0, 0)
        );
    plots->recoTrackEta_true.push_back(
        result->AddHist1D("recoTrackEta_true_"+trackerID, "true Reco track eta", "", "", 100, -5, 5, 0, 0)
        );

    plots->recoTrackPtResolution.push_back(
        result->AddHist1D("recoTrackPtResolution_"+trackerID, "Reconstructed track pT resolution", "", "", 400, -10, 10, 0, 0)
        );

    plots->z_phi12_true.push_back(
        result->AddHist1D("z_phi12_true_"+trackerID, "signed z component of vector of cross between phi1 and phi2", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi13_true.push_back(
        result->AddHist1D("z_phi13_true_"+trackerID, "signed z component of vector of cross between phi1 and phi3", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi23_true.push_back(
        result->AddHist1D("z_phi23_true_"+trackerID, "signed z component of vector of cross between phi2 and phi3", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi12m23_true.push_back(
        result->AddHist1D("z_phi12m23_true"+trackerID, "", "", "", 100, -0.001, 0.001)
        );
    plots->z_phi12p23_true.push_back(
        result->AddHist1D("z_phi12p23_true"+trackerID, "", "", "", 100, -0.2, 0.2)
        );
    plots->z_phi12_fake.push_back(
        result->AddHist1D("z_phi12_fake_"+trackerID, "signed z component of vector of cross between phi1 and phi2", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi13_fake.push_back(
        result->AddHist1D("z_phi13_fake_"+trackerID, "signed z component of vector of cross between phi1 and phi3", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi23_fake.push_back(
        result->AddHist1D("z_phi23_fake_"+trackerID, "signed z component of vector of cross between phi2 and phi3", "", "", 100, -0.2, 0.2, 0, 0)
        );
    plots->z_phi12m23_fake.push_back(
        result->AddHist1D("z_phi12m23_fake"+trackerID, "", "", "", 100, -0.001, 0.001)
        );
    plots->z_phi12p23_fake.push_back(
        result->AddHist1D("z_phi12p23_fake"+trackerID, "", "", "", 100, -0.2, 0.2)
        );




  }

  plots->nDelphesTracks = result->AddHist1D("nDelphesTracks", "TrueTracks", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0);
  plots->nDelphesTracks1GeV = result->AddHist1D("nDelphesTracks1GeV", "TrueTracks > 1 GeV", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0);
  plots->nDelphesTracks10GeV = result->AddHist1D("nDelphesTracks10GeV", "TrueTracks > 10 GeV", "Number of tracks", "Number of events", hitMultiplicityBins, 0, hitMultiplicity,0,0);

  plots->curvatureDifference_true = result->AddHist1D("curvatureDifference_true", "", "Curvature difference", "Number of tracks", 100, -0.1, 0.1, 0, 0);
  plots->curvatureDifference_fake = result->AddHist1D("curvatureDifference_fake", "", "Curvature difference", "Number of tracks", 100, -0.1, 0.1, 0, 0);

  //plots->curvatureDifference_curvature_true = new TGraph();
  plots->curvatureDifference_curvature_true = result->AddTGraph("curvatureDifference_curvature_true");
  plots->curvatureDifference_curvature_fake = result->AddTGraph("curvatureDifference_curvature_fake");


  plots->fractionOfFakeTracks = result->AddHist1D("fractionOfFakeTracks", "", "", "", 100, 0, 1, 0, 0);

  plots->trueParticlePt_numParticles = result->AddHist1D("trueParticlePt_numParticles", "", "", "", 100, 0, 100, 0, 0);
  plots->trueParticlePt_numRecoTracks = result->AddHist1D("trueParticlePt_numRecoTracks", "", "", "", 100, 0, 100, 0, 0);

}

