#define trackLoop_cxx

#include "trackLoop.h"
#include "ExRootAnalysis/ExRootResult.h"

#include <algorithm>
#include <set>

class cutHolder{
  private:
    std::string m_cutName;
  public:

    int Fake;
    int True;

    cutHolder(std::string cutName){
      m_cutName = cutName;
      Fake =0;
      True =0;
    }

};

inline TVector3 fillTV3(float r, float phi, float z){

  // calculate cartesian coordinates
  //float theta = 2*atan( exp(-eta) );
  //float z = r*sin(theta); 
  float x = r*cos(phi);
  float y = r*sin(phi);

  TVector3 tempVec(x, y, z);
  return tempVec; 
}

template<typename T>
inline bool isnt_in_vector(std::vector<T>& theVector, T& item){
  return (std::find(theVector.begin(), theVector.end(), item) == theVector.end());
}

inline float calcDeltaPhiLimit(float r1, float r2){
  float twiceRadius = 1000 * 2 / 1.99; // bending radius [mm] for track of 2 GeV in 4 T magnetic field
  twiceRadius *= 2; // twice the radius  
  return fabs(acos(r1 / twiceRadius) - acos(r2 / twiceRadius));
}

inline void printTrackStats(TString info, float n1, float n2){
  std::cout << info << n1 << " / " << n2 << " \t(" << 100*n1/n2 << "%)" << std::endl;
}


void trackLoop::Loop(Plots* plots, int branchCounter, std::vector<float> tripletLayers, float zresiduumCut)
{

  // Calculate some DeltaPhi angles for 
  float deltaPhi12Limit = calcDeltaPhiLimit(tripletLayers.at(0), tripletLayers.at(1));
  float deltaPhi23Limit = calcDeltaPhiLimit(tripletLayers.at(1), tripletLayers.at(2));
  float deltaPhi13Limit = calcDeltaPhiLimit(tripletLayers.at(0), tripletLayers.at(2));

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "About to loop over " << nentries << std::endl;
  Long64_t nbytes = 0, nb = 0;

  /*********
  std::map<std::string, int> cutMapPt2;
  cutHolder("kappa");
  cutMapPt2["kappa"] = 0;
  cutMapPt2["deltaPhi12"] = 0;
  cutMapPt2["deltaPhi21"] = 0;
  cutMapPt2["deltaPhi13"] = 0;
  cutMapPt2["deltaPhi"] = 0;
  cutMapPt2["zresiduum"] = 0;
  **********/

  // All tracks
  int numRecoTracksPt0 = nentries;
  int numRecoTracksPt2(0);
  int numRecoTracksPt0Surviving(0);
  int numRecoTracksPt2Surviving(0);

  // True tracks
  int numTrueTracksPt0(0);
  int numTrueTracksPt2(0);
  int numTrueTracksPt0Surviving(0);
  int numTrueTracksPt2Surviving(0);

  // Fake tracks 
  int numFakeTracksPt0(0);
  int numFakeTracksPt2(0);
  int numFakeTracksPt0Surviving(0);
  int numFakeTracksPt2Surviving(0);

  std::vector<float> hitMomenta;
  hitMomenta.reserve(nentries); // is this too much?  

  std::vector<int> hitMomentaUniqueID;
  std::set<long long int> uniqueHitMomenta; 

  //std::vector<float> trueTrackHitPt;
  //std::vector<float> fakeTrackHitPt; 
  std::set<float> trueTrackHitPt;
  std::set<float> fakeTrackHitPt; 

  // Some definitions for the kappa cut
  //std::vector<float> pTthresholds    = {0.0,   2.0,   10.0,  50.0};
  //std::vector<float> kappaThresholds = {0.005, 0.004, 0.002, 0.001}; 
  std::vector<float> pTthresholds    = {0.0,   2.0,   10.0,   20.0,  30.0,  40.0,  50.0};
  std::vector<float> kappaThresholds = {0.005, 0.002, 0.0015, 0.001, 0.001, 0.001, 0.0005}; 


  // Loop over tracks 
  // Note that there are more tracks than hits in the outer layer, since there are fake tracks 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%50000 == 0) std::cout << "Entry " << jentry << "/" << nentries << std::endl;

    bool keepTrack = true; 

    // Useful delta kappa 
    float deltaKappa = kappa_123 - kappa_013;

    // Calculate sign of angular differences between points (in 3D)
    TVector3 p1 = fillTV3(hit1rho, hit1phi, 0); 
    TVector3 p2 = fillTV3(hit2rho, hit2phi, 0); 
    TVector3 p3 = fillTV3(hit3rho, hit3phi, 0); 
    TVector3 zhad = fillTV3(0, 0, 1);

    float sin12 = p1.Cross(p2) / (p1.Mag() * p2.Mag()) // does this even make any sense, as to calculate the signed angle difference 


    ////////////////////////////////////////////////
    // Counters for numbers of tracks (before cuts)
    ////////////////////////////////////////////////
    
    // All tracks (true or fake) 
    if(pT > 2){
      numRecoTracksPt2++;
    }
    if(isFake){
      numFakeTracksPt0++;
      if(pT > 2) numFakeTracksPt2++;
    }
    else{
      numTrueTracksPt0++;
      if(pT > 2) numTrueTracksPt2++;
    }

    // Plots for tracks without selection
    plots->nHitsPt.at(branchCounter)->Fill(hit3pT);

    /////////////////////////
    // Ideas to cut out fakes
    /////////////////////////
    
    // Kappa 
    for(int i=0; i<pTthresholds.size(); ++i){
      if(pT > pTthresholds.at(i)){
        if(fabs(deltaKappa) > kappaThresholds.at(i)) continue; // cut the track if kappa is too large
      }
    }


    // Delta-phi constraints
    float deltaPhi12 = fabs(hit1phi - hit2phi);
    float deltaPhi23 = fabs(hit2phi - hit3phi);
    float deltaPhi13 = fabs(hit1phi - hit3phi);
    if(deltaPhi12 > deltaPhi12Limit) continue;
    if(deltaPhi23 > deltaPhi23Limit) continue;
    if(deltaPhi13 > deltaPhi13Limit) continue;
  

    // z residuum 
    if(fabs(zresiduum) > zresiduumCut) continue; 

    // z residuum over eta (extra information)
    if(fabs( zresiduum/fabs(eta) ) > 0.4) continue; 

    //long long int hit3PtID = static_cast<long long int>( hit3pT*10E8 );
    //uniqueHitMomenta.insert(hit3PtID);

    /*********
      if( isnt_in_vector( hitMomenta, hit3pT ) ){ 
      hitMomenta.push_back(hit3pT); 
      uniqueHitMomenta++;
    // If so, add to hit pT 
    }
     *********/

    plots->recoTrackPt.at(branchCounter)->Fill(pT);

    // Count number of surviving tracks 
    numRecoTracksPt0Surviving++;
    if(pT > 2) {
      numRecoTracksPt2Surviving++;
    }


    if(isFake){
      // FAKE TRACKS 
      plots->recoTrackPt_fake.at(branchCounter)->Fill(pT); 
      //fakeTrackHitPt.push_back(hit3pT);
      fakeTrackHitPt.insert(hit3pT);
      numFakeTracksPt0Surviving++;
      if(pT>2) numFakeTracksPt2Surviving++;
    }
    else{
      // TRUE TRACKS
      //trueTrackHitPt.push_back(hit3pT);
      trueTrackHitPt.insert(hit3pT);
      plots->recoTrackHitPt_true.at(branchCounter)->Fill(hit3pT);
      plots->recoTrackHitEta_true.at(branchCounter)->Fill(hit3eta);
      plots->recoTrackPt_true.at(branchCounter)->Fill(pT); 
      numTrueTracksPt0Surviving++;
      if(pT > 2){
        numTrueTracksPt2Surviving++;
      }
    }

  } // End loop over tracks

  std::cout << "Results: " << std::endl;

  // All tracks
  printTrackStats("Number of tracks: ", numRecoTracksPt0Surviving, numRecoTracksPt0);
  printTrackStats("Number of tracks pT > 2: ", numRecoTracksPt2Surviving, numRecoTracksPt2); 

  // Fake tracks
  printTrackStats("Fake tracks: ", numFakeTracksPt0Surviving, numFakeTracksPt0);
  printTrackStats("Fake tracks pT > 2: ", numFakeTracksPt2Surviving, numFakeTracksPt2);

  // True tracks
  printTrackStats("True tracks: ", numTrueTracksPt0Surviving, numTrueTracksPt0);
  printTrackStats("True tracks pT > 2: ", numTrueTracksPt2Surviving, numTrueTracksPt2);

  std::cout << "Average fake rate (pT > 2): " << 100*float(numFakeTracksPt2Surviving)/float(numRecoTracksPt2Surviving) << " %" << std::endl;

  std::cout << "Unique hits momenta: "   << uniqueHitMomenta.size() << std::endl;


  // Fill histograms (these probably wont work, since they should be event quantities)  
  plots->nDelphesHits.at(branchCounter)->Fill(numRecoTracksPt0);
  plots->nDelphesHitsPt2.at(branchCounter)->Fill(numRecoTracksPt2);
  plots->nRecoTracksMatched.at(branchCounter)->Fill(numTrueTracksPt0);
  plots->nRecoTracksMatchedPt2.at(branchCounter)->Fill(numTrueTracksPt2);
  plots->nRecoTracksPt2.at(branchCounter)->Fill(numRecoTracksPt2);

  std::cout << "Finished Loop()" << std::endl;
}

# ifndef __CINT__
int main(int argc, char* argv[]){

  TString sampleName = argv[1]; 
  TString outputFile = argv[2];

  ExRootResult *result = new ExRootResult();
  Plots        *plots  = new Plots;

  BookHistograms(result, plots);

  std::vector<TString> branchNames = { "Tracks10", "Tracks20", "Tracks30", "Tracks40", "Tracks50" };

  std::map<TString, std::vector<float>> tripletLayouts;
  tripletLayouts["Tracks10"] = {572, 582, 592}; // 10 mm
  tripletLayouts["Tracks20"] = {562, 582, 602}; // 20 mm
  tripletLayouts["Tracks30"] = {552, 582, 612}; // 30 mm
  tripletLayouts["Tracks40"] = {542, 582, 622}; // 40 mm
  tripletLayouts["Tracks50"] = {532, 582, 632}; // 50 mm

  std::map<TString, float> zresduumCuts;
  zresduumCuts["Tracks10"] = 0.05; // 10 mm
  zresduumCuts["Tracks20"] = 0.075; // 20 mm
  zresduumCuts["Tracks30"] = 0.1; // 30 mm
  zresduumCuts["Tracks40"] = 0.2; // 40 mm
  zresduumCuts["Tracks50"] = 0.3; // 50 mm

  std::map<TString, float> zresiduumOverEtaCut;
  zresduumCuts["Tracks10"] = 0.1; // 10 mm
  zresduumCuts["Tracks20"] = 0.4; // 20 mm
  zresduumCuts["Tracks30"] =0.4; // 30 mm
  zresduumCuts["Tracks40"] =0.4; // 40 mm
  zresduumCuts["Tracks50"] = 0.4; // 50 mm


  int branchCounter(0);
  for(auto branchName : branchNames){
    trackLoop theLoop(branchName, sampleName);
    theLoop.Loop(plots, branchCounter, tripletLayouts[branchName], zresduumCuts[branchName]);
    branchCounter++;
  }

  // write  histograms 
  std::cout << "Writing to file: " << outputFile << std::endl;
  result->Write(outputFile);

  // cleanup
  delete plots;
  delete result;

  return 0;
}
# endif
