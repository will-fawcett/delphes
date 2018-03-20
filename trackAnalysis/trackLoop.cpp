#define trackLoop_cxx

#include "TVector3.h"

#include "trackLoop.h"
#include "ExRootAnalysis/ExRootResult.h"

#include <algorithm>
#include <set>

#include "IClassifierReader.h"
#include "BDTG_depth3_10mm.h"
#include "BDTG_depth3_20mm.h"
#include "BDTG_depth3_30mm.h"
#include "BDTG_depth3_40mm.h"
#include "BDTG_depth3_50mm.h"

#include "json.h"
#include <fstream>


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


json trackLoop::Loop(Plots* plots, int branchCounter, std::vector<float> tripletLayers, float zresiduumCut, int cutsOption, int nEventsMax, float BDTCut, float deltaKappaCut)
{
  // triplet spacing string
  std::string tripletSpacingStr = std::to_string( (branchCounter+1)*10 );

  // Calculate some DeltaPhi angles for 
  float deltaPhi12Limit = calcDeltaPhiLimit(tripletLayers.at(0), tripletLayers.at(1));
  float deltaPhi23Limit = calcDeltaPhiLimit(tripletLayers.at(1), tripletLayers.at(2));
  float deltaPhi13Limit = calcDeltaPhiLimit(tripletLayers.at(0), tripletLayers.at(2));

  if (fChain == 0){
    json empty; // hack
    return empty; 
  }

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
  int numRecoTracksHitPt2(0);
  int numRecoTracksPt3(0);
  int numRecoTracksPt5(0);

  int numRecoTracksPt0Surviving(0);
  int numRecoTracksPt2Surviving(0);
  int numRecoTracksHitPt2Surviving(0);
  int numRecoTracksPt3Surviving(0);
  int numRecoTracksPt5Surviving(0);

  // True tracks
  int numTrueTracksPt0(0);
  int numTrueTracksPt2(0);
  int numTrueTracksHitPt2(0);
  int numTrueTracksPt3(0);
  int numTrueTracksPt5(0);
  int numTrueTracksPt0Surviving(0);
  int numTrueTracksPt2Surviving(0);
  int numTrueTracksHitPt2Surviving(0);
  int numTrueTracksPt3Surviving(0);
  int numTrueTracksPt5Surviving(0);

  // Fake tracks 
  int numFakeTracksPt0(0);
  int numFakeTracksPt2(0);
  int numFakeTracksHitPt2(0);
  int numFakeTracksPt3(0);
  int numFakeTracksPt5(0);
  int numFakeTracksPt0Surviving(0);
  int numFakeTracksPt2Surviving(0);
  int numFakeTracksHitPt2Surviving(0);
  int numFakeTracksPt3Surviving(0);
  int numFakeTracksPt5Surviving(0);

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


  // Definition for the BDT
  std::vector<std::string> inputNames = { 
    "abs(tracks"+tripletSpacingStr+".kappa_123-tracks"+tripletSpacingStr+".kappa_013)",
    "tracks"+tripletSpacingStr+".pT",
    "abs(tracks"+tripletSpacingStr+".zresiduum)",
    "abs(tracks"+tripletSpacingStr+".beamlineIntersect)",
    "abs(tracks"+tripletSpacingStr+".hit1phi-tracks"+tripletSpacingStr+".hit2phi)",
    "abs(tracks"+tripletSpacingStr+".hit2phi-tracks"+tripletSpacingStr+".hit3phi)",
    "abs(tracks"+tripletSpacingStr+".hit1phi-tracks"+tripletSpacingStr+".hit3phi)",
    "tracks"+tripletSpacingStr+".z_phi12*tracks"+tripletSpacingStr+".z_phi23",
    "abs(tracks"+tripletSpacingStr+".zresiduum/tracks"+tripletSpacingStr+".eta)"
  };

  IClassifierReader * BDTCalculator;
  if(branchCounter == 0 ) BDTCalculator  = new ReadBDTG10( inputNames ); 
  if(branchCounter == 1 ) BDTCalculator  = new ReadBDTG20( inputNames ); 
  if(branchCounter == 2 ) BDTCalculator  = new ReadBDTG30( inputNames ); 
  if(branchCounter == 3 ) BDTCalculator  = new ReadBDTG40( inputNames ); 
  if(branchCounter == 4 ) BDTCalculator  = new ReadBDTG50( inputNames ); 


  // Loop over tracks 
  // Note that there are more tracks than hits in the outer layer, since there are fake tracks 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%50000 == 0) std::cout << "Entry " << jentry << "/" << nentries << std::endl;

    if(nEventsMax > 0){
      if(jentry > nEventsMax) break; // cancel loop after nEventsMax
    }


    // Useful delta kappa 
    float deltaKappa = kappa_123 - kappa_013;

    // Calculate sign of angular differences between points (in 3D)
    TVector3 p1 = fillTV3(hit1rho, hit1phi, 0); 
    TVector3 p2 = fillTV3(hit2rho, hit2phi, 0); 
    TVector3 p3 = fillTV3(hit3rho, hit3phi, 0); 
    TVector3 zhat(0, 0, 1); // z direction vector

    // z component of the cross product of two angles
    float z_phi12 = p1.Cross(p2)*zhat / (p1.Mag() * p2.Mag()); 
    float z_phi13 = p1.Cross(p3)*zhat / (p1.Mag() * p3.Mag());
    float z_phi23 = p2.Cross(p3)*zhat / (p2.Mag() * p3.Mag());
    float z_12m23 = z_phi12 * z_phi23;
    float z_12p23 = z_phi12 + z_phi13; 

    if(isFake){
      fakeTrackHitPt.insert(hit3pT);
      plots->z_phi12_fake.at(branchCounter)->Fill(z_phi12);
      plots->z_phi13_fake.at(branchCounter)->Fill(z_phi13);
      plots->z_phi23_fake.at(branchCounter)->Fill(z_phi23);
      plots->z_phi12m23_fake.at(branchCounter)->Fill(z_12m23);
      plots->z_phi12p23_fake.at(branchCounter)->Fill(z_12p23);
    }
    else{
      trueTrackHitPt.insert(hit3pT);
      plots->z_phi12_true.at(branchCounter)->Fill(z_phi12);
      plots->z_phi13_true.at(branchCounter)->Fill(z_phi13);
      plots->z_phi23_true.at(branchCounter)->Fill(z_phi23);
      plots->z_phi12m23_true.at(branchCounter)->Fill(z_12m23);
      plots->z_phi12p23_true.at(branchCounter)->Fill(z_12p23);
    }


    ////////////////////////////////////////////////
    // Counters for numbers of tracks (before cuts)
    ////////////////////////////////////////////////
    
    // All tracks (true or fake) 
    if(hit3pT > 2) numRecoTracksHitPt2++;
    if(pT > 2){
      numRecoTracksPt2++;
      if(pT > 3) numRecoTracksPt3++;
      if(pT > 5) numRecoTracksPt5++;
    }
    if(isFake){
      numFakeTracksPt0++;
      if(hit3pT > 2) numFakeTracksHitPt2++;
      if(pT > 2) numFakeTracksPt2++;
      if(pT > 3) numFakeTracksPt3++;
      if(pT > 5) numFakeTracksPt5++;
    }
    else{
      numTrueTracksPt0++;
      // plot for efficiency (assume selection 100% efficient for true tracks, prior to the selections in this code) 
      plots->nHitsPt.at(branchCounter)->Fill(hit3pT);
      if(hit3pT > 2) numTrueTracksHitPt2++;
      if(pT > 2) numTrueTracksPt2++;
      if(pT > 3) numTrueTracksPt3++;
      if(pT > 5) numTrueTracksPt5++;
    }



    /////////////////////////
    // Use of BDTg to remove fakes
    /////////////////////////
    if(cutsOption == 1 || cutsOption == 2){
      std::vector<double> inputVariables = {
        fabs(kappa_123-kappa_013),
        pT,
        fabs(zresiduum),
        fabs(beamlineIntersect),
        fabs(hit1phi-hit2phi),
        fabs(hit2phi-hit3phi),
        fabs(hit1phi-hit3phi),
        z_phi12*z_phi23,
        fabs(zresiduum/eta)
      };

      // get the BDT response  
      double BDTResponse = BDTCalculator->GetMvaValue( inputVariables );
      if(BDTResponse < BDTCut) continue; 
    }
    if(cutsOption == 0 || cutsOption == 2){

      /////////////////////////
      // Ideas to cut out fakes
      /////////////////////////
      
      // Kappa 
      bool rejectTrack(false);
      for(int i=0; i<pTthresholds.size(); ++i){
        if(pT > pTthresholds.at(i)){
          if(fabs(deltaKappa) > kappaThresholds.at(i)) rejectTrack = true; // cut the track if kappa is too large
        }
      }
      if(rejectTrack) continue; 
      //if(pT > 0 && fabs(deltaKappa)    > 0.005) continue;
      if(pT > 0 && fabs(deltaKappa)    > 0.003) continue;
      if(pT > 2.0 && fabs(deltaKappa)  > 0.002) continue; // 0.002 is super efficient!
      if(pT > 2.0 && fabs(deltaKappa)  > 0.0015) continue; // new, not harsh  

      if(pT > 3.5 && fabs(deltaKappa)  > 0.001) continue; // new 

      //if(pT > 5.0 && fabs(deltaKappa) > 0.001) continue; // new 
      //if(pT > 2.0 && fabs(deltaKappa)  > 0.0005) continue;
      //if(pT > 10.0 && fabs(deltaKappa) > 0.0015) continue;
      //if(pT > 10.0 && fabs(deltaKappa) > 0.0005) continue;
      //if(pT > 20.0 && fabs(deltaKappa) > 0.001) continue;
      if(pT > 50.0 && fabs(deltaKappa) > 0.0005) continue;

      // Apply scanned deltaKappa cut
      if(fabs(deltaKappa) > deltaKappaCut) continue; 



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

      // Phi constraint, such that the differences in phi angles are in the same sense when moving from one tracking layer to the next
      if(z_12m23 < 0) continue;

      // beamline intersect
      if(fabs(beamlineIntersect) > 200) continue; 

    } // end if else(cutsOption)


    //long long int hit3PtID = static_cast<long long int>( hit3pT*10E8 );

  
    plots->recoTrackPt.at(branchCounter)->Fill(pT);

    // Count number of surviving tracks 
    numRecoTracksPt0Surviving++;
    if(hit3pT > 2) numRecoTracksHitPt2Surviving++;
    if(pT > 2) numRecoTracksPt2Surviving++;
    if(pT > 3) numRecoTracksPt3Surviving++;
    if(pT > 5) numRecoTracksPt5Surviving++;

    if(isFake){
      // FAKE TRACKS 
      plots->recoTrackPt_fake.at(branchCounter)->Fill(pT); 
      //fakeTrackHitPt.push_back(hit3pT);
      numFakeTracksPt0Surviving++;
      if(hit3pT > 2) numFakeTracksHitPt2Surviving++;
      if(pT>2) numFakeTracksPt2Surviving++;
      if(pT>3) numFakeTracksPt3Surviving++;
      if(pT>5) numFakeTracksPt5Surviving++;
    }
    else{
      // TRUE TRACKS
      //trueTrackHitPt.push_back(hit3pT);
      plots->recoTrackPt_true.at(branchCounter)->Fill(pT); 
      numTrueTracksPt0Surviving++;
      if(hit3pT > 2) numTrueTracksHitPt2Surviving++;
      if(pT > 2)numTrueTracksPt2Surviving++;
      if(pT > 3)numTrueTracksPt3Surviving++;
      if(pT > 5)numTrueTracksPt5Surviving++;
      plots->recoTrackHitPt_true.at(branchCounter)->Fill(hit3pT); // for efficiency
      plots->recoTrackHitEta_true.at(branchCounter)->Fill(hit3eta); // for efficiency
    }

  } // End loop over tracks

  std::cout << "Results: " << std::endl;
  std::cout << "unique fake hit pT: " << fakeTrackHitPt.size() << std::endl;
  std::cout << "unique true hit pT: " << trueTrackHitPt.size() << std::endl;

  // All tracks
  printTrackStats("Number of tracks: ", numRecoTracksPt0Surviving, numRecoTracksPt0);
  printTrackStats("Number of tracks pT > 2: ", numRecoTracksPt2Surviving, numRecoTracksPt2); 
  printTrackStats("Number of tracks hit pT > 2: ", numRecoTracksHitPt2Surviving, numRecoTracksHitPt2); 
  printTrackStats("Number of tracks pT > 3: ", numRecoTracksPt3Surviving, numRecoTracksPt3); 
  printTrackStats("Number of tracks pT > 5: ", numRecoTracksPt5Surviving, numRecoTracksPt5); 

  // Fake tracks
  printTrackStats("Fake tracks: ", numFakeTracksPt0Surviving, numFakeTracksPt0);
  printTrackStats("Fake tracks pT > 2: ", numFakeTracksPt2Surviving, numFakeTracksPt2);
  printTrackStats("Fake tracks hit pT > 2: ", numFakeTracksHitPt2Surviving, numFakeTracksHitPt2);
  printTrackStats("Fake tracks pT > 3: ", numFakeTracksPt3Surviving, numFakeTracksPt3);
  printTrackStats("Fake tracks pT > 5: ", numFakeTracksPt5Surviving, numFakeTracksPt5);

  // True tracks
  printTrackStats("True tracks: ", numTrueTracksPt0Surviving, numTrueTracksPt0);
  printTrackStats("True tracks pT > 2: ", numTrueTracksPt2Surviving, numTrueTracksPt2);
  printTrackStats("True tracks hit pT > 2: ", numTrueTracksHitPt2Surviving, numTrueTracksHitPt2);
  printTrackStats("True tracks pT > 3: ", numTrueTracksPt3Surviving, numTrueTracksPt3);
  printTrackStats("True tracks pT > 5: ", numTrueTracksPt5Surviving, numTrueTracksPt5);

  std::cout << "Average fake rate (pT > 2): " << 100*float(numFakeTracksPt2Surviving)/float(numRecoTracksPt2Surviving) << " %" << std::endl;
  std::cout << "Average fake rate (pT > 3): " << 100*float(numFakeTracksPt3Surviving)/float(numRecoTracksPt3Surviving) << " %" << std::endl;
  std::cout << "Average fake rate (pT > 5): " << 100*float(numFakeTracksPt5Surviving)/float(numRecoTracksPt5Surviving) << " %" << std::endl;

  std::cout << "Unique hits momenta: "   << uniqueHitMomenta.size() << std::endl;

  // Fill histograms (these probably wont work, since they should be event quantities)  
  plots->nDelphesHits.at(branchCounter)->Fill(numRecoTracksPt0);
  plots->nDelphesHitsPt2.at(branchCounter)->Fill(numRecoTracksPt2);
  plots->nRecoTracksMatched.at(branchCounter)->Fill(numTrueTracksPt0);
  plots->nRecoTracksMatchedPt2.at(branchCounter)->Fill(numTrueTracksPt2);
  plots->nRecoTracksPt2.at(branchCounter)->Fill(numRecoTracksPt2);


  // Store results in a json
  json output;
  output["FakeOriginal"]  = { {"Pt0", numFakeTracksPt0},          {"Pt2", numFakeTracksPt2},          {"Pt3", numFakeTracksPt3},          {"Pt5", numFakeTracksPt5} };
  output["FakeSurviving"] = { {"Pt0", numFakeTracksPt0Surviving}, {"Pt2", numFakeTracksPt2Surviving}, {"Pt3", numFakeTracksPt3Surviving}, {"Pt5", numFakeTracksPt5Surviving} }; 

  // 
  output["TrueOriginal"]  = { {"Pt0", numTrueTracksPt0},          {"Pt2", numTrueTracksPt2},          {"Pt3", numTrueTracksPt3},          {"Pt5", numTrueTracksPt5} };
  output["TrueSurviving"] = { {"Pt0", numTrueTracksPt0Surviving}, {"Pt2", numTrueTracksPt2Surviving}, {"Pt3", numTrueTracksPt3Surviving}, {"Pt5", numTrueTracksPt5Surviving} }; 

  // 
  output["TrueOriginalHit"]  = { {"Pt2", numTrueTracksHitPt2} };
  output["TrueSurvivingHit"] = { {"Pt2", numTrueTracksHitPt2Surviving} };

  output["FakeOriginalHit"]  = { {"Pt2", numFakeTracksHitPt2} };
  output["FakeSurvivingHit"] = { {"Pt2", numFakeTracksHitPt2Surviving} };

  std::cout << "Finished Loop()\n" << std::endl;
  return output; 


}

# ifndef __CINT__
int main(int argc, char* argv[]){

  TString sampleName = argv[1]; 
  TString outputFile = argv[2];

  int cutsOption     = atoi(argv[3]);
  int nEventsMax     = atoi(argv[4]);

  float BDTCut        = atof(argv[5]);
  float deltaKappaCut = atof(argv[6]);

  // print input arguments
  std::cout << "Input arguments:" << std::endl;

  std::cout << "sampleName: " << sampleName << std::endl;
  std::cout << "outputFile: " << outputFile << std::endl;
  std::cout << "cutsOption: " << cutsOption << std::endl;
  std::cout << "nEventsMax: " << nEventsMax << std::endl;
  std::cout << "BDTCut: " << BDTCut << std::endl;
  std::cout << "deltaKappaCut: " << deltaKappaCut << std::endl;

  if(!(cutsOption == 0 || cutsOption == 1 || cutsOption == 2)){
    std::cerr << "ERROR: arg3 must be 0, 1 or 2" << std::endl;
    std::cout << "You entered: " << cutsOption << std::endl;
    return 0;
  }

  if(cutsOption == 0){
    std::cout << "User selected to use rectangular cuts" << std::endl;
    std::cout << "deltaKappa cut: " << deltaKappaCut << std::endl;
  }
  if(cutsOption == 1){ 
    std::cout << "User selected to use BDT" << std::endl;
    std::cout << "BDT cut: " << BDTCut << std::endl;
  }
  if(cutsOption == 2){
    std::cout << "You selected both BDT and rectangular cuts" << std::endl;
    std::cout << "deltaKappa cut: " << deltaKappaCut << std::endl;
    std::cout << "BDT cut: " << BDTCut << std::endl;
  }

  if(nEventsMax == -1){
    std::cout << "Will run over all events" << std::endl;
  }
  else{
    std::cout << "Will run over: " << nEventsMax << " events." << std::endl;
  }

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
  zresiduumOverEtaCut["Tracks10"] = 0.1; // 10 mm
  zresiduumOverEtaCut["Tracks20"] = 0.4; // 20 mm
  zresiduumOverEtaCut["Tracks30"] = 0.4; // 30 mm
  zresiduumOverEtaCut["Tracks40"] = 0.4; // 40 mm
  zresiduumOverEtaCut["Tracks50"] = 0.4; // 50 mm


  json jStore;

  int branchCounter(0);
  for(auto branchName : branchNames){
    trackLoop theLoop(branchName, sampleName);
    json output = theLoop.Loop(plots, branchCounter, tripletLayouts[branchName], zresduumCuts[branchName], cutsOption, nEventsMax, BDTCut, deltaKappaCut);

    jStore[static_cast<std::string>(branchName)] = output;
    branchCounter++;
  }


  // write  histograms 
  std::cout << "Writing to file: " << outputFile << std::endl;
  result->Write(outputFile);

  // Write a nice json file to store stuff :)
  TString jFileName = outputFile.ReplaceAll("root", "json");
  std::cout << "Writing to json file: " << jFileName << std::endl;
  std::ofstream ojsonfile(jFileName);
  ojsonfile << std::setw(4) << jStore << std::endl;

  // cleanup
  delete plots;
  delete result;

  return 0;
}
# endif
