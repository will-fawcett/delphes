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

// global
bool m_debug = false;

// Fast way to test if a file exists
inline bool fileExists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

//------------------------------------------------------------------------------

std::string PrintTLorentz(TLorentzVector &v){
  std::ostringstream s;
  s <<  " (" << v.Pt() << ", " << v.Eta() << ", " << v.Phi() << ", " << v.M() << ")";
  return s.str();
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


// Rule for sorting vector of tracks 
bool reverse(const Track* i, const Track* j){
  return i->PT > j->PT; 
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

//------------------------------------------------------------------------------


struct Plots
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

void BookHistograms(ExRootResult *result, Plots *plots)
{

  // Track occupancy
  plots->trackOccupancy = result->AddHist2D("trackOccupancy", "", "Track multi per tower", "Track pT [GeV]", 100, 0, 100, 100, 0, 300);
  plots->trackOccupancyProf = result->AddProfile("trackOccupancyProf", "", "Track multi per tower", "Number of towers", 100, 0, 100);




  plots->misidPVLogx = result->AddHist1D("misidPVLogx", "Misidentified primary vertices", "Distance between PB and true PV [mm]", "Number of events", 200, 0, 400, 1, 0); // probably will be some overflow
  plots->misidPV = result->AddHist1D("misidPV", "Misidentified primary vertices", "Distance between PB and true PV [mm]", "Number of events", 200, 0, 200, 0, 0); // probably will be some overflow

  // primary bin 
  //plots->trackMultiplicityInPB    = result->AddHist1D("trackMultiplicityInPB", "Track multiplicity in PB", "", "", 100, 0, 100);
  //plots->vertexMultiplicityInPB   = result->AddHist1D("vertexMultiplicityInPB", "Track multiplicity in PB", "", "", 100, 0, 100);
  //plots->nVertexVsZVertexPosition = result->AddHist2D("nVertexVsZVertexPosition", "Number of vertices in PB", "z position of true PV", 100, 0, 100, 100, -300, 300);
  //plots->nVertexVsZPBPosition     = result->AddHist2D("nVertexVsZPBPosition", "Number of vertices in PB", "z position of PB", 100, 0, 100, 100, -300, 300);
  //plots->PBvZVertexPosition       = result->AddHist2D("PBvZVertexPosition", "z position of PB", "z position of PV", 100, -300, 300, 100, -300, 300);



  plots->nJets           = result->AddHist1D( "nJets", "nJets", "Number of Jets", "", 100, 0, 100, 0, 0 );
  plots->nAssociatedJets = result->AddHist1D( "nAssociatedJets", "nAssociatedJets", "Number of Jets", "", 100, 0, 100, 0, 0 );
  plots->allJetPt        = result->AddHist1D("allJetPt", "", "Jet p_{T} (all jets) [GeV]", "", 100, 0, 1000);

  // z pt
  plots->binnedZpT = result->AddHist1D( "binnedZpT", "", "z position [mm]", "Sum(p_{T}) [GeV]", 600, -300, 300);



}


//------------------------------------------------------------------------------

void AnalyseEvents(const int nEvents, ExRootTreeReader *treeReader, Plots *plots)
{


  // Define branches
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  //TClonesArray *branchTruthTrack = treeReader->UseBranch("TruthTrack");
  //TClonesArray *branchTrack      = treeReader->UseBranch("Track");
  //TClonesArray *branchTrackJet   = treeReader->UseBranch("TrackJet");

  // Normal vertex
  TClonesArray *branchVertex     = treeReader->UseBranch("Vertex");
  // primary bin from normal delphes tracks
  TClonesArray *branchPrimaryBin = treeReader->UseBranch("PrimaryBin");

  TClonesArray *branchPrimaryBinFromHits = treeReader->UseBranch("PrimaryBinHits");

  //TClonesArray *branchTower      = treeReader->UseBranch("Tower");
  //TClonesArray *branchPileupParticle = treeReader->UseBranch("PileupParticle");



  // Loop over all events
  Long64_t allEntries = treeReader->GetEntries();
  int nEventsCorrectlyIdentifiedVertex(0);
  int vertexIsInPB(0);
  int vertexIsInPBHits(0);
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


    // scan parameters
    int beamMinZ    = -300; // mm
    int beamMaxZ    = 300; // mm
    float binWidth  = 1.0; // mm
    float slideStep = 0.1; // mm



    //std::map<std::string, float> PBInfo = findPrimaryBin(branchTrack, binWidth, slideStep, beamMinZ, beamMaxZ);
    //float zMin   = PBInfo["zMin"];
    //float zMax   = PBInfo["zMax"];
    //int nTracksInPB = PBInfo["nTracks"];
    //float zWidth = PBInfo["zWidth"];
    //float PBCentroid = PBInfo["zPosition"]; 

  
    ////////////////////////////////////////////////////////////////
    // Find location of the primary bin (from normal Delphes tracks)
    ////////////////////////////////////////////////////////////////
    int numPrimaryBins = branchPrimaryBin->GetEntriesFast();
    float PB_Z(0);
    float PB_Z_max(0), PB_Z_min(0);
    //std::cout << "There are " << numPrimaryBins << " nominal primary bins" << std::endl;
    for(int i=0; i<numPrimaryBins; ++i){
      Vertex * primaryBin = (Vertex*) branchPrimaryBin->At(i);
      PB_Z = primaryBin->Z; 
    }
    if(numPrimaryBins > 1){
      std::cout << "ERROR: more than 1 primary bin" << std::endl;
    }
    // 1mm bin width
    PB_Z_max = PB_Z + 1.0;
    PB_Z_min = PB_Z - 1.0;

    ////////////////////////////////////////////////////////////////
    // Find location of the primary bin (from tracks reconstructed from hits)
    ////////////////////////////////////////////////////////////////
    int numPrimaryBinsFromHits = branchPrimaryBinFromHits->GetEntriesFast();
    float hitsPB_Z(0);
    float hitsPB_Z_max(0), hitsPB_Z_min(0);
    //std::cout << "There are " << numPrimaryBinsFromHits << " nominal primary bins" << std::endl;
    for(int i=0; i<numPrimaryBinsFromHits; ++i){
      Vertex * primaryBin = (Vertex*) branchPrimaryBinFromHits->At(i);
      hitsPB_Z = primaryBin->Z; 
    }
    if(numPrimaryBinsFromHits > 1){
      std::cout << "ERROR: more than 1 primary bin" << std::endl;
    }
    // 1mm bin width
    hitsPB_Z_max = hitsPB_Z + 1.0;
    hitsPB_Z_min = hitsPB_Z - 1.0;


    ///////////////////////////////////
    // Find the location of the real PV
    // count how accurate locating this is 
    ///////////////////////////////////
    float vertexZ(0.0);
    int nVerticesInPB(0);
    int nVerticesInHitsPB(0);
    int nPrimaryVertices(0);
    //std::cout << "Number of primary vertices: " << branchVertex->GetEntriesFast() << std::endl;
    for(int i=0; i<branchVertex->GetEntriesFast(); ++i){
      Vertex * vertex = (Vertex*) branchVertex->At(i);
      if(vertex->Z < PB_Z_max && vertex->Z > PB_Z_min) nVerticesInPB++;
      if(vertex->Z < hitsPB_Z_max && vertex->Z > hitsPB_Z_min) nVerticesInHitsPB++;
      if(vertex->IsPU) continue;
      else{
        vertexZ = vertex->Z;
        nPrimaryVertices++;
        //std::cout << "Event: " << entry << " vertex z: " << vertexZ << std::endl;
      }
      // is _the_ primary vertex inside the primary bin?
      if(vertex->Z < PB_Z_max && vertex->Z > PB_Z_min) vertexIsInPB++;
      if(vertex->Z < hitsPB_Z_max && vertex->Z > hitsPB_Z_min) vertexIsInPBHits++;
      
    }
    if(nPrimaryVertices > 1){
      std::cout << "ERROR: more than one true PV: " << nPrimaryVertices << std::endl; 
    }

    //trackMultiplicityInPB->Fill(nTracksInPB);
    //vertexMultiplicityInPB->Fill(nVerticesInPB);
    //nVertexVsZVertexPosition->(nVerticesInPB, vertexZ);
    //nVertexVsZPBPosition->(nVertexVsZPBPosition, PBCentroid);
    //PBvZVertexPosition->(PBCentroid, vertexZ)
    


    /********************************
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
    // ******* (start block comment)
      std::cout << "Vertex isn't inside PB" << std::endl;
      std::cout << "vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]\t" << fabs(vertexZ - (zMax+zMin)/2) << "\t nVertices: " << nVerticesInPB << std::endl;
      std::cout << "nTracks in PB " << nTracksInPB << std::endl;
      plots->misidPV->Fill( fabs(vertexZ - (zMax+zMin)/2) );
      plots->misidPVLogx->Fill( fabs(vertexZ - (zMax+zMin)/2) );
      //std::cout << "vertex z: " << vertexZ << " PB range: [" << zMin << ", " << zMax << "]" << std::endl;
      // *********** (end block comment)
    }
  **************************/







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

  std::cout << "Of " << allEntries << " events: " << std::endl;
  std::cout << "\t" <<  vertexIsInPB << " had the vertex correctly identified (nominal), i.e. " << static_cast<float>(vertexIsInPB)/static_cast<float>(allEntries) << " of events."  << std::endl;
  std::cout << "\t" << vertexIsInPBHits << " had the vertex correctly identified (hits), i.e. " << static_cast<float>(vertexIsInPBHits)/static_cast<float>(allEntries) << " of events."  << std::endl;
} // end AnalyseEvents



//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, Plots *plots)
{
  result->Print("pdf");
}


//------------------------------------------------------------------------------




//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{

  gROOT->SetBatch(1);
  //gROOT->ProcessLine("#include <vector>");

  std::string appName = "primaryVertexEfficiency";
  std::cout << "Will execute " << appName << std::endl;

  std::string outputFile;
  std::string inputGlob; bool inputGlobSet(false);
  std::string inputFile; bool inputFileSet(false);

  // declare input arguments
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
  BookHistograms(result, plots);
  AnalyseEvents(nEvents, treeReader, plots);

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
