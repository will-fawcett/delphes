//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  7 17:52:02 2018 by ROOT version 6.10/04
// from TTree Tracks10/A tree with tracks 10
// found on file: hits_ttbar_pu0_multiGeometry_tracks.root
//////////////////////////////////////////////////////////

#ifndef trackLoop_h
#define trackLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include "analysis/plotting.cpp"

// Header file for the classes stored in the TTree if any.

class trackLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         pT;
   Float_t         d0;
   Float_t         z0;
   Float_t         eta;
   Float_t         phi;
   Float_t         kappa_013;
   Float_t         kappa_123;
   Int_t           isFake;
   Float_t         zresiduum;
   Float_t         beamlineIntersect;
   Float_t         hit3pT;
   Float_t         hit3rho;
   Float_t         hit3eta;
   Float_t         hit3phi;
   Float_t         hit2pT;
   Float_t         hit2rho;
   Float_t         hit2eta;
   Float_t         hit2phi;
   Float_t         hit1pT;
   Float_t         hit1rho;
   Float_t         hit1eta;
   Float_t         hit1phi;

   // List of branches
   TBranch        *b_trackBranch;   //!

   trackLoop(TString branchName, TString sampleName, TTree *tree=0);
   virtual ~trackLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree* tree, TString branchIdentifier);
   virtual void     Loop(Plots* plots, int counter, std::vector<float> layerRadii, float zresiduumCut, int useBDT);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  private:
     TString m_branchName;
     TString m_sampleName; 
};

#endif

#ifdef trackLoop_cxx
trackLoop::trackLoop(TString branchName, TString sampleName, TTree *tree) : fChain(0), m_branchName(branchName), m_sampleName(sampleName)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.

  /*TString resultsPath = "/atlas/data4/userdata/wfawcett/delphes/results/hits_phiEtaSeg_tolerance05mm_phi2GeV_curvature0005_nVertexSigma5";*/
  /*TString filePath = resultsPath+"/hits_ttbar_pu0_multiGeometry_tracks.root";*/
  std::cout << "Opening file: " << m_sampleName << std::endl;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(m_sampleName);
    if (!f || !f->IsOpen()) {
      f = new TFile(m_sampleName);
    }
    std::cout << "With branch: " << m_branchName << std::endl;
    f->GetObject(m_branchName,tree);

  }
  m_branchName.ToLower(); // retarded ROOT function, modifies TString .
  Init(tree, m_branchName);
}

trackLoop::~trackLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trackLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trackLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void trackLoop::Init(TTree *tree, TString branchIdentifier)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // WJF: not sure why only one variable set in this way :/ 
   fChain->SetBranchAddress(branchIdentifier, &pT, &b_trackBranch); 
   Notify();
}

Bool_t trackLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trackLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trackLoop::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trackLoop_cxx
