#ifndef VertexTrackAssociator_h
#define VertexTrackAssociator_h

/** \class VertexTrackAssociator
 *
 *  Cluster vertices from tracks
 *
 *  \authors A. Hart, M. Selvaggi
 *
 */


#include "classes/DelphesModule.h"

#include <string>
#include <vector>
#include <map>

class TObjArray;
class TIterator;

class VertexTrackAssociator: public DelphesModule
{
public:

  VertexTrackAssociator();
  ~VertexTrackAssociator();

  void Init();
  void Process();
  void Finish();

private:

  void createSeeds ();
  void growCluster (const UInt_t);
  Double_t weight (const UInt_t);
  void addTrackToCluster (const UInt_t, const UInt_t);
  void removeTrackFromCluster (const UInt_t, const UInt_t);

  Double_t fSigma;
  Double_t fMinPT;
  Double_t fMaxEta;
  Double_t fSeedMinPT;
  Int_t fMinNDF;
  Int_t fGrowSeeds;

  TObjArray *fTrackInputArray;
  TIterator *fItTrackInputArray;
  TObjArray *fVertexInputArray;
  TIterator *fItVertexInputArray;

  TObjArray *fOutputArray;
  TObjArray *fVertexOutputArray;

  std::map<UInt_t, std::map<std::string, Double_t> > trackIDToDouble;
  std::map<UInt_t, std::map<std::string, Int_t> > trackIDToInt;
  std::map<UInt_t, std::map<std::string, Bool_t> > trackIDToBool;

  std::map<UInt_t, std::map<std::string, Double_t> > clusterIDToDouble;
  std::map<UInt_t, std::map<std::string, Int_t> > clusterIDToInt;
  std::map<UInt_t, std::map<std::string, Bool_t> > clusterIDToBool;
  std::vector<std::pair<UInt_t, Double_t> > trackPT;
  std::vector<std::pair<UInt_t, Double_t> > clusterSumPT2;

  ClassDef(VertexTrackAssociator, 1)
};

#endif
