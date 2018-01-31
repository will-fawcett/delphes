#ifndef HitFinder_h
#define HitFinder_h

/** \class HitFinder
 *
 *  Associate tracks to a vertex  
 *
 *  \authors W. Fawcett
 *
 */


#include "classes/DelphesModule.h"
#include "classes/Barrel.h"

#include <string>
#include <vector>
#include <map>
#include "TLorentzVector.h"

class TObjArray;
class TIterator;
class Barrel;


class HitFinder: public DelphesModule
{
public:

  HitFinder();
  ~HitFinder();

  void Init();
  void Process();
  void Finish();
  void ParticlePropagator(float, float, int, bool, bool);

private:

  Int_t fNVertexToAssociate;

  const TObjArray *fInputArray;
  const TObjArray *fBeamSpotInputArray; //!
  TIterator *fItInputArray;

  Double_t fBz; // magnetic field in z direction [T]
  Double_t fBarrelLength; // length of the tracker barrel [m] (only in the positive z direction, so it's the "half length") 
  Double_t fTrackPtMin; // minimum track pT to be considered 
  Double_t fEndCapRadius; // raduis of endcap discs [m] 

  std::vector<float> fBarrelLayerRadii;
  std::vector<float> fEndcapZPositions; 


  TObjArray *fHitOutputArray;

  ClassDef(HitFinder, 1)
};

#endif
