#ifndef TrackReconstructor_h
#define TrackReconstructor_h

/** \class TrackReconstructor
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


class TrackReconstructor: public DelphesModule
{
public:

  TrackReconstructor();
  ~TrackReconstructor();

  void Init();
  void Process();
  void Finish();
  std::vector<TLorentzVector> ParticlePropagator(float, float, bool, bool);

private:

  Int_t fNVertexToAssociate;

  const TObjArray *fInputArray;
  const TObjArray *fBeamSpotInputArray; //!
  TIterator *fItInputArray;

  Double_t fBz; // magnetic field in z direction [T]
  Double_t fBarrelLength; // length of the tracker barrel [m] (only in the positive z direction, so it's the "half length") 
  Double_t fTrackPtMin; // minimum track pT to be considered 
  Double_t fEndCapRadius; // raduis of endcap discs [m] 

  std::vector<Barrel> fBarrelLayers;
  std::vector<float> fEndcapZPositions; 


  TObjArray *fHitOutputArray;

  ClassDef(TrackReconstructor, 1)
};

#endif
