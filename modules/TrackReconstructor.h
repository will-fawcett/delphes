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
  std::vector<TLorentzVector> ParticlePropagator(float, float);

private:

  Int_t fNVertexToAssociate;

  const TObjArray *fInputArray;
  const TObjArray *fBeamSpotInputArray; //!
  TIterator *fItInputArray;

  Double_t fBz; // magnetic field in z direction 
  Double_t fBarrelLength; 
  Double_t fDiscHeight;
  Double_t fTrackPtMin;

  std::vector<Barrel> fBarrelLayers;


  TObjArray *fOutputArray;

  ClassDef(TrackReconstructor, 1)
};

#endif
