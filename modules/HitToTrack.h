#ifndef HitToTrack_h
#define HitToTrack_h

/** \class HitToTrack
 *
 * Takes hits produced by the HitFinder class
 * Performs track reconstruction
 * Returns Delphes tracks 
 *
 *  \authors W. Fawcett
 *
 */


#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include "classes/Location.h"

#include <string>
#include <vector>
#include <map>

class TObjArray;
class TIterator;

struct TrackParameterSet{

   Float_t pT;
   Float_t d0;
   Float_t z0;
   Float_t eta;
   Float_t theta;
   Float_t phi;
   Float_t kappa_013;
   Float_t kappa_123;
   Int_t isFake;
   Int_t Charge;
   Int_t isPU;

   //Float_t zresiduum;
   //Float_t beamlineIntersect;

};


class HitToTrack: public DelphesModule
{
public:

  HitToTrack();
  ~HitToTrack();

  void Init();
  void Process();
  void Finish();

private:

  bool m_debug;

  float fPhiSeedingWindowSize;
  float fEtaSeedingWindowSize;
  float fSeedingTrackPtMinGeV; 
  float fBField;

  int fLuminousRegionSize; 
  float fNVertxSigma; 
  float fZResiduumTolerance;

  std::vector< std::vector<Candidate*> > FindSeedsTriplet( 
    std::map<int, std::vector<Candidate*> >& hc,
    std::map<std::string, std::vector<Candidate*> >& hitMap, 
    Location& loc, 
    std::vector<int> layerIDs
    ) const;

  TrackParameterSet CalculateTrackParametersTriplet(std::vector<Candidate*>&) const; 
  float CalculateD0(float, float, float) const;
  bool isFake(std::vector<Candidate*>&) const;
  bool isTrackPU(std::vector<Candidate*>&) const;

  const TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fTrackOutputArray;

  ClassDef(HitToTrack, 1)
};

#endif
