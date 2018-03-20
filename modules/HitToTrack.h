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

#include <string>
#include <vector>
#include "TLorentzVector.h"

class TObjArray;
class TIterator;


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

  const TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fTrackOutputArray;

  ClassDef(HitToTrack, 1)
};

#endif
