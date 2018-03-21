#ifndef FakeTrackRemover_h
#define FakeTrackRemover_h

/** \class FakeTrackRemover
 *
 *  Associate tracks to a vertex  
 *
 *  \authors W. Fawcett
 *
 */

#include "classes/DelphesModule.h"


class TObjArray;
class TIterator;

class FakeTrackRemover: public DelphesModule
{
public:

  FakeTrackRemover();
  ~FakeTrackRemover();

  void Init();
  void Process();
  void Finish();

private:

  // max DeltaPhi between layers i and j based on a 2 GeV track 
  float m_deltaPhi12Limit;
  float m_deltaPhi23Limit;
  float m_deltaPhi13Limit;

  bool m_debug;

  TObjArray *fTrackInputArray;
  TIterator *fItTrackInputArray;

  TObjArray *fOutputArray;

  ClassDef(FakeTrackRemover, 1)
};

#endif
