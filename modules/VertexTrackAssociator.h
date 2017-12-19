#ifndef VertexTrackAssociator_h
#define VertexTrackAssociator_h

/** \class VertexTrackAssociator
 *
 *  Associate tracks to a vertex  
 *
 *  \authors W. Fawcett
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

  Int_t fNVertexToAssociate;

  TObjArray *fTrackInputArray;
  TObjArray *fVertexInputArray;
  TIterator *fItTrackInputArray;
  TIterator *fItVertexInputArray;

  TObjArray *fOutputArray;

  ClassDef(VertexTrackAssociator, 1)
};

#endif
