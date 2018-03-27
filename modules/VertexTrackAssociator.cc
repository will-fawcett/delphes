/** \class VertexTrackAssociator
 *
 *  Associate tracks to a vertex  
 *
 *  \authors W. Fawcett
 *
 */


#include "modules/VertexTrackAssociator.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include "TVector3.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;


//------------------------------------------------------------------------------

VertexTrackAssociator::VertexTrackAssociator() :
  fNVertexToAssociate(0)
{
}

//------------------------------------------------------------------------------

VertexTrackAssociator::~VertexTrackAssociator()
{
}

//------------------------------------------------------------------------------

void VertexTrackAssociator::Init()
{

  m_debug = GetBool("debug", false);

  // This module will only associate tracks to a given vertex
  // - there may be many vertices in the event
  // - the primary vertex will be the one with the largest sumPT, with index 0 (the default) 
  // - any other vertices can be selected with an index > 0, and will have a lower sum(pT) 
  fNVertexToAssociate = GetInt("IndexOfVertexToAssociate", 0);

  // Input arrays
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "VertexFinder/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fVertexInputArray = ImportArray(GetString("VertexInputArray", "VertexFinder/vertices"));
  fItVertexInputArray = fVertexInputArray->MakeIterator();

  // Output arrays 
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void VertexTrackAssociator::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItVertexInputArray) delete fItVertexInputArray;
}

//------------------------------------------------------------------------------

void VertexTrackAssociator::Process()
{
  ////////////////////////////////////////////////
  // loops over all input tracks
  // Tests if the track DZ is matches with the vertex z position (within uncertainties)
  // (performs geometrical matching between a track and a vertex)
  ////////////////////////////////////////////////
  
  // Extract vertex information
  if(fVertexInputArray->GetEntriesFast() < fNVertexToAssociate+1){
    throw runtime_error("VertexTrackAssociator: vertex index larger than number of identified vertices");
  }
  Candidate *vertex = static_cast<Candidate*>(fVertexInputArray->At(fNVertexToAssociate)); // assumes vertices already sorted by pT
  float vertexZ = vertex->Position.Z();
  float vertexZerror = vertex->PositionError.Z();
  float vertexZMax = vertexZ + vertexZerror;
  float vertexZMin = vertexZ - vertexZerror;
  

  // Match tracks to vertex 
  Candidate *track;
  fItTrackInputArray->Reset();
  int itrack(0);
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {

    itrack++;
    if(track->DZ < vertexZMax && track->DZ > vertexZMin)
    {
      fOutputArray->Add(track);
    }

  }


}

//------------------------------------------------------------------------------
