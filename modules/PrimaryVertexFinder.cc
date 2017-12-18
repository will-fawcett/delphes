/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class PrimaryVertexFinder
 *
 *  Selects vertices from the (track) InputArray
 *  Returns 1 vertex per event, and the tracks that have been associated to that vertex 
 *
 *  \author W. Fawcett 
 *
 */

#include "modules/PrimaryVertexFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

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

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

PrimaryVertexFinder::PrimaryVertexFinder() :
  fItTrackInputArray(0)
{
}

//------------------------------------------------------------------------------

PrimaryVertexFinder::~PrimaryVertexFinder()
{
}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Init()
{

  // read parameters
  // 1: primary bin (sliding window)
  // 2: could be something else (to be implemented)  
  fSearchAlgorithm = GetInt("SearchAlgorithm", 1);

  // Parameters for sliding window algorithm 
  fBeamMaxWidth = GetInt("BeamWidth", 600); // [mm] maximum length of the beam to be considered 
  fBinWidth = GetInt("BinWidth", 1.0); // [mm] width of bins to be considered in sliding window algorithm 
  fBinStepSize = GetDouble("BinStepSize", 0.1); // [mm] step size of bins for sliding window 

  // Calculate properties for sliding window based on input properties 
  fNBins = fBeamMaxWidth/fBinWidth; // number of histograms bins
  fNHistograms = fBinWidth/fBinStepSize; // number of histograms for sliding window alg
  fBeamMaxZ = fBeamMaxWidth/2;
  fBeamMinZ = -1*fBeamMaxZ; 

  // import input array(s)
  fTrackInputArray   = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output array(s)
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
  fTrackOutputArray  = ExportArray(GetString("TrackOutputArray", "track"));

}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
}

//------------------------------------------------------------------------------

void PrimaryVertexFinder::Process()
{
  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  // Define histograms for sliding window alg. 
  std::vector<TH1F*> windowHists;
  for(int i=0; i< fNHistograms; ++i){
    windowHists.push_back( new TH1F( (std::to_string(i)+"window").c_str(), "", fNBins, fBeamMinZ+fBinStepSize*i, fBeamMaxZ+fBinStepSize*i) );
  }
  

  // loop over all input candidates
  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;


    /*****************

    // define histograms for sliding window algorithm
    //std::vector<TH1F*> windowHists;
    //std::vector<TH1F*> nTrackHists;
    std::vector< std::pair< TH1F*, TH1F*> > hists;
    for(int i=0; i< binWidth/slideStep; ++i){

      //windowHists.push_back( new TH1F( (std::to_string(i)+"window").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i) );
      //nTrackHists.push_back( new TH1F( (std::to_string(i)+"nTrack").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i) );
      TH1F * nTrackHist = new TH1F( (std::to_string(i)+"nTrack").c_str(), "", nBins, beamMinZ+slideStep*i, beamMaxZ+slideStep*i );
      hists.push_back( std::make_pair(windowHist, nTrackHist) );
    }

    // loop over all tracks, fill histograms
    for(auto itTrack = branchTrack->begin(); itTrack != branchTrack->end(); ++itTrack){
      float trackPt   = dynamic_cast<Track*>(*itTrack)->PT;
      float zPosition = dynamic_cast<Track*>(*itTrack)->DZ;
      //for(auto hist : windowHists)  hist->Fill(zPosition, trackPt);
      for(auto histPair : hists){
        histPair.first->Fill(zPosition, trackPt);
        histPair.second->Fill(zPosition, 1); // number of tracks, maybe could have used a TH2D?
      }
    }

    // Extract the largest sum(pT)
    float previousMaxPt(0), zBinLow(0), zBinHigh(0), zBinWidth(0);
    int nTracks(0);
    //for(auto hist : windowHists){
    for(auto histPair : hists){
      Int_t binWithPtMax = histPair.first->GetMaximumBin();
      TAxis * xaxis = static_cast<TAxis*>(histPair.first->GetXaxis());
      float maxPt = histPair.first->GetBinContent(binWithPtMax);
      if(maxPt > previousMaxPt){
        previousMaxPt = maxPt;
        zBinLow   = xaxis->GetBinLowEdge(binWithPtMax);
        zBinHigh  = xaxis->GetBinUpEdge(binWithPtMax);
        zBinWidth = xaxis->GetBinWidth(binWithPtMax);
        nTracks = histPair.second->GetBinContent(binWithPtMax);
      }
    }

    // store results
    std::map<std::string, float> results;
    results["zMin"]   = zBinLow;
    results["zMax"]   = zBinHigh;
    results["zWidth"] = zBinWidth;
    results["nTracks"] = nTracks;

    // delete histograms
    //for(auto hist : windowHists){
      //delete hist;
    // }
    windowHists.clear();
    for(auto histPair : hists){
      delete histPair.first;
      delete histPair.second;
    }
    hists.clear();

    return results;
    ********/


    fVertexOutputArray->Add(candidate);
    fTrackOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
