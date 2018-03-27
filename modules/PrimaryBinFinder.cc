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


/** \class PrimaryBinFinder
 *
 *  Selects vertices from the (track) InputArray
 *  Returns 1 vertex per event, and the tracks that have been associated to that vertex 
 *
 *  \author W. Fawcett 
 *
 */

#include "modules/PrimaryBinFinder.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TH1.h"
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

PrimaryBinFinder::PrimaryBinFinder() :
  fItTrackInputArray(0)
{
}

//------------------------------------------------------------------------------

PrimaryBinFinder::~PrimaryBinFinder()
{
}

//------------------------------------------------------------------------------

void PrimaryBinFinder::Init()
{

  m_debug = GetBool("debug", false);

  // read parameters
  // 1: primary bin (sliding window)
  // 2: could be something else (to be implemented)  
  fSearchAlgorithm = GetInt("SearchAlgorithm", 1);

  // Parameters for sliding window algorithm 
  fBeamMaxWidth = GetInt("BeamWidth", 600); // [mm] maximum length of the beam to be considered 
  fBinWidth = GetInt("BinWidth", 1.0); // [mm] width of bins to be considered in sliding window algorithm 
  fBinStepSize = GetDouble("BinStepSize", 0.1); // [mm] step size of bins for sliding window 
  fNumPrimaryBins = GetInt("NumPrimaryBins", 1); // number of primary bins to return

  // Calculate properties for sliding window based on input properties 
  fNBins = fBeamMaxWidth/fBinWidth; // number of histograms bins
  fNHistograms = fBinWidth/fBinStepSize; // number of histograms for sliding window alg
  fBeamMaxZ = fBeamMaxWidth/2;
  fBeamMinZ = -1*fBeamMaxZ; 


  // import input array(s)
  fTrackInputArray   = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output array(s)
  fVertexOutputArray = ExportArray(GetString("OutputArray", "vertices"));
  //fTrackOutputArray  = ExportArray(GetString("TrackOutputArray", "track"));

}

//------------------------------------------------------------------------------

void PrimaryBinFinder::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
}

//------------------------------------------------------------------------------

// some functions/structs for use in Process 

struct bin
{
  Double_t binLowEdge;
  Double_t binUpEdge;
  Double_t binCenter;
  Double_t binWidth;
};

std::string formatBin(const bin& ibin){
  // useful for debugging 
  std::ostringstream s;
  s << "bin center: " << ibin.binCenter << " bin range: [" << ibin.binLowEdge << ", " << ibin.binUpEdge << "]";
  return s.str();
}


bool reverseSort(const std::pair<float, bin> i, const std::pair<float, bin> j){
  return i.first > j.first; 
}

void printVec(std::vector<std::pair<float, bin> >& vec){
  // useful for debugging 
  for(auto p : vec){
    std::cout << "pt: " << p.first << " " << formatBin(p.second) <<  std::endl;
  }
}

//------------------------------------------------------------------------------

void PrimaryBinFinder::Process()
{
  Candidate *candidate;

  // Define histograms for sliding window alg. 
  std::vector<TH1F*> windowHists;
  for(int i=0; i< fNHistograms; ++i){
    windowHists.push_back( new TH1F( (std::to_string(i)+"window").c_str(), "", fNBins, fBeamMinZ+fBinStepSize*i, fBeamMaxZ+fBinStepSize*i) );
  }

  // loop over all input candidates, fill histograms
  fItTrackInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next())))
  {
    for(int i=0; i<fNHistograms; ++i) windowHists.at(i)->Fill(candidate->DZ, candidate->PT);
  }

  // loop over all candidates (tracks), extract sum(pT) information from each bin 
  std::vector< std::pair< float, bin > > binArray; 
  binArray.clear();
  if(fNumPrimaryBins > 1){
    // What if only one vertex in an isolated region, many bins will have the same pT ... ? 
    // Extract all the bins from all the histograms, sort into pT order. 
    for(int i=0; i<fNHistograms; ++i){
      for(int ibin=1; ibin<windowHists.at(i)->GetNbinsX()+1; ++ibin){ // start counting at 1, bin 0 is underflow bin. GetNbinsX() is the last bin, hence the +1. 
        
        TAxis * xaxis = static_cast<TAxis*>(windowHists.at(i)->GetXaxis());
        float binPt = windowHists.at(i)->GetBinContent(ibin);
        bin thisBin = {
          xaxis->GetBinLowEdge(ibin),
          xaxis->GetBinUpEdge(ibin),
          xaxis->GetBinCenter(ibin),
          xaxis->GetBinWidth(ibin)
        };

        binArray.push_back( std::make_pair(binPt, thisBin) );

      }
    }

    // sort bin array 
    std::sort(binArray.begin(), binArray.end(), reverseSort); 

  }
  else{

    // if there is just 1 bin requested (faster algorithm) 
    // Extract the largest sum(pT)
    float previousMaxPt(0), zBinLow(0), zBinHigh(0), zBinWidth(0), zBinCenter(0);
    for(int i=0; i<fNHistograms; ++i){
      Int_t binWithPtMax = windowHists.at(i)->GetMaximumBin();
      TAxis * xaxis = static_cast<TAxis*>(windowHists.at(i)->GetXaxis());
      float maxPt = windowHists.at(i)->GetBinContent(binWithPtMax);
      if(maxPt > previousMaxPt){
        previousMaxPt = maxPt;
        zBinLow   = xaxis->GetBinLowEdge(binWithPtMax);
        zBinHigh  = xaxis->GetBinUpEdge(binWithPtMax);
        zBinWidth = xaxis->GetBinWidth(binWithPtMax);
        zBinCenter = xaxis->GetBinCenter(binWithPtMax);
      }
    }

    bin maxPtBin = {
      zBinLow, 
      zBinHigh, 
      zBinCenter,
      zBinWidth
    };
    

    // only ever "pushed back" once, 
    binArray.push_back( std::make_pair( previousMaxPt, maxPtBin) );


  } // end of case with only 1 primary bin 

  // Create a new ``vertex'' that corresponds to the primary bin(s)
  DelphesFactory *factory = GetFactory();
  for(int ivertex=0; ivertex<fNumPrimaryBins; ++ivertex){
    Candidate *newVertex;
    newVertex = factory->NewCandidate();
    newVertex->DZ      = binArray.at(ivertex).second.binCenter;
    newVertex->ErrorDZ = binArray.at(ivertex).second.binWidth/2; // divide by 2, since error is +/-  
    newVertex->Position.SetXYZT(0.0, 0.0, newVertex->DZ, 0.0);
    newVertex->PositionError.SetXYZT(0.0, 0.0, newVertex->ErrorDZ, 0.0);
    fVertexOutputArray->Add(newVertex);
  }

  // Check that fVertexOutputArray has the same size as the number of primary bins requested
  if(fVertexOutputArray->GetEntriesFast() != fNumPrimaryBins){
    throw runtime_error("PrimaryBinFinder: Mismatch between number of identified primary bins and the number requested by the user");
  }

  // Delete histograms
  for(int i=0;i<fNHistograms;++i){
    delete windowHists.at(i);
  }
  windowHists.clear();

}

//------------------------------------------------------------------------------
