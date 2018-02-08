#ifndef TrackFitter_h
#define TrackFitter_h

#include "classes/DelphesClasses.h"
#include "classes/HitCollection.h"
#include "classes/lineParameters.h"
#include "classes/myTrack.h"


#include <vector>
#include <map>
#include <iostream>
#include <string>

using hitContainer = std::map<int, std::vector<Hit*> >; 

// enum for different fit types (can add more later) 
enum fitTypes{
  linearOutToIn,
  linearInToOut,
  linearInnerAndOuterMatching, // Anna's suggestion
  simpleLinear,
  MAX
};



// class to 
class TrackFitter{

  private:

    bool m_debug;

    // private member variables
    std::vector< myTrack > m_tracks;
    fitTypes fitType;
    std::vector<float> m_parameters;
    std::vector<HitCollection> m_associatedHitCollection;
    std::vector<int> m_layerIDs;

    // algorithms to associate hits 
    bool associateHitsLinearOutToIn(hitContainer, float, float);
    bool associateHitsLinearInToOut(hitContainer, float, float);
    bool associateHitsSimple(hitContainer, float, float);

    // functions to return tracks from hit collections
    bool combineHitsToTracksInToOut(); 
    bool combineHitsToTracksOutToIn(); 
    bool combineHitsToTracksMatchingInnerAndOutermost();

    // functions to calculate search windows
    float calculateZWindowForNextLevel(float, float, float, float);
    bool calculateRPhiWindowOutToIn(const float, const float, const float);
    float calculateRPhiWindowInToOut(const float, const float, const float);

    // track fitting functions
    /*myTrack simpleLinearLeastSquaresFit(std::vector< Hit* > hits) const;*/
    /*lineParameters simpleLinearLeastSquaresFit(std::vector< std::pair< float, float> >) const; */

  public:
    
    // constructor 
    TrackFitter(const fitTypes ftIn, std::vector<float> paramIn, std::vector<int> layerIDs){
      fitType=ftIn;
      m_parameters=paramIn;
      m_debug = false;
      m_layerIDs = layerIDs;
    };

    TrackFitter(){
      fitType=MAX;
      std::vector<float> emptyVec;
      m_parameters = emptyVec;
      m_debug = false;
    }

    void debug(){ m_debug = true; }

    // public functions
    bool AssociateHits(hitContainer hc);
    std::vector <myTrack> GetTracks();

};
#endif // TrackFitter_h
