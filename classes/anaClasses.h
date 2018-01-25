#ifndef anaClasses_h
#define anaClasses_h

#include "classes/DelphesClasses.h"

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>

using hitContainer = std::map<int, std::vector<Hit*> >; 

// enum for different fit types (can add more later) 
enum fitTypes{
  linearOutToIn,
  linearInToOut,
  MAX
};


// object to store track parameters (TO DO)! 
class myTrack{
  public:
    float coord;
};

// Class to store essentially a vector of hits, but the hits know about their relationship to one another
class HitCollection{
  private:
    Hit * m_hit; 
    TLorentzVector m_position;

    std::vector<HitCollection*> m_assignedHits; 

  public:
    int SurfaceID;

    // constructor
    HitCollection(Hit * hitIn){
      m_hit = hitIn; 
      m_position = hitIn->Position();
      SurfaceID = hitIn->SurfaceID;
    };

    // default constructor
    HitCollection(){
      m_hit = 0;
      m_position = TLorentzVector();
      SurfaceID = -1; 
    };

    // other hits matched to this hit 
    void addHit(HitCollection* collectionIn){ m_assignedHits.push_back( collectionIn ); }

    int countAssignedHits(){ return m_assignedHits.size(); }

    // print functions 
    void printAssignedHits();
    void printHit();
    std::string hitInfo();

    // access some of the TLorentzVector functions
    float X() {return m_hit->X;}
    float Y() {return m_hit->Y;}
    float Z() {return m_hit->Z;}
    float T() {return m_hit->T;}
    float Perp() {return m_position.Perp();}
};

// class to 
class TrackFitter{

  private:

    // private member variables
    std::vector< myTrack > tracks; 
    fitTypes fitType;
    std::vector<float> m_parameters;
    /*std::vector< std::vector< Hit* >> associatedHitCollection; */
    std::vector<HitCollection> m_associatedHitCollection;

    // algorithms to associate hits 
    bool associateHitsLinearOutToIn(hitContainer, float, float);
    bool associateHitsLinearInToOut(hitContainer, float, float);

    bool combineHitsToTracksInToOut(); 
    bool combineHitsToTracksOutToIn(); 

    float calculateZWindowForNextLevel(float, float, float, float);

  public:
    
    // constructor 
    TrackFitter(const fitTypes ftIn, std::vector<float> paramIn){
      fitType=ftIn;
      m_parameters=paramIn;
    };

    TrackFitter(){
      fitType=MAX;
      std::vector<float> emptyVec;
      m_parameters = emptyVec;
    }

    // public functions
    bool AssociateHits(hitContainer hc);
    std::vector <myTrack> GetTracks();

};




#endif // anaClasses_h
