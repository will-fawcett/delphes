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

    void printAssignedHitPointers(){
      for(HitCollection * h : m_assignedHits){
        std::cout << h << std::endl;
      }
    }

    // other hits matched to this hit 
    void addHit(HitCollection* collectionIn){ m_assignedHits.push_back( collectionIn ); }

    int countAssignedHits() const { return m_assignedHits.size(); }

    // print functions 
    void printHit() const;
    std::string hitInfo() const;
    void printMatchedHits(int) const;
    void printMatchedHits() const;

    // access some of the TLorentzVector functions
    float X() const {return m_hit->X;}
    float Y() const {return m_hit->Y;}
    float Z() const {return m_hit->Z;}
    float T() const {return m_hit->T;}
    float Perp() const {return m_position.Perp();}
    float Phi() const {return m_position.Phi();} // returns angle from [-pi, pi]
    float DeltaPhi(HitCollection&) const;
    TLorentzVector GetPosition() const {return m_position;}

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

    // functions to calculate search windows
    float calculateZWindowForNextLevel(float, float, float, float);
    bool calculateRPhiWindowOutToIn(const float, const float, const float);
    float calculateRPhiWindowInToOut(const float, const float, const float);

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
