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

// store parameters of a line (could be put inside TrackFitter class?)
struct lineParameters {
  float gradient;
  float intercept;
};

// store cartesian coordinate 
struct cartesianCoordinate {
  float x;
  float y;
  float z;
};


// Class to store essentially a vector of hits, but the hits know about their relationship to one another
class HitCollection{
  private:
    Hit * m_hit; 
    TLorentzVector m_position;
    std::vector<HitCollection*> m_assignedHits; 
    bool m_debug;

  public:
    int SurfaceID;

    // constructor
    HitCollection(Hit * hitIn){
      m_hit = hitIn; 
      m_position = hitIn->Position();
      SurfaceID = hitIn->SurfaceID;
      m_debug = false;
    };

    // default constructor
    HitCollection(){
      m_hit = 0;
      m_position = TLorentzVector();
      SurfaceID = -1; 
      m_debug = false;
    };

    void printAssignedHitPointers(){
      for(HitCollection * h : m_assignedHits){
        std::cout << h << std::endl;
      }
    }

    void debug(){m_debug=true;}

    // other hits matched to this hit 
    void addHit(HitCollection* collectionIn){ m_assignedHits.push_back( collectionIn ); }

    int countAssignedHits() const { return m_assignedHits.size(); }

    std::vector< std::vector<Hit*> > makeHitCollection();

    Hit * getHit() const {return m_hit; }

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

// object to store (or calculate?) track parameters
// what do we want a track to be defined by? 
// we need the track parameters: z0, d0, 
class myTrack{
  private:

    float m_gradient;
    float m_intercept;
    std::vector<Hit*> m_associatedHits; 

  public:

    // track parameters
    float d0;
    float z0;
    float phi;
    float theta;
    float qOverP;

    myTrack(){
      m_gradient=0;
      m_intercept=0;
      d0=0;
      z0=0;
      phi=0;
      theta=0;
      qOverP=0;
    } // default constructor


    myTrack(lineParameters params, std::vector<Hit*> hits){
      m_gradient = params.gradient;
      m_intercept = params.intercept;
      m_associatedHits = hits;
      calculateTrackParameters();
    }

    void calculateTrackParameters(){
      // calculate the track parameters relative to the detector origin
      cartesianCoordinate coordinate;
      coordinate.x = 0.0;
      coordinate.y = 0.0;
      coordinate.z = 0.0;
      calculateTrackParameters(coordinate);
    }
      
    // Calculate the track parameters relative to some specified coordinate
    void calculateTrackParameters( cartesianCoordinate ); 

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

    // algorithms to associate hits 
    bool associateHitsLinearOutToIn(hitContainer, float, float);
    bool associateHitsLinearInToOut(hitContainer, float, float);

    bool combineHitsToTracksInToOut(); 
    bool combineHitsToTracksOutToIn(); 

    // functions to calculate search windows
    float calculateZWindowForNextLevel(float, float, float, float);
    bool calculateRPhiWindowOutToIn(const float, const float, const float);
    float calculateRPhiWindowInToOut(const float, const float, const float);

    // track fitting functions
    myTrack simpleLinearLeastSquaresFit(std::vector< Hit* > hits) const;
    lineParameters simpleLinearLeastSquaresFit(std::vector< std::pair< float, float> >) const; 

  public:
    
    // constructor 
    TrackFitter(const fitTypes ftIn, std::vector<float> paramIn){
      fitType=ftIn;
      m_parameters=paramIn;
      m_debug = false;
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
#endif // anaClasses_h
