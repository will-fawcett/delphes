#ifndef myTrack_h
#define myTrack_h

#include "classes/DelphesClasses.h"

// store cartesian coordinate 
struct cartesianCoordinate {
  float x;
  float y;
  float z;
};

// enumeration for different track parameter algorithms 
enum trackParamAlgo{
  beamlineConstraint,
  noBeamlineConstraint,
  MAXIMUM
};




// object to store (or calculate?) track parameters
// what do we want a track to be defined by? 
// we need the track parameters: z0, d0, 
class myTrack{
  private:

    /*float m_gradient;*/
    /*float m_intercept;*/
    std::vector<Hit*> m_associatedHits; 

    // track parameters
    float m_d0;
    float m_z0;
    float m_phi;
    float m_theta;
    float m_eta;
    /*float m_qOverP;*/
    float m_pT;

    bool m_initialised;

    // functions to calculate track parameters 
    bool trackParametersBeamlineConstraint();
    bool trackParametersNoBeamlineConstraint();

  public:



    myTrack(){
      m_d0=0;
      m_z0=0;
      m_phi=0;
      m_theta=0;
      /*m_qOverP=0;*/
      m_pT=0;
      m_initialised = false;
      m_associatedHits = {};
    } // default constructor


    // constructor 
    myTrack(std::vector<Hit*> hits){
      m_associatedHits = hits;
      m_initialised = calculateTrackParameters(); // call default track parameter calculation
    }

    // Calculate the track parameters relative to the detector origin, uses default algorithm
    bool calculateTrackParameters(){
      cartesianCoordinate coordinate;
      coordinate.x = 0.0;
      coordinate.y = 0.0;
      coordinate.z = 0.0;
      m_initialised = calculateTrackParameters(coordinate, beamlineConstraint); // note the default algorithm is beamlineConstraint (!)
    }

    // Calculate track parameters relative to the detector origin, using a specified algorithm
    bool calculateTrackParameters(trackParamAlgo algo){
      cartesianCoordinate coordinate;
      coordinate.x = 0.0;
      coordinate.y = 0.0;
      coordinate.z = 0.0;
      m_initialised = calculateTrackParameters(coordinate, algo);
    }
      
    // Calculate the track parameters relative to some specified coordinate, with a specified algorithm
    bool calculateTrackParameters( cartesianCoordinate, trackParamAlgo ); 

    // accessor functions
    float Pt() const {return m_pT;}
    float Z0() const {return m_z0;}
    float Theta() const {return m_theta;}
    float Phi() const {return m_phi;}
    float D0() const {return m_d0;}
    float Eta() const {return m_eta;} 

    // test if track is a fake track 
    bool isFake() const;
    bool isNotFake() const;

};

#endif // myTrack_h
