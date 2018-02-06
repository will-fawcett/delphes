#ifndef myTrack_h
#define myTrack_h

#include "classes/DelphesClasses.h"
#include "classes/lineParameters.h"

// store cartesian coordinate 
struct cartesianCoordinate {
  float x;
  float y;
  float z;
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

    bool isNotFake() const;

};

#endif // myTrack_h
