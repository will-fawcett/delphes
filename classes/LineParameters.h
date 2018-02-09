#ifndef LineParameters_h
#define LineParameters_h

#include "classes/DelphesClasses.h"

// store parameters of a line 
class LineParameters {
  
  private:
    std::vector<std::pair<float, float> > m_coordinates;
    float m_gradient;
    float m_y_intercept;
    float m_x_intercept;

  public: 

    // default constructor
    LineParameters(){
      gradient = 0;
      y_intercept = 0;
      x_intercept =0;
    }

    // constructor
    LineParameters(std::vector<std::pair<float, float> > coordinates){
      m_coordinates = coordinates; 
    }

    // constructor 
    LineParameters(std::vector<Hit*> hits){
      for(Hit* hit : hits){
        m_coordinates.push_back( std::make_pair(hit->Z, hit->Perp()) ); // coordinates in (r, z)
      }
    }

    void simpleLinearLeastSquaresFit();



};

#endif // lineParameters_h
