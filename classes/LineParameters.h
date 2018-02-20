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
      m_gradient = 0;
      m_y_intercept = 0;
      m_x_intercept =0;
    }

    // constructor
    LineParameters(std::vector<std::pair<float, float> > coordinates){
      m_coordinates = coordinates; 
    }

    // constructor 
    LineParameters(std::vector<Hit*> hits){
      for(Hit* hit : hits){
        m_coordinates.push_back( std::make_pair(hit->Z, hit->HitRadius ) ); // coordinates in (z, r)
      }
    }

    // functions to calculate line parameters
    void simpleLinearLeastSquaresFit();
    void calculateLineParameters(float x0, float y0, float x1, float y1);

    // accessor functions
    float gradient() const{ return m_gradient;}
    float y_intercept() const { return m_y_intercept;}
    float x_intercept() const { return m_x_intercept;}


};

#endif // LineParameters_h
