#include "classes/LineParameters.h"

void LineParameters::calculateLineParameters(float x0, float y0, float x1, float y1){
  // Calculate the parameters of the straight line passing through the coordinates (x0, y0), (x1, y1)

  m_gradient = (y0 - y1) / (x0 - x1); 
  m_y_intercept  = (x1*y0 - x0*y1) / (x1 - x0);
  //m_y_intercept = y1 - m_gradient*x1;
  m_x_intercept = -1*m_y_intercept/m_gradient; 
}


void LineParameters::simpleLinearLeastSquaresFit() {
  // Function to do simple linear least squares fitting (for a straight line)
  // with the parameters y = mx + c. 
  // Extracts the best fit for m and c. 
  // Takes a vector of the coordinates {xi, yi} 

  float X(0), Y(0), XX(0), XY(0);
  float n = m_coordinates.size();

  for(const auto& coord : m_coordinates){
    float xi = coord.first;
    float yi = coord.second;
    X += xi;
    Y += yi;
    XX += xi*xi;
    XY += xi*yi;
  }

  // gradient
  m_gradient = (XY*n - X*Y) / (XX*n -X*X);

  // y-intercept
  m_y_intercept = (Y - m_gradient*X) / n; 
  
  // x-intercept
  m_x_intercept = -1*m_y_intercept/m_gradient;

}

