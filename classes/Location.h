#ifndef Location_h
#define Location_h

#include "classes/DelphesClasses.h"
#include <string>
#include <iostream>

/*#define M_PI 3.14159265359*/

class Location{

  private: 
    std::vector< std::string > m_setOfLocationStrings; 
    bool m_hasLocationsSet;

    float m_phiWindowSize;
    float m_etaWindowSize;
    
    int m_nPhiBins;
    int m_nEtaBins;

    // Eta ranges from -inf, inf. Practically however this is from -2.5, 2.5. Extended to -3, 3 for safety 
    const float m_etaMin = -3.0;
    const float m_etaMax = 3.0;

    // Phi ranges from -pi, pi
    const float m_phiMin = -M_PI; 
    const float m_phiMax = M_PI;

    // maximum eta and phi bin numbers
    int m_minPhiBin, m_maxPhiBin, m_minEtaBin, m_maxEtaBin;

    std::string formatLocation(int, int, int) const;


  public:

    Location(){
      m_phiWindowSize = 0.06; 
      m_etaWindowSize = 0.1;  

      m_nPhiBins = 11;
      m_nEtaBins = 60; 
      m_hasLocationsSet = false;
    }

    Location(float phiWindowSize, float etaWindowSize){
      m_phiWindowSize = phiWindowSize;
      m_etaWindowSize = etaWindowSize;

      m_nPhiBins = ceil(2*M_PI / phiWindowSize); 
      m_nEtaBins = ceil( (m_etaMax + fabs(m_etaMin)) / m_etaWindowSize);

      m_minPhiBin = 0; 
      m_maxPhiBin = m_nPhiBins;

      m_minEtaBin = -1*m_nEtaBins;
      m_maxEtaBin = m_nEtaBins;

      m_hasLocationsSet = false; 

    }

    std::string locationFromHit(Hit*) const;
    std::string locationFromEtaPhi(int, float, float) const; 

    std::vector<std::string> listOfLocationsInLayer(std::string, int) ;

    void addSetOfLocationsStrings(std::vector<std::string> ); 
    void printProperties() const;

};

#endif // Location_h
