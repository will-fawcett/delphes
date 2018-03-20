#include "classes/Location.h"
#include <sstream>
#include <algorithm>


// Function to split string along delimiter, from: https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
std::vector<int> split(const std::string& s, char delimiter)
{
   std::vector<int> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(stoi(token));
   }
   return tokens;
}

// Format the location string
std::string Location::formatLocation(int layerID, int phiBin, int etaBin) const{
  return std::to_string(layerID)+"_"+std::to_string(phiBin)+"_"+std::to_string(etaBin);
}

// Calculate location string from a Candidate object
std::string Location::locationFromHit(Candidate* hit) const{
  return this->locationFromEtaPhi(hit->SurfaceID, hit->Eta, hit->Phi);
}

// Calculate location string from a Hit object
std::string Location::locationFromHit(Hit* hit) const{
  return this->locationFromEtaPhi(hit->SurfaceID, hit->Eta, hit->Phi);
}

// Calculate location string from a surface ID and the eta and phi coordinates 
std::string Location::locationFromEtaPhi(int surfaceID, float eta, float phi) const{

  float hitPhi = phi + M_PI; // convert phi from [-pi, pi] into [0, 2pi]  
  
  int thisPhiBin = floor(hitPhi/(2*M_PI/m_nPhiBins));
  int thisEtaBin = floor(eta/(3.0/m_nEtaBins)); // tolerate negative eta coordinates

  //std::string loc = std::to_string(hit->SurfaceID)+"_"+std::to_string(thisPhiBin)+"_"+std::to_string(thisEtaBin);
  //std::cout << "eta: " << hit->Eta << " phi: " << hit->Phi << " loc: " << loc << std::endl;
  return this->formatLocation(surfaceID, thisPhiBin, thisEtaBin);
}

bool is_in_vector(std::vector<std::string>& v, std::string x){

  if(std::find(v.begin(), v.end(), x) != v.end()) {
    return true;
  } else {
    return false;
  }
}


std::vector< std::string > Location::listOfLocationsInLayer(std::string location, int layerID) {

  // Based on matching criteria defined here, return a list of the locations in layer layerID
  // that are matched

  // matching criteria is the set of squares in eta phi that surround the square given by the location
  
  //std::cout << "Location::listOfLocationsInLayer()" << std::endl;
  //std::cout << "\tInputs: location " << location << " layer ID " << layerID << std::endl;
  
  // split the string
  std::vector<int> tokens = split(location, '_');
  if(tokens.size() != 3){
    std::cerr << "ERROR: unable to split location string" << std::endl;
  }
  int layer = tokens.at(0);
  int phiBin = tokens.at(1);
  int etaBin = tokens.at(2); 

  // locations to return are those in the square around the input location
  std::vector<std::string> newLocations;
  for(int iPhi=phiBin-1; iPhi<=phiBin+1; ++iPhi){
    if(iPhi > m_maxPhiBin || iPhi < m_minPhiBin) continue; // reject combinations outside of range
    for(int iEta=etaBin-1; iEta<=etaBin+1; ++iEta){
      if(iEta > m_maxEtaBin || iEta < m_minEtaBin) continue; // reject combinations outside of range
      std::string newLocation = this->formatLocation(layerID, iPhi, iEta);
      
      // check location is in list of possible locations
      //if(m_hasLocationsSet){
      //if( is_in_vector( m_setOfLocationStrings, newLocation) ){
          // if so add location 
          newLocations.push_back(newLocation);
          //}
          //}

    }
  }

  return newLocations;
}

void Location::addSetOfLocationsStrings(std::vector<std::string> set){
  m_setOfLocationStrings = set; 
  m_hasLocationsSet = true;
}

void Location::printProperties() const{
  std::cout << "Max eta bin" << m_maxEtaBin << std::endl;
  std::cout << "Min eta bin" << m_minEtaBin << std::endl;

  std::cout << "Max phi bin" << m_maxPhiBin << std::endl;
  std::cout << "Min phi bin" << m_minPhiBin << std::endl;

  std::cout << "number of eta bins: " << m_nEtaBins << std::endl;
  std::cout << "number of phi bins: " << m_nPhiBins << std::endl;
}

