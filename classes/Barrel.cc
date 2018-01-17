#include "classes/Barrel.h"

Barrel::Barrel(float pradius, float pmoduleLength, float pmoduleWidth){
 fRadius=pradius;
 fModuleLength=pmoduleLength;
 fModuleWidth=pmoduleWidth;
}

void Barrel::print(){
  std::cout << "Radius: " << fRadius << std::endl;
  std::cout << "Module Length: " << fModuleLength << std::endl;
  std::cout << "Module Width: " << fModuleWidth << std::endl;
}

float Barrel::GetRadius(){
  return fRadius;
}

Barrel::~Barrel(){
}
