#ifndef Barrel_h
#define Barrel_h

#include <iostream>


class Barrel {
   private:
     float fRadius;
     float fModuleLength;
     float fModuleWidth;
 
   public:

     // default constructor
     Barrel() {fRadius=0; fModuleLength=0; fModuleWidth=0;}
  
     // parameterized constructor
     Barrel(float pradius, float pmoduleLength, float pmoduleWidth);
 
     // destructor
     ~Barrel();

     // member functions
     void print();
     float GetRadius();
 };

#endif /* Barrel_h */
