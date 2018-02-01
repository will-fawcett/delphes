#ifndef HitCollection_h
#define HitCollection_h

#include "classes/DelphesClasses.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>

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

    std::vector< std::vector<Hit*> > makeHitCollection() const;

    Hit * getHit() const {return m_hit; }

    // get hit map (includes "this" hit)
    std::map<int, std::vector<Hit*> > getLayerToHitMap() const;

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

#endif // HitCollection_h
