#include "classes/HitCollection.h"

// -----------------------------------------------------------------------------

std::vector< std::vector<Hit*> > allHitCombinations(const std::vector< std::vector<Hit*> > collection, std::vector<Hit*> hitsInLayer){
  /**************
   * Multiply two sets in different bases, and increase the dimentinoality of the vector space 
   * Say the input "collection" is set A and hitsInLayer is set B
   * Then this function will return the product set AxB 
   * 
   * Parameters:
   *  collection: A vector of collections of hits
   *  hitsInLayer: A vector of the hits in the next layer to be matched 
   *
   * Returns: 
   *  A vector of collections of hits
   * ************/
  bool debug=false;

  if(debug){
    std::cout << "Input vector has " << collection.size() << " pairs" << std::endl;
    std::cout << "Will be matched with " << hitsInLayer.size() << " hits" << std::endl;
    std::cout << "Will then have " << collection.size()*hitsInLayer.size() << "pairs" << std::endl;
  }

  std::vector< std::vector<Hit* > > newCollection;
  for(Hit* newHit : hitsInLayer){
    // for each hit, create a copy of the previous collection, and add the hit to every of the vector or hits
    // The new collection is then the sum of all these collections
    
    //std::vector< std::vector<Hit*> > collectionCopy = collection; 
    for(auto existingGroup : collection){
      existingGroup.push_back(newHit);
      newCollection.push_back(existingGroup);
    }
  }
  if(debug){
    std::cout << "Resulting collection has size " << newCollection.size() << std::endl;
  }
  return newCollection;
}


// -----------------------------------------------------------------------------


std::map<int, std::vector<Hit*>> HitCollection::getLayerToHitMap() const{

  std::map<int, std::vector<Hit*>> layerToHits; 
  for(const auto matchedHit : m_assignedHits){
    layerToHits[matchedHit->SurfaceID].push_back(matchedHit->getHit());
  }
  // add "this" hit
  layerToHits[this->SurfaceID].push_back(this->getHit());

  return layerToHits;
}


// -----------------------------------------------------------------------------


std::vector< std::vector< Hit* > > HitCollection::makeHitCollection() const{
  // From a filled HitCollection object, create a vector of hits that correspond to the hits in each layer
  // There may be many combination, so then there is a vector of these vectors
  // Function should only be called on relevant HitCollections
  //
  // To be honest, I'm not sure if this function should even be a member of HitCollection 
  
  std::map<int, std::vector<Hit*> > layerToHits;

  // loop over pointers to matched HitCollections
  // assign to specific layer
  for(const auto matchedHit : m_assignedHits){
    if( matchedHit->SurfaceID != this->SurfaceID){
      layerToHits[matchedHit->SurfaceID].push_back(matchedHit->getHit());
    }
  }

  // list of unique layers 
  std::vector<int> layers;
  for(auto const& key : layerToHits){
    layers.push_back(key.first);
  }
  std::sort(layers.begin(), layers.end()); // should count in ascending order 

  // Seed the combinations with the first hit
  std::vector< std::vector< Hit* > > allCombinations;
  std::vector<Hit*> firstHit;
  firstHit.push_back(this->getHit());
  allCombinations.push_back(firstHit); 

  // Now return all possible combinations 
  for(int layer : layers){
      allCombinations = allHitCombinations(allCombinations, layerToHits[layer]);
  }

  return allCombinations;
}


// -----------------------------------------------------------------------------


float HitCollection::DeltaPhi(HitCollection& h) const  {
  return acos( cos ( this->Phi() - h.Phi() ));
}


// -----------------------------------------------------------------------------


void HitCollection::printMatchedHits(int level) const{
  std::cout << this->hitInfo() << std::endl; 
  int iHit(1);

  // loop over pointers of subhits  
  for(HitCollection* subHit : m_assignedHits){

    for(int i=0; i<level+1; ++i) std::cout << "\t";
    std::cout << "h" << iHit << " " << subHit->hitInfo() <<  "deltaPhi " << this->DeltaPhi(*subHit) << std::endl;
    iHit++;

    // recursively print out hits belonging to any subhits
    if(subHit->countAssignedHits() > 0){
      subHit->printMatchedHits(level+1);
    }

  }
}

// -----------------------------------------------------------------------------

void HitCollection::printMatchedHits() const{
  this->printMatchedHits(0);
}

// -----------------------------------------------------------------------------

std::string HitCollection::hitInfo() const{
  std::ostringstream s;
  s << "position: (" << m_hit->X << ", " << m_hit->Y << ", " << m_hit->Z << ")"  // cartesian
    << " (" << m_hit->HitRadius << "," << m_hit->Phi << ")" // (r, phi, z) cylindrical polars (z already printed, so not shown here)
    << "\t surface: " << m_hit->SurfaceID
    << "\t has " << this->countAssignedHits() << " hits."
    << " Particle pT=" << dynamic_cast<GenParticle*>(m_hit->Particle.GetObject())->PT << " GeV, "
    << "ref=" << dynamic_cast<GenParticle*>(m_hit->Particle.GetObject())->GetUniqueID() << ". ";
    //<< "ref=" << m_hit->Particle << ". ";
  return s.str();
}

// -----------------------------------------------------------------------------

void HitCollection::printHit() const{
  std::cout << this->hitInfo() << std::endl;
}
