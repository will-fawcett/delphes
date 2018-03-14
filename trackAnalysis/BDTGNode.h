#ifndef BDTGNode__def
#define BDTGNode__def
   
class BDTGNode {
   
public:
   
   // constructor of an essentially "empty" node floating in space
   BDTGNode ( BDTGNode* left,BDTGNode* right,
                          int selector, double cutValue, bool cutType, 
                          int nodeType, double purity, double response ) :
   fLeft         ( left         ),
   fRight        ( right        ),
   fSelector     ( selector     ),
   fCutValue     ( cutValue     ),
   fCutType      ( cutType      ),
   fNodeType     ( nodeType     ),
   fPurity       ( purity       ),
   fResponse     ( response     ){
   }

   virtual ~BDTGNode();

   // test event if it descends the tree at this node to the right
   virtual bool GoesRight( const std::vector<double>& inputValues ) const;
   BDTGNode* GetRight( void )  {return fRight; };

   // test event if it descends the tree at this node to the left 
   virtual bool GoesLeft ( const std::vector<double>& inputValues ) const;
   BDTGNode* GetLeft( void ) { return fLeft; };   

   // return  S/(S+B) (purity) at this node (from  training)

   double GetPurity( void ) const { return fPurity; } 
   // return the node type
   int    GetNodeType( void ) const { return fNodeType; }
   double GetResponse(void) const {return fResponse;}

private:

   BDTGNode*   fLeft;     // pointer to the left daughter node
   BDTGNode*   fRight;    // pointer to the right daughter node
   int                     fSelector; // index of variable used in node selection (decision tree)   
   double                  fCutValue; // cut value applied on this node to discriminate bkg against sig
   bool                    fCutType;  // true: if event variable > cutValue ==> signal , false otherwise
   int                     fNodeType; // Type of node: -1 == Bkg-leaf, 1 == Signal-leaf, 0 = internal 
   double                  fPurity;   // Purity of node from training
   double                  fResponse; // Regression response value of node
}; 
   
//_______________________________________________________________________
   BDTGNode::~BDTGNode()
{
   if (fLeft  != NULL) delete fLeft;
   if (fRight != NULL) delete fRight;
}; 
   
//_______________________________________________________________________
bool BDTGNode::GoesRight( const std::vector<double>& inputValues ) const
{
   // test event if it descends the tree at this node to the right
   bool result;
     result = (inputValues[fSelector] > fCutValue );
   if (fCutType == true) return result; //the cuts are selecting Signal ;
   else return !result;
}
   
//_______________________________________________________________________
bool BDTGNode::GoesLeft( const std::vector<double>& inputValues ) const
{
   // test event if it descends the tree at this node to the left
   if (!this->GoesRight(inputValues)) return true;
   else return false;
}
   
#endif
