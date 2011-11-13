#ifndef topdown_h
#define topdown_h
#include "tree.h"

//! Topdown tree constructor class
template<Equation kernelName>
class TopDown : virtual public TreeStructure<kernelName> {
private:
//! Nodes are primitive cells
  struct Node {
    int LEVEL;                                                  //!< Level of node
    int ICHILD;                                                 //!< Flag of empty child nodes
    int NLEAF;                                                  //!< Number of leafs in node
    bigint I;                                                   //!< Cell index
    bigint CHILD[8];                                            //!< Iterator offset of child nodes
    B_iter LEAF[NCRIT];                                         //!< Iterator for leafs
    vect X;                                                     //!< Node center
    real R;                                                     //!< Node radius
  };
  std::vector<Node> nodes;                                      //!< Nodes in the tree

//! Calculate octant from position
  int getOctant(const vect X, int i) {
    int octant = 0;                                             // Initialize octant
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      octant += (X[d] > nodes[i].X[d]) << d;                    //  interleave bits and accumulate octant
    }                                                           // End loop over dimensions
    return octant;                                              // Return octant
  }

//! Add child node and link it
  void addChild(const int octant, int i) {
    bigint pOff = ((1 << 3* nodes[i].LEVEL   ) - 1) / 7;        // Parent cell index offset
    bigint cOff = ((1 << 3*(nodes[i].LEVEL+1)) - 1) / 7;        // Current cell index offset
    vect x = nodes[i].X;                                        // Initialize new center position with old center
    real r = nodes[i].R/2;                                      // Initialize new size
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      x[d] += r * (((octant & 1 << d) >> d) * 2 - 1);           //  Calculate new center position
    }                                                           // End loop over dimensions
    Node node;                                                  // Node structure
    node.NLEAF = node.ICHILD = 0;                               // Initialize child node counters
    node.X = x;                                                 // Initialize child node center
    node.R = r;                                                 // Initialize child node radius
    node.LEVEL = nodes[i].LEVEL + 1;                            // Level of child node
    node.I = ((nodes[i].I-pOff) << 3) + octant + cOff;          // Cell index of child node
    nodes[i].ICHILD |= (1 << octant);                           // Flip bit of octant
    nodes[i].CHILD[octant] = nodes.end()-nodes.begin();         // Link child node to its parent
    nodes.push_back(node);                                      // Push child node into vector
  }

//! Add leaf to node
  void addLeaf(B_iter b, int i) {
    nodes[i].LEAF[nodes[i].NLEAF] = b;                          // Assign body iterator to leaf
    nodes[i].NLEAF++;                                           // Increment leaf counter
  }

//! Split node and reassign leafs to child nodes
  void splitNode(int i) {
    for( int l=0; l!=NCRIT; ++l ) {                             // Loop over leafs in parent node
      int octant = getOctant(nodes[i].LEAF[l]->X,i);            //  Find the octant where the body belongs
      if( !(nodes[i].ICHILD & (1 << octant)) ) {                //  If child doesn't exist in this octant
        addChild(octant,i);                                     //   Add new child to list
      }                                                         //  Endif for octant
      int c = nodes[i].CHILD[octant];                           //  Set counter for child node
      addLeaf(nodes[i].LEAF[l],c);                              //  Add leaf to child
      if( nodes[c].NLEAF >= NCRIT ) {                           //  If there are still too many leafs
        splitNode(c);                                           //   Split the node into smaller ones
      }                                                         //  Endif for too many leafs
    }                                                           // End loop over leafs
  }

//! Traverse tree
  void traverse(typename std::vector<Node>::iterator N) {
    if( N->NLEAF >= NCRIT ) {                                   // If node has children
      for( int i=0; i!=8; ++i ) {                               // Loop over children
        if( N->ICHILD & (1 << i) ) {                            //  If child exists in this octant
          traverse(nodes.begin()+N->CHILD[i]);                  //   Recursively search child node
        }                                                       //  Endif for octant
      }                                                         // End loop over children
    } else {                                                    //  If child doesn't exist
      for( int i=0; i!=N->NLEAF; ++i ) {                        //   Loop over leafs
        N->LEAF[i]->ICELL = N->I;                               //    Store cell index in bodies
      }                                                         //   End loop over leafs
    }                                                           //  Endif for child existence
  }

public:
//! Constructor
  TopDown() : TreeStructure<kernelName>() {}
//! Destructor
  ~TopDown() {}

//! Grow tree from root
  void grow(Bodies &bodies) {
    this->startTimer("Grow tree    ");                          // Start timer
    int octant;                                                 // In which octant is the body located?
    Node node;                                                  // Node structure
    node.LEVEL = node.NLEAF = node.ICHILD = node.I = 0;         // Initialize root node counters
    node.X = this->X0;                                          // Initialize root node center
    node.R = this->R0;                                          // Initialize root node radius
    nodes.push_back(node);                                      // Push child node into vector
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int i = 0;                                                //  Reset node counter
      while( nodes[i].NLEAF >= NCRIT ) {                        //  While the nodes have children
        nodes[i].NLEAF++;                                       //   Increment the cumulative leaf counter
        octant = getOctant(B->X,i);                             //   Find the octant where the body belongs
        if( !(nodes[i].ICHILD & (1 << octant)) ) {              //   If child doesn't exist in this octant
          addChild(octant,i);                                   //    Add new child to list
        }                                                       //   Endif for child existence
        i = nodes[i].CHILD[octant];                             //    Update node iterator to child
      }                                                         //  End while loop
      addLeaf(B,i);                                             //  Add body to node as leaf
      if( nodes[i].NLEAF >= NCRIT ) {                           //  If there are too many leafs
        splitNode(i);                                           //   Split the node into smaller ones
      }                                                         //  Endif for splitting
    }                                                           // End loop over bodies
    this->stopTimer("Grow tree    ",this->printNow);            // Stop timer
  }

//! Store cell index of all bodies
  void setIndex() {
    this->startTimer("Set index    ");                          // Start timer
    traverse(nodes.begin());                                    // Traverse tree
    this->stopTimer("Set index    ",this->printNow);            // Stop timer 
  }
};

#endif
