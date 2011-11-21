/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef topdown_h
#define topdown_h
#include "tree.h"

//! Topdown tree constructor class
template<Equation equation>
class TopDown : public TreeStructure<equation> {
public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell

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
  TopDown() : TreeStructure<equation>() {}
//! Destructor
  ~TopDown() {}

//! Grow tree from root
  void grow(Bodies &bodies) {
    startTimer("Grow tree    ");                                // Start timer
    int octant;                                                 // In which octant is the body located?
    Node node;                                                  // Node structure
    node.LEVEL = node.NLEAF = node.ICHILD = node.I = 0;         // Initialize root node counters
    node.X = X0;                                                // Initialize root node center
    node.R = R0;                                                // Initialize root node radius
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
    stopTimer("Grow tree    ",printNow);                        // Stop timer
  }

//! Store cell index of all bodies
  void setIndex() {
    startTimer("Set index    ");                                // Start timer
    traverse(nodes.begin());                                    // Traverse tree
    stopTimer("Set index    ",printNow);                        // Stop timer 
  }
};

#endif
