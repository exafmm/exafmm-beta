#ifndef topdown_h
#define topdown_h
#include "tree.h"

class TopDownTreeConstructor : public TreeStructure {
private:

  struct node {
    typedef std::vector<node>::iterator N_iter;                 // Iterator for node vectors

    int LEVEL;                                                  // Level of node
    int NLEAF;                                                  // Number of leafs in node
    int ISCHILD;                                                // Flag of empty child nodes
    bigint I;                                                   // Morton index
    B_iter LEAF[NCRIT];                                         // Iterator for leafs
    N_iter CHILD[8];                                            // Pointer to child nodes
    vect X;                                                     // Node center
    real R;                                                     // Node radius

    void init(vect x, real r) {
      NLEAF = 0;                                                // Initialize leaf counter
      ISCHILD = 0;                                              // Initialize empty child flag
      X = x;                                                    // Assign node center
      R = r;                                                    // Assign node radius
    }

    void addLeaf(B_iter b) {
      LEAF[NLEAF] = b;                                          // Assign body iterator to leaf
      NLEAF++;                                                  // Increment leaf counter
    }

    int getOctant(vect const pos) {
      int octant(0);                                            // Initialize octant
      for( int d=0; d!=3; ++d )                                 // Loop over dimensions
        octant += (pos[d] > X[d]) << d;                         //  interleave bits and accumulate octant
      return octant;                                            // Return octant
    }

    void addChild(int const i, N_iter &N) {
      int pOffset = ((1 << 3*LEVEL) - 1)/7;                     // Parent Morton offset
      int cOffset = ((1 << 3*(LEVEL+1)) - 1)/7;                 // Current Morton offset
      vect x(X);                                                // Initialize new center position with old center
      real r(R/2);                                              // Initialize new size
      for( int d=0; d!=3; ++d )                                 // Loop over dimensions
        x[d] += r * (((i & 1 << d) >> d) * 2 - 1);              //  Calculate new center position
      CHILD[i] = ++N;                                           // Increment node pointer and assign to child
      CHILD[i]->init(x,r);                                      // Initialize child node
      CHILD[i]->LEVEL=LEVEL+1;                                  // Level of child node
      CHILD[i]->I = ((I-pOffset) << 3) + i + cOffset;           // Morton index of child node
      ISCHILD |= (1 << i);                                      // Flip bit of octant
    }

    void splitNode(N_iter &N) {
      for( int i=0; i!=NCRIT; ++i ) {                           // Loop over all leafs in parent node
        int octant = getOctant(LEAF[i]->pos);                   //  Find the octant where the body belongs
        if( !(ISCHILD & (1 << octant)) )                        //  If child doesn't exist in this octant
          addChild(octant,N);                                   //   Add new child to list
        CHILD[octant]->addLeaf(LEAF[i]);                        //  Add leaf to child
        if( CHILD[octant]->NLEAF >= NCRIT )                     //  If there are still too many leafs
          CHILD[octant]->splitNode(N);                          //   Split the node into smaller ones
      }                                                         // End loop over leafs
    }
  };

  typedef std::vector<node>           Nodes;                    // Vector of nodes
  typedef std::vector<node>::iterator N_iter;                   // Iterator for node vectors
  Nodes nodes;                                                  // Nodes in the tree
  N_iter N0,NN;                                                 // Iterators for nodes

public:
  TopDownTreeConstructor(Bodies &b) : TreeStructure(b){         // Constructor
    nodes.resize(bodies.size());                                // Resize node vector
  }
  ~TopDownTreeConstructor() {}                                  // Destructor


  void traverse(N_iter N, int &boxes, int &leafs, B_iter begin) {
    if( N->NLEAF >= NCRIT ) {
      for( int i=0; i!=8; ++i )
        if( N->ISCHILD & (1 << i) )
          traverse(N->CHILD[i],boxes,leafs,begin);
    } else {
      for( int i=0; i!=N->NLEAF; ++i ) {
        std::cout << boxes << " " << N->I << " " << leafs << " " << N->LEAF[i]-begin << std::endl;
        leafs++;
      }
      boxes++;
    }
  }

  void grow() {
    int octant;                                                 // In which octant is the body located?
    N0 = nodes.begin();                                         // Set iterator to first node
    N0->init(X0,R0);                                            // Initialize root node
    N0->I = 0;                                                  // Morton index of root node
    N0->LEVEL = 0;                                              // Level of root node
    NN = N0;                                                    // Keep copy for node counter
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      N_iter N=N0;                                              //  Reset node iterator
      while( N->NLEAF >= NCRIT ) {                              //  While the nodes have children
        N->NLEAF++;                                             //   Increment the cumulative leaf counter
        octant = N->getOctant(B->pos);                          //   Find the octant where the body belongs
        if( !(N->ISCHILD & (1 << octant)) ) {                   //   If child doesn't exist in this octant
          N->addChild(octant,NN);                               //    Add new child to list
        }                                                       //   Endif for child existence
        N = N->CHILD[octant];                                   //    Update node iterator to child
      }                                                         //  End while loop
      N->addLeaf(B);                                            //  Add body to node as leaf
      if( N->NLEAF >= NCRIT ) {                                 //  If there are too many leafs
        N->splitNode(NN);                                       //   Split the node into smaller ones
      }                                                         //  Endif for splitting
    }                                                           // End loop over all bodies
    int boxes(0),leafs(0);
    traverse(nodes.begin(),boxes,leafs,bodies.begin());
    std::cout << boxes << " "<< leafs << std::endl;
  }

};

#endif
