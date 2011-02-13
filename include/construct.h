#ifndef construct_h
#define construct_h
#include "tree.h"

class TopDown : virtual public TreeStructure {
private:
  struct node {                                                 // Nodes are primitive cells
    typedef std::vector<node>::iterator N_iter;                 // Iterator for node vectors

    int LEVEL;                                                  // Level of node
    int NLEAF;                                                  // Number of leafs in node
    int ICHILD;                                                 // Flag of empty child nodes
    bigint I;                                                   // Cell index
    B_iter LEAF[NCRIT];                                         // Iterator for leafs
    N_iter CHILD[8];                                            // Pointer to child nodes
    vect X;                                                     // Node center
    real R;                                                     // Node radius

    void init(vect x, real r) {                                 // Initialize node information
      NLEAF = 0;                                                // Initialize leaf counter
      ICHILD = 0;                                               // Initialize empty child flag
      X = x;                                                    // Assign node center
      R = r;                                                    // Assign node radius
    }

    void addLeaf(B_iter b) {                                    // Add leaf to node
      LEAF[NLEAF] = b;                                          // Assign body iterator to leaf
      NLEAF++;                                                  // Increment leaf counter
    }

    int getOctant(vect const pos) {                             // Calculate octant from position
      int octant(0);                                            // Initialize octant
      for( int d=0; d!=3; ++d ) {                               // Loop over dimensions
        octant += (pos[d] > X[d]) << d;                         //  interleave bits and accumulate octant
      }                                                         // End loop over dimensions
      return octant;                                            // Return octant
    }

    void addChild(int const i, N_iter &N) {                     // Add child node and link it
      bigint pOff = ((1 << 3* LEVEL   ) - 1) / 7;               // Parent cell index offset
      bigint cOff = ((1 << 3*(LEVEL+1)) - 1) / 7;               // Current cell index offset
      vect x = X;                                               // Initialize new center position with old center
      real r = R/2;                                             // Initialize new size
      for( int d=0; d!=3; ++d ) {                               // Loop over dimensions
        x[d] += r * (((i & 1 << d) >> d) * 2 - 1);              //  Calculate new center position
      }                                                         // End loop over dimensions
      CHILD[i] = ++N;                                           // Increment node pointer and assign to child
      CHILD[i]->init(x,r);                                      // Initialize child node
      CHILD[i]->LEVEL = LEVEL + 1;                              // Level of child node
      CHILD[i]->I = ((I-pOff) << 3) + i + cOff;                 // Cell index of child node
      ICHILD |= (1 << i);                                       // Flip bit of octant
    }

    void splitNode(N_iter &N) {                                 // Split node and reassign leafs to child nodes
      for( int i=0; i!=NCRIT; ++i ) {                           // Loop over all leafs in parent node
        int octant = getOctant(LEAF[i]->pos);                   //  Find the octant where the body belongs
        if( !(ICHILD & (1 << octant)) ) {                       //  If child doesn't exist in this octant
          addChild(octant,N);                                   //   Add new child to list
        }                                                       //  Endif for octant
        CHILD[octant]->addLeaf(LEAF[i]);                        //  Add leaf to child
        if( CHILD[octant]->NLEAF >= NCRIT ) {                   //  If there are still too many leafs
          CHILD[octant]->splitNode(N);                          //   Split the node into smaller ones
        }                                                       //  Endif for too many leafs
      }                                                         // End loop over leafs
    }
  };

  typedef std::vector<node>           Nodes;                    // Vector of nodes
  typedef std::vector<node>::iterator N_iter;                   // Iterator for node vectors
  Nodes nodes;                                                  // Nodes in the tree
  N_iter N0,NN;                                                 // Iterators for nodes

public:
  TopDown(Bodies &b) : TreeStructure(b){                        // Constructor
    nodes.resize(bodies.size()/NCRIT*8);                        // Resize node vector
  }
  ~TopDown() {}                                                 // Destructor

  void grow() {                                                 // Grow tree from root
    int octant;                                                 // In which octant is the body located?
    N0 = nodes.begin();                                         // Set iterator to first node
    N0->init(X0,R0);                                            // Initialize root node
    N0->I = 0;                                                  // Cell index of root node
    N0->LEVEL = 0;                                              // Level of root node
    NN = N0;                                                    // Keep copy for node counter
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      N_iter N=N0;                                              //  Reset node iterator
      while( N->NLEAF >= NCRIT ) {                              //  While the nodes have children
        N->NLEAF++;                                             //   Increment the cumulative leaf counter
        octant = N->getOctant(B->pos);                          //   Find the octant where the body belongs
        if( !(N->ICHILD & (1 << octant)) ) {                    //   If child doesn't exist in this octant
          N->addChild(octant,NN);                               //    Add new child to list
        }                                                       //   Endif for child existence
        N = N->CHILD[octant];                                   //    Update node iterator to child
      }                                                         //  End while loop
      N->addLeaf(B);                                            //  Add body to node as leaf
      if( N->NLEAF >= NCRIT ) {                                 //  If there are too many leafs
        N->splitNode(NN);                                       //   Split the node into smaller ones
      }                                                         //  Endif for splitting
    }                                                           // End loop over all bodies
  }

  void traverse(N_iter N) {                                     // Traverse tree
    if( N->NLEAF >= NCRIT ) {                                   // If node has children
      for( int i=0; i!=8; ++i ) {                               // Loop over children
        if( N->ICHILD & (1 << i) ) {                            //  If child exists in this octant
          traverse(N->CHILD[i]);                                //   Recursively search child node
        }                                                       //  Endif for octant
      }                                                         // End loop over children
    } else {                                                    //  If child doesn't exist
      for( int i=0; i!=N->NLEAF; ++i ) {                        //   Loop over leafs
        N->LEAF[i]->I = N->I;                                   //    Store cell index in bodies
      }                                                         //   End loop over leafs
    }                                                           //  Endif for child existence
  }

  void setIndex() {                                             // Store cell index of all bodies
    traverse(nodes.begin());                                    // Traverse tree
  }
};

class BottomUp : virtual public TreeStructure {
public:
  BottomUp(Bodies &b) : TreeStructure(b){}                      // Constructor
  ~BottomUp() {}                                                // Destructor

  int getMaxLevel() {                                           // Max level for bottom up tree build
    int const N = bodies.size();                                // Number of bodies
    int level;                                                  // Define max level
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;   // Decide max level from N/Ncrit
    return level;                                               // Return max level
  }

  void setIndex(int level=0, int begin=0, int end=0 ) {         // Set cell index of all bodies
    bigint i;                                                   // Levelwise cell index
    if( level == 0 ) level = getMaxLevel();                     // Decide max level
    bigint off = ((1 << 3*level) - 1) / 7;                      // Offset for each level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // Define 3-D cell index
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        nx[d] = int( ( bodies[b].pos[d] - (X0[d]-R0) ) / r );   //   3-D cell index
      }                                                         //  End loop over dimension
      i = 0;                                                    //  Initialize cell index
      for( int l=0; l!=level; ++l ) {                           //  Loop over all levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d] % 2 << (3 * l + d);                        //    Accumulate cell index
          nx[d] >>= 1;                                          //    Bitshift 3-D cell index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      bodies[b].I = i+off;                                      //  Store index in bodies
    }                                                           // End loop over all bodies
  }

  void prune() {                                                // Prune tree by merging cells
    int maxLevel = getMaxLevel();                               // Max level for bottom up tree build
    for( int l=maxLevel; l>0; --l ) {                           // Loop upwards from bottom level
      int level = getLevel(bodies[0].I);                        //  Current level
      bigint cOff = ((1 << 3 * level) - 1) / 7;                 //  Current ce;; index offset
      bigint pOff = ((1 << 3 * (l-1)) - 1) / 7;                 //  Parent cell index offset
      bigint index = ((bodies[0].I-cOff) >> 3*(level-l+1)) + pOff;//  Current cell index
      int begin = 0;                                            //  Begin cell index for bodies in cell
      int size = 0;                                             //  Number of bodies in cell
      int b = 0;                                                //  Current body index
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {//  Loop over all bodies
        level = getLevel(B->I);                                 //   Level of twig
        cOff = ((1 << 3*level) - 1) / 7;                        //   Offset of twig
        bigint p = ((B->I-cOff) >> 3*(level-l+1)) + pOff;       //   Cell index of parent cell
        if( p != index ) {                                      //   If it's a new parent cell
          if( size < NCRIT ) {                                  //    If parent cell has few enough bodies
            for( int i=begin; i!=begin+size; ++i ) {            //     Loop over all bodies in that cell
              bodies[i].I = index;                              //      Renumber cell index to parent cell
            }                                                   //     End loop over all bodies in cell
          }                                                     //    Endif for merging
          begin = b;                                            //    Set new begin index
          size = 0;                                             //    Reset number of bodies
          index = p;                                            //    Set new parent cell
        }                                                       //   Endif for new cell
        size++;                                                 //   Increment body counter
      }                                                         //  End loop over all bodies
      if( size < NCRIT ) {                                      //  If last parent cell has few enough bodies
        for( int i=begin; i!=begin+size; ++i ) {                //   Loop over all bodies in that cell
          bodies[i].I = index;                                  //    Renumber cell index to parent cell
        }                                                       //   End loop over all bodies in cell
      }                                                         //  Endif for merging
    }                                                           // End loop over levels
  }

  void grow(int level=0, int begin=0, int end=0) {              // Grow tree by splitting cells
    bigint index = bodies[begin].I;                             // Initialize cell index
    int off=begin, size=0;                                      // Initialize offset, and size
    if( level == 0 ) level = getMaxLevel();                     // Max level for bottom up tree build
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over all bodies under consideration
      if( bodies[b].I != index ) {                              //  If it's a new cell
        if( size >= NCRIT ) {                                   //   If the cell has too many bodies
          level++;                                              //    Increment level
          setIndex(level,off,off+size);                         //    Set new cell index considering new level
          sort(bodies,buffer,true,off,off+size);                //    Sort new cell index
          grow(level,off,off+size);                             //    Recursively grow tree
          level--;                                              //    Go back to previous level
        }                                                       //   Endif for splitting
        off = b;                                                //   Set new offset
        size = 0;                                               //   Reset number of bodies
        index = bodies[b].I;                                    //   Set new cell
      }                                                         //  Endif for new cell
      size++;                                                   //  Increment body counter
    }                                                           // End loop over bodies
    if( size >= NCRIT ) {                                       // If last cell has too many bodies
      level++;                                                  //  Increment level
      setIndex(level,off,off+size);                             //  Set new cell index considering new level
      sort(bodies,buffer,true,off,off+size);                    //  Sort new cell index
      grow(level,off,off+size);                                 //  Recursively grow tree
      level--;                                                  //  Go back to previous level
    }                                                           // Endif for splitting
  }
};

class TreeConstructor : public TopDown, public BottomUp {
public:
  TreeConstructor(Bodies &b) : TreeStructure(b), TopDown(b), BottomUp(b) {}// Constructor
  ~TreeConstructor() {}                                         // Destructor

  void topdown(bool print=true) {
    double tic,toc;
    tic = get_time();
    TopDown::grow();
    toc = get_time();
    if(print) std::cout << "Grow tree     : " << toc-tic << std::endl;

    tic = get_time();
    TopDown::setIndex();
    toc = get_time();
    if(print) std::cout << "Set index     : " << toc-tic << std::endl;

    tic = get_time();
    sort(bodies,buffer);
    toc = get_time();
    if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

    tic = get_time();
    bodies2twigs();
    toc = get_time();
    if(print) std::cout << "Bodies2twigs  : " << toc-tic << std::endl;

    tic = get_time();
    twigs2cells();
    toc = get_time();
    if(print) std::cout << "Twigs2cells   : " << toc-tic << std::endl;
  }

  void bottomup(bool print=true) {
    double tic,toc;
    tic = get_time();
    BottomUp::setIndex();
    toc = get_time();
    if(print) std::cout << "Set index     : " << toc-tic << std::endl;

    tic = get_time();
    sort(bodies,buffer);
    toc = get_time();
    if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

    tic = get_time();
    prune();
    toc = get_time();
    if(print) std::cout << "Prune tree    : " << toc-tic << std::endl;

    tic = get_time();
    BottomUp::grow();
    toc = get_time();
    if(print) std::cout << "Grow tree     : " << toc-tic << std::endl;

    tic = get_time();
    sort(bodies,buffer);
    toc = get_time();
    if(print) std::cout << "Sort descend  : " << toc-tic << std::endl;

    tic = get_time();
    bodies2twigs();
    toc = get_time();
    if(print) std::cout << "Bodies2twigs  : " << toc-tic << std::endl;

    tic = get_time();
    twigs2cells();
    toc = get_time();
    if(print) std::cout << "Twigs2cells   : " << toc-tic << std::endl;
  }

};

#endif
