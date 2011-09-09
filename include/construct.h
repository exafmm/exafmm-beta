#ifndef construct_h
#define construct_h
#include "tree.h"

class TopDown : virtual public TreeStructure {                  // Topdown tree constructor class
private:
  struct Node {                                                 // Nodes are primitive cells
    int LEVEL;                                                  // Level of node
    int ICHILD;                                                 // Flag of empty child nodes
    int NLEAF;                                                  // Number of leafs in node
    bigint I;                                                   // Cell index
    bigint CHILD[8];                                            // Iterator offset of child nodes
    B_iter LEAF[NCRIT];                                         // Iterator for leafs
    vect X;                                                     // Node center
    real R;                                                     // Node radius
  };
  std::vector<Node> nodes;                                      // Nodes in the tree

  int getOctant(const vect X, int i) {                          // Calculate octant from position
    int octant = 0;                                             // Initialize octant
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      octant += (X[d] > nodes[i].X[d]) << d;                    //  interleave bits and accumulate octant
    }                                                           // End loop over dimensions
    return octant;                                              // Return octant
  }

  void addChild(const int octant, int i) {                      // Add child node and link it
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

  void addLeaf(B_iter b, int i) {                               // Add leaf to node
    nodes[i].LEAF[nodes[i].NLEAF] = b;                          // Assign body iterator to leaf
    nodes[i].NLEAF++;                                           // Increment leaf counter
  }

  void splitNode(int i) {                                       // Split node and reassign leafs to child nodes
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

  void traverse(std::vector<Node>::iterator N) {                // Traverse tree
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
  TopDown() : TreeStructure() {}                                // Constructor
  ~TopDown() {}                                                 // Destructor

  void grow(Bodies &bodies) {                                   // Grow tree from root
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
  }

  void setIndex() {                                             // Store cell index of all bodies
    traverse(nodes.begin());                                    // Traverse tree
  }
};

class BottomUp : virtual public TreeStructure {                 // Bottomup tree constructor
protected:
  int getMaxLevel(Bodies &bodies) {                             // Max level for bottom up tree build
    const long N = bodies.size() * MPISIZE;                     // Number of bodies
    int level;                                                  // Max level
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;   // Decide max level from N/Ncrit
    int MPIlevel = int(log(MPISIZE-1) / M_LN2 / 3) + 1;         // Level of local root cell
    if( MPISIZE == 1 ) MPIlevel = 0;                            // For serial execution local root cell is root cell
    if( MPIlevel > level ) {                                    // If process hierarchy is deeper than tree
      std::cout << "Process hierarchy is deeper than tree @ rank" << MPIRANK << std::endl;
      level = MPIlevel;
    }
    return level;                                               // Return max level
  }

public:
  BottomUp() : TreeStructure() {}                               // Constructor
  ~BottomUp() {}                                                // Destructor

  void setIndex(Bodies &bodies, int level=-1, int begin=0, int end=0, bool update=false) {// Set cell index of all bodies
    bigint i;                                                   // Levelwise cell index
    if( level == -1 ) level = getMaxLevel(bodies);              // Decide max level
    bigint off = ((1 << 3*level) - 1) / 7;                      // Offset for each level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // 3-D cell index
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        nx[d] = int( ( bodies[b].X[d] - (X0[d]-R0) ) / r );     //   3-D cell index
      }                                                         //  End loop over dimension
      i = 0;                                                    //  Initialize cell index
      for( int l=0; l!=level; ++l ) {                           //  Loop over levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d] % 2 << (3 * l + d);                        //    Accumulate cell index
          nx[d] >>= 1;                                          //    Bitshift 3-D cell index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      if( !update ) {                                           //  If this not an update
        bodies[b].ICELL = i+off;                                //   Store index in bodies
      } else if( i+off > bodies[b].ICELL ) {                    //  If the new cell index is larger
        bodies[b].ICELL = i+off;                                //   Store index in bodies
      }                                                         //  Endif for update
    }                                                           // End loop over bodies
  }

  void prune(Bodies &bodies) {                                  // Prune tree by merging cells
    int maxLevel = getMaxLevel(bodies);                         // Max level for bottom up tree build
    for( int l=maxLevel; l>0; --l ) {                           // Loop upwards from bottom level
      int level = getLevel(bodies[0].ICELL);                    //  Current level
      bigint cOff = ((1 << 3 * level) - 1) / 7;                 //  Current ce;; index offset
      bigint pOff = ((1 << 3 * (l-1)) - 1) / 7;                 //  Parent cell index offset
      bigint index = ((bodies[0].ICELL-cOff) >> 3*(level-l+1)) + pOff;// Current cell index
      int begin = 0;                                            //  Begin cell index for bodies in cell
      int size = 0;                                             //  Number of bodies in cell
      int b = 0;                                                //  Current body index
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {//  Loop over bodies
        level = getLevel(B->ICELL);                             //   Level of twig
        cOff = ((1 << 3*level) - 1) / 7;                        //   Offset of twig
        bigint p = ((B->ICELL-cOff) >> 3*(level-l+1)) + pOff;   //   Cell index of parent cell
        if( p != index ) {                                      //   If it's a new parent cell
          if( size < NCRIT ) {                                  //    If parent cell has few enough bodies
            for( int i=begin; i!=begin+size; ++i ) {            //     Loop over bodies in that cell
              bodies[i].ICELL = index;                          //      Renumber cell index to parent cell
            }                                                   //     End loop over bodies in cell
          }                                                     //    Endif for merging
          begin = b;                                            //    Set new begin index
          size = 0;                                             //    Reset number of bodies
          index = p;                                            //    Set new parent cell
        }                                                       //   Endif for new cell
        size++;                                                 //   Increment body counter
      }                                                         //  End loop over bodies
      if( size < NCRIT ) {                                      //  If last parent cell has few enough bodies
        for( int i=begin; i!=begin+size; ++i ) {                //   Loop over bodies in that cell
          bodies[i].ICELL = index;                              //    Renumber cell index to parent cell
        }                                                       //   End loop over bodies in cell
      }                                                         //  Endif for merging
    }                                                           // End loop over levels
  }

  void grow(Bodies &bodies, int level=0, int begin=0, int end=0) {// Grow tree by splitting cells
    bigint index = bodies[begin].ICELL;                         // Initialize cell index
    int off=begin, size=0;                                      // Initialize offset, and size
    if( level == 0 ) level = getMaxLevel(bodies);               // Max level for bottom up tree build
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over bodies under consideration
      if( bodies[b].ICELL != index ) {                          //  If it's a new cell
        if( size >= NCRIT ) {                                   //   If the cell has too many bodies
          level++;                                              //    Increment level
          setIndex(bodies,level,off,off+size);                  //    Set new cell index considering new level
          sortBodies(bodies,buffer,false,off,off+size);         //    Sort new cell index
          grow(bodies,level,off,off+size);                      //    Recursively grow tree
          level--;                                              //    Go back to previous level
        }                                                       //   Endif for splitting
        off = b;                                                //   Set new offset
        size = 0;                                               //   Reset number of bodies
        index = bodies[b].ICELL;                                //   Set new cell
      }                                                         //  Endif for new cell
      size++;                                                   //  Increment body counter
    }                                                           // End loop over bodies
    if( size >= NCRIT ) {                                       // If last cell has too many bodies
      level++;                                                  //  Increment level
      setIndex(bodies,level,off,off+size);                      //  Set new cell index considering new level
      sortBodies(bodies,buffer,false,off,off+size);             //  Sort new cell index
      grow(bodies,level,off,off+size);                          //  Recursively grow tree
      level--;                                                  //  Go back to previous level
    }                                                           // Endif for splitting
  }
};

class TreeConstructor : public TopDown, public BottomUp {       // General tree constructor interface
public:
  TreeConstructor() : TopDown(), BottomUp() {                   // Constructor
    preCalculation();
  }
  ~TreeConstructor() {                                          // Destructor
    postCalculation();
  }

  void topdown(Bodies &bodies, Cells &cells) {                  // Topdown tree constructor interface
    startTimer("Grow tree    ");                                // Start timer
    TopDown::grow(bodies);                                      // Grow tree structure topdown
    stopTimer("Grow tree    ",printNow);                        // Stop timer & print

    startTimer("Set index    ");                                // Start timer
    TopDown::setIndex();                                        // Set index of cells
    stopTimer("Set index    ",printNow);                        // Stop timer & print

    startTimer("Sort bodies  ");                                // Start timer
    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    stopTimer("Sort bodies  ",printNow);                        // Stop timer & print

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes that are not twigs in the current process
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

  void bottomup(Bodies &bodies, Cells &cells) {                 // Bottomup tree constructor interface
    startTimer("Set index    ");                                // Start timer
    BottomUp::setIndex(bodies);                                 // Set index of cells
    stopTimer("Set index    ",printNow);                        // Stop timer & print

    startTimer("Sort bodies  ");                                // Start timer
    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    stopTimer("Sort bodies  ",printNow);                        // Stop timer & print

/*
    startTimer("Prune tree   ");                                // Start timer
    prune(bodies);                                              // Prune tree structure bottomup
    stopTimer("Prune tree   ",printNow);                        // Stop timer & print

    startTimer("Grow tree    ");                                // Start timer
    BottomUp::grow(bodies);                                     // Grow tree structure at bottom if necessary
    stopTimer("Grow tree    ",printNow);                        // Stop timer & print

    startTimer("Sort bodies  ");                                // Start timer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    stopTimer("Sort bodies  ",printNow);                        // Stop timer & print
*/

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes not twigs here
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }
};

#endif
