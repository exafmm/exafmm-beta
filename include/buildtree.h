#ifndef buildtree_h
#define buildtree_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class BuildTree : public Logger {
 private:
  typedef vec<8,int> ivec8;                                     //!< Vector of 8 integer types
//! Binary tree is used for counting number of bodies with a recursive approach
  struct BinaryTreeNode {
    ivec8            NBODY;                                     //!< Number of descendant bodies
    BinaryTreeNode * LEFT;                                      //!< Pointer to left child
    BinaryTreeNode * RIGHT;                                     //!< Pointer to right child
    BinaryTreeNode * BEGIN;                                     //!< Pointer to begining of memory space
    BinaryTreeNode * END;                                       //!< Pointer to end of memory space
  };

//! Octree is used for building the FMM tree structure as "nodes", then transformed to "cells" data structure
  struct OctreeNode {
    int          BODY;                                          //!< Index offset for first body in node
    int          NBODY;                                         //!< Number of descendant bodies
    int          NNODE;                                         //!< Number of descendant nodes
    OctreeNode * CHILD[8];                                      //!< Pointer to child node
    vec3         X;                                             //!< Coordinate at center
  };

  int          ncrit;                                           //!< Number of bodies per leaf cell
  int          nspawn;                                          //!< Threshold of NBODY for spawning new threads
  int          maxlevel;                                        //!< Maximum level of tree
  B_iter       B0;                                              //!< Iterator of first body
  OctreeNode * N0;                                              //!< Octree root node

 private:
//! Get number of binary tree nodes for a given number of bodies
  inline int getNumBinNode(int n) const {
    if (n <= nspawn) return 1;                                  // If less then threshold, use only one node
    else return 4 * ((n - 1) / nspawn) - 1;                     // Else estimate number of binary tree nodes
  }

//! Get maximum number of binary tree nodes for a given number of bodies
  inline int getMaxBinNode(int n) const {
    return (4 * n) / nspawn;                                    // Conservative estimate of number of binary tree nodes
  }

//! Exclusive scan with offset
  inline ivec8 exclusiveScan(ivec8 input, int offset) const {
    ivec8 output;                                               // Output vector
    for (int i=0; i<8; i++) {                                   // Loop over elements
      output[i] = offset;                                       //  Set value
      offset += input[i];                                       //  Increment offset
    }                                                           // End loop over elements
    return output;                                              // Return output vector
  }

//! Create an octree node
  OctreeNode * makeOctNode(int begin, int end, vec3 X, bool nochild) const {
    OctreeNode * octNode = new OctreeNode();                    // Allocate memory for single node
    octNode->BODY = begin;                                      // Index of first body in node
    octNode->NBODY = end - begin;                               // Number of bodies in node
    octNode->NNODE = 1;                                         // Initialize counter for decendant nodes
    octNode->X = X;                                             // Center position of node
    if (nochild) {                                              // If node has no children
      for (int i=0; i<8; i++) octNode->CHILD[i] = NULL;         //  Initialize pointers to children
    }                                                           // End if for node children
    return octNode;                                             // Return node
  }

//! Count bodies in each octant using binary tree recursion
  void countBodies(Bodies& bodies, int begin, int end, vec3 X, BinaryTreeNode * binNode) {
    assert(getNumBinNode(end - begin) <= binNode->END - binNode->BEGIN + 1);
    if (end - begin <= nspawn) {                                // If number of bodies is less than threshold
      for (int i=0; i<8; i++) binNode->NBODY[i] = 0;            //  Initialize number of bodies in octant
      binNode->LEFT = binNode->RIGHT = NULL;                    //  Initialize pointers to left and right child node
      for (int i=begin; i<end; i++) {                           //  Loop over bodies in node
        vec3 x = bodies[i].X;                                   //   Position of body
        int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to
        binNode->NBODY[octant]++;                               //   Increment body count in octant
      }                                                         //  End loop over bodies in node
    } else {                                                    // Else if number of bodies is larger than threshold
      int mid = (begin + end) / 2;                              //  Split range of bodies in half
      int numLeftNode = getNumBinNode(mid - begin);             //  Number of binary tree nodes on left branch
      int numRightNode = getNumBinNode(end - mid);              //  Number of binary tree nodes on right branch
      assert(numLeftNode + numRightNode <= binNode->END - binNode->BEGIN);
      binNode->LEFT = binNode->BEGIN;                           //  Assign first memory address to left node pointer
      binNode->LEFT->BEGIN = binNode->LEFT + 1;                 //  Assign next memory address to left begin pointer
      binNode->LEFT->END = binNode->LEFT + numLeftNode;         //  Keep track of last memory address used by left
      binNode->RIGHT = binNode->LEFT->END;                      //  Assign that same address to right node pointer
      binNode->RIGHT->BEGIN = binNode->RIGHT + 1;               //  Assign next memory address to right begin pointer
      binNode->RIGHT->END = binNode->RIGHT + numRightNode;      //  Keep track of last memory address used by right
      spawn_tasks {                                             //  Initialize tasks
	spawn_task1(bodies, {                                   //  Spawn new task using Intel TBB or MTHREAD
	    countBodies(bodies, begin, mid, X, binNode->LEFT);  //   Recursive call for left branch
	  });                                                   //  Close lambda expression
	countBodies(bodies, mid, end, X, binNode->RIGHT);       //  Recursive call for right branch
	sync_tasks;                                             //  Synchronize tasks
	binNode->NBODY = binNode->LEFT->NBODY + binNode->RIGHT->NBODY;// Sum contribution from both branches
      }
    }                                                           // End if for number of bodies
  }

//! Sort bodies according to octant (Morton order)
  void moveBodies(Bodies& bodies, Bodies& buffer, int begin, int end,
                  BinaryTreeNode * binNode, ivec8 octantOffset, vec3 X) const {
    if (binNode->LEFT == NULL) {                                // If there are no more child nodes
      for (int i=begin; i<end; i++) {                           //  Loop over bodies
        vec3 x = bodies[i].X;                                   //   Position of body
        int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to`
        buffer[octantOffset[octant]] = bodies[i];               //   Permute bodies out-of-place according to octant
        octantOffset[octant]++;                                 //   Increment body count in octant
      }                                                         //  End loop over bodies
    } else {                                                    // Else if there are child nodes
      int mid = (begin + end) / 2;                              //  Split range of bodies in half
      spawn_tasks {                                             //  Initialize tasks
	spawn_task2(bodies, buffer, {                           //  Spawn new task using Intel TBB or MTHREAD
	    moveBodies(bodies, buffer, begin, mid, binNode->LEFT, octantOffset, X);// Recursive call for left branch
	  });                                                   //  Close lambda expression
	octantOffset += binNode->LEFT->NBODY;                   //  Increment the octant offset for right branch
	moveBodies(bodies, buffer, mid, end, binNode->RIGHT, octantOffset, X);// Recursive call for right branch
	sync_tasks;                                             //  Synchronize tasks
      }
    }                                                           // End if for child existance
  }

//! Build nodes of octree adaptively using a top-down approach based on recursion (uses task based thread parallelism)
  OctreeNode * buildNodes(Bodies& bodies, Bodies& buffer, int begin, int end,
                          BinaryTreeNode * binNode, vec3 X, real_t R0, int level=0, bool direction=false) {
    assert(getMaxBinNode(end - begin) <= binNode->END - binNode->BEGIN);
    if (begin == end) return NULL;                              // If no bodies left, return null pointer
    if (end - begin <= ncrit) {                                 // If number of bodies is less than threshold
      if (direction)                                            //  If direction of data is from bodies to buffer
        for (int i=begin; i<end; i++) buffer[i] = bodies[i];    //   Copy bodies to buffer
      return makeOctNode(begin, end, X, true);                  //  Create an octree node and return it's pointer
    }                                                           // End if for number of bodies
    OctreeNode * octNode = makeOctNode(begin, end, X, false);   // Create an octree node with child nodes
    countBodies(bodies, begin, end, X, binNode);                // Count bodies in each octant using binary recursion
    ivec8 octantOffset = exclusiveScan(binNode->NBODY, begin);  // Exclusive scan to obtain offset from octant count
    moveBodies(bodies, buffer, begin, end, binNode, octantOffset, X);// Sort bodies according to octant
    BinaryTreeNode * binNodeOffset = binNode->BEGIN;            // Initialize pointer offset for binary tree nodes
    spawn_tasks {                                               // Initialize tasks
      for (int i=0; i<8; i++) {                                 // Loop over children
	int maxBinNode = getMaxBinNode(binNode->NBODY[i]);      //  Get maximum number of binary tree nodes
	assert(binNodeOffset + maxBinNode <= binNode->END);
	spawn_task2(buffer, bodies, {                           //  Spawn new task using Intel TBB or MTHREAD
	    vec3 Xchild = X;                                    //   Initialize center position of child node
	    real_t r = R0 / (1 << (level + 1));                 //   Radius of cells for child's level
	    for (int d=0; d<3; d++) {                           //   Loop over dimensions
	      Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);   //    Shift center position to that of child node
	    }                                                   //   End loop over dimensions
	    BinaryTreeNode binNodeChild[1];                     //   Allocate new root for this branch
	    binNodeChild->BEGIN = binNodeOffset;                //   Assign first memory address from offset
	    binNodeChild->END = binNodeOffset + maxBinNode;     //   Keep track of last memory address
	    octNode->CHILD[i] = buildNodes(buffer, bodies,      //   Recursive call for each child
	      octantOffset[i], octantOffset[i] + binNode->NBODY[i],//   Range of bodies is calcuated from octant offset
	      binNodeChild, Xchild, R0, level+1, !direction);   //   Alternate copy direction bodies <-> buffer
	  });                                                   //  Close lambda expression
	binNodeOffset += maxBinNode;                            //  Increment offset for binNode memory address
      }                                                         // End loop over children
      sync_tasks;                                               // Synchronize tasks
    }
    for (int i=0; i<8; i++) {                                   // Loop over children
      if (octNode->CHILD[i]) octNode->NNODE += octNode->CHILD[i]->NNODE;// If child exists increment child node counter
    }                                                           // End loop over chlidren
    return octNode;                                             // Return octree node
  }

//! Get cell index
  long long getKey(vec3 X, vec3 Xmin, real_t diameter, int level) {
    int iX[3] = {0, 0, 0};                                      // Initialize 3-D index
    for (int d=0; d<3; d++) iX[d] = int((X[d] - Xmin[d]) / diameter);// 3-D index
    long long index = ((1 << 3 * level) - 1) / 7;               // Levelwise offset
    for (int l=0; l<level; l++) {                               // Loop over levels
      for (int d=0; d<3; d++) index += (iX[d] & 1) << (3 * l + d); //  Interleave bits into Morton key
      for (int d=0; d<3; d++) iX[d] >>= 1;                      //  Bitshift 3-D index
    }                                                           // End loop over levels
    return index;                                               // Return Morton key
  }

//! Create cell data structure from nodes
  void nodes2cells(OctreeNode * octNode, C_iter C, C_iter C0, C_iter CN, vec3 X0, real_t R0, int level=0, int iparent=0) {
    C->PARENT = iparent;                                        // Index of parent cell
    C->R      = R0 / (1 << level);                              // Cell radius
    C->X      = octNode->X;                                     // Cell center
    C->NBODY  = octNode->NBODY;                                 // Number of decendant bodies
    C->IBODY  = octNode->BODY;                                  // Index of first body in cell
    C->BODY   = B0 + C->IBODY;                                  // Iterator of first body in cell
    C->ICELL  = getKey(C->X, X0-R0, 2*C->R, level);             // Get Morton key
    if (octNode->NNODE == 1) {                                  // If node has no children
      C->ICHILD  = 0;                                           //  Set index of first child cell to zero
      C->NCHILD = 0;                                            //  Number of child cells
      assert(C->NBODY > 0);                                     //  Check for empty leaf cells
      maxlevel = std::max(maxlevel, level);                     //  Update maximum level of tree
    } else {                                                    // Else if node has children
      int nchild = 0;                                           //  Initialize number of child cells
      int octants[8];                                           //  Map of child index to octants (for when nchild < 8)
      for (int i=0; i<8; i++) {                                 //  Loop over octants
        if (octNode->CHILD[i]) {                                //   If child exists for that octant
          octants[nchild] = i;                                  //    Map octant to child index
          nchild++;                                             //    Increment child cell counter
        }                                                       //   End if for child existance
      }                                                         //  End loop over octants
      C_iter Ci = CN;                                           //  CN points to the next free memory address
      C->ICHILD = Ci - C0;                                      //  Set Index of first child cell
      C->NCHILD = nchild;                                       //  Number of child cells
      assert(C->NCHILD > 0);                                    //  Check for childless non-leaf cells
      CN += nchild;                                             //  Increment next free memory address
      spawn_tasks {                                             //  Initialize tasks
	for (int i=0; i<nchild; i++) {                          //  Loop over children
	  int octant = octants[i];                              //   Get octant from child index
	  spawn_task0_if(octNode->NNODE > 1000,                 //   Spawn task if number of sub-nodes is large
	    nodes2cells(octNode->CHILD[octant], Ci, C0, CN, X0, R0, level+1, C-C0));// Recursive call for each child
	  Ci++;                                                 //   Increment cell iterator
	  CN += octNode->CHILD[octant]->NNODE - 1;              //   Increment next free memory address
	}                                                       //  End loop over children
	sync_tasks;                                             //  Synchronize tasks
      }
      for (int i=0; i<nchild; i++) {                            //  Loop over children
        int octant = octants[i];                                //   Get octant from child index
        delete octNode->CHILD[octant];                          //   Free child pointer to avoid memory leak
      }                                                         //  End loop over children
      maxlevel = std::max(maxlevel, level+1);                   //  Update maximum level of tree
    }                                                           // End if for child existance
  }

  //! Transform Xmin & Xmax to X (center) & R (radius)
  Box bounds2box(Bounds bounds) {
    vec3 Xmin = bounds.Xmin;                                    // Set local Xmin
    vec3 Xmax = bounds.Xmax;                                    // Set local Xmax
    Box box;                                                    // Bounding box
    for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2; // Calculate center of domain
    box.R = 0;                                                  // Initialize localRadius
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      box.R = std::max(box.X[d] - Xmin[d], box.R);              //  Calculate min distance from center
      box.R = std::max(Xmax[d] - box.X[d], box.R);              //  Calculate max distance from center
    }                                                           // End loop over dimensions
    box.R *= 1.00001;                                           // Add some leeway to radius
    return box;                                                 // Return box.X and box.R
  }

//! Grow tree structure top down
  void growTree(Bodies &bodies, vec3 X0, real_t R0) {
    assert(R0 > 0);                                             // Check for bounds validity
    Bodies buffer = bodies;                                     // Copy bodies to buffer
    startTimer("Grow tree");                                    // Start timer
    B0 = bodies.begin();                                        // Bodies iterator
    BinaryTreeNode binNode[1];                                  // Allocate root node of binary tree
    int maxBinNode = getMaxBinNode(bodies.size());              // Get maximum size of binary tree
    binNode->BEGIN = new BinaryTreeNode[maxBinNode];            // Allocate array for binary tree nodes
    binNode->END = binNode->BEGIN + maxBinNode;                 // Set end pointer
    N0 = buildNodes(bodies, buffer, 0, bodies.size(), binNode, X0, R0);// Build tree recursively
    delete[] binNode->BEGIN;                                    // Deallocate binary tree array
    stopTimer("Grow tree");                                     // Stop timer
  }

//! Link tree structure
  Cells linkTree(vec3 X0, real_t R0) {
    startTimer("Link tree");                                    // Start timer
    Cells cells;                                                // Initialize cell array
    if (N0 != NULL) {                                           // If he node tree is empty
      cells.resize(N0->NNODE);                                  //  Allocate cells array
      C_iter C0 = cells.begin();                                //  Cell begin iterator
      nodes2cells(N0, C0, C0, C0+1, X0, R0);                    //  Convert nodes to cells recursively
      delete N0;                                                //  Deallocate nodes
    }                                                           // End if for empty node tree
    stopTimer("Link tree");                                     // Stop timer
    return cells;                                               // Return cells array
  }

 public:
  BuildTree(int _ncrit, int _nspawn) : ncrit(_ncrit), nspawn(_nspawn), maxlevel(0) {}

//! Build tree structure top down
  Cells buildTree(Bodies &bodies, Bounds bounds) {
    Box box = bounds2box(bounds);                               // Get box from bounds
    growTree(bodies,box.X,box.R);                               // Grow tree from root
    return linkTree(box.X,box.R);                               // Form parent-child links in tree
  }

//! Print tree structure statistics
  void printTreeData(Cells &cells) {
    if (verbose) {                                              // If verbose flag is true
      printTitle("Tree stats");                                 //  Print title
      std::cout  << std::setw(stringLength) << std::left        //  Set format
	        << "Bodies"     << " : " << cells.front().NBODY << std::endl// Print number of bodies
	        << std::setw(stringLength) << std::left         //  Set format
	        << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
	        << std::setw(stringLength) << std::left         //  Set format
	        << "Tree depth" << " : " << maxlevel << std::endl;//  Print number of levels
    }                                                           // End if for verbose flag
  }
};
#endif
