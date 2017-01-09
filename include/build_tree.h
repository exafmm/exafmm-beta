#ifndef build_tree_tbb_h
#define build_tree_tbb_h
#include "logger.h"
#include "thread.h"
#include "types.h"

namespace exafmm {
  template<typename Kernel>
  class BuildTree {
    typedef typename Kernel::Bodies Bodies;                     //!< Vector of bodies
    typedef typename Kernel::Cells Cells;                       //!< Vector of cells
    typedef typename Kernel::B_iter B_iter;                     //!< Iterator of body vector
    typedef typename Kernel::C_iter C_iter;                     //!< Iterator of cell vector
    typedef vec<8,int> ivec8;                                   //!< Vector of 8 integer types

  private:
    //! Binary tree is used for counting number of bodies with a recursive approach
    struct BinaryTreeNode {
      ivec8            NBODY;                                   //!< Number of descendant bodies
      BinaryTreeNode * LEFT;                                    //!< Pointer to left child
      BinaryTreeNode * RIGHT;                                   //!< Pointer to right child
      BinaryTreeNode * BEGIN;                                   //!< Pointer to begining of memory space
      BinaryTreeNode * END;                                     //!< Pointer to end of memory space
    };

    //! Octree is used for building the FMM tree structure as "nodes", then transformed to "cells" data structure
    struct OctreeNode {
      int          IBODY;                                       //!< Index offset for first body in node
      int          NBODY;                                       //!< Number of descendant bodies
      int          NNODE;                                       //!< Number of descendant nodes
      OctreeNode * CHILD[8];                                    //!< Pointer to child node
      vec3         X;                                           //!< Coordinate at center
    };

    const int    ncrit;                                         //!< Number of bodies per leaf cell
    const int    nspawn;                                        //!< Threshold of NBODY for spawning new threads
    int          numLevels;                                     //!< Number of levels in tree
    B_iter       B0;                                            //!< Iterator of first body
    OctreeNode * N0;                                            //!< Pointer to octree root node

  private:
    //! Counting bodies in each octant
    static void countBodies(Bodies & bodies, int begin, int end, vec3 X,
                     BinaryTreeNode * binNode) {
      for (int i=0; i<8; i++) binNode->NBODY[i] = 0;            // Initialize number of bodies in octant
      binNode->LEFT = binNode->RIGHT = NULL;                    // Initialize pointers to left and right child node
      for (int i=begin; i<end; i++) {                           // Loop over bodies in node
        vec3 x = bodies[i].X;                                   //  Coordinates of body
        if (bodies[i].ICELL < 0)                                //  If using residual index
          x = bodies[i+bodies[i].ICELL].X;                      //   Use coordinates of first body in residual group
        int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to
        binNode->NBODY[octant]++;                               //  Increment body count in octant
      }                                                         // End loop over bodies in node
    }

    //! Sorting bodies according to octant (Morton order)
    static void moveBodies(Bodies & bodies, Bodies & buffer, int begin, int end,
                    BinaryTreeNode * binNode, ivec8 octantOffset, vec3 X) {
      for (int i=begin; i<end; i++) {                           // Loop over bodies
        vec3 x = bodies[i].X;                                   //  Coordinates of body
        if (bodies[i].ICELL < 0)                                //  If using residual index
          x = bodies[i+bodies[i].ICELL].X;                      //   Use coordinates of first body in residual group
        int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to`
        buffer[octantOffset[octant]] = bodies[i];               //   Permute bodies out-of-place according to octant
        octantOffset[octant]++;                                 //  Increment body count in octant
      }                                                         // End loop over bodies
    }

    //! Create an octree node
    static OctreeNode * makeOctNode(OctreeNode *& octNode, int begin, int end, vec3 X, bool nochild) {
      octNode = new OctreeNode();                               // Allocate memory for single node
      octNode->IBODY = begin;                                   // Index of first body in node
      octNode->NBODY = end - begin;                             // Number of bodies in node
      octNode->NNODE = 1;                                       // Initialize counter for decendant nodes
      octNode->X = X;                                           // Center coordinates of node
      if (nochild) {                                            // If node has no children
        for (int i=0; i<8; i++) octNode->CHILD[i] = NULL;       //  Initialize pointers to children
      }                                                         // End if for node children
      return octNode;                                           // Return node
    }

    //! Exclusive scan with offset
    static inline ivec8 exclusiveScan(ivec8 input, int offset) {
      ivec8 output;                                             // Output vector
      for (int i=0; i<8; i++) {                                 // Loop over elements
        output[i] = offset;                                     //  Set value
        offset += input[i];                                     //  Increment offset
      }                                                         // End loop over elements
      return output;                                            // Return output vector
    }

    //! Get maximum number of binary tree nodes for a given number of bodies
    static inline int getMaxBinNode(int n) {
      return (4 * n) / 5000;                                  // Conservative estimate of number of binary tree nodes
    }

    //! Recursive functor for building nodes of an octree adaptively using a top-down approach
    static void buildNodes(OctreeNode *& octNode, Bodies & bodies, Bodies & buffer,
                           int begin, int end, BinaryTreeNode * binNode,
                           vec3 X, real_t R0, int ncrit,
                           int level=0, bool direction=false) {
      assert(getMaxBinNode(end - begin) <= binNode->END - binNode->BEGIN);// Bounds checking for node range
      if (begin == end) {                                       // If no bodies are left
        octNode = NULL;                                         //  Assign null pointer
        return;                                                 //  End buildNodes()
      }                                                         // End if for no bodies
      if (end - begin <= ncrit) {                               // If number of bodies is less than threshold
        if (direction)                                          //  If direction of data is from bodies to buffer
          for (int i=begin; i<end; i++) buffer[i] = bodies[i];  //   Copy bodies to buffer
        octNode = makeOctNode(octNode,begin,end,X,true);        //  Create an octree node and assign it's pointer
        return;                                                 //  End buildNodes()
      }                                                         // End if for number of bodies
      octNode = makeOctNode(octNode,begin,end,X,false);         // Create an octree node with child nodes
      countBodies(bodies, begin, end, X, binNode);              // Count bodies in each octant
      ivec8 octantOffset = exclusiveScan(binNode->NBODY, begin);// Exclusive scan to obtain offset from octant count
      moveBodies(bodies, buffer, begin, end, binNode, octantOffset, X);// Sort bodies according to octant
      BinaryTreeNode * binNodeOffset = binNode->BEGIN;          // Initialize pointer offset for binary tree nodes
      BinaryTreeNode binNodeChild[8];                           // Allocate new root for this branch
      for (int i=0; i<8; i++) {                                 // Loop over children
        int maxBinNode = getMaxBinNode(binNode->NBODY[i]);      //  Get maximum number of binary tree nodes
        assert(binNodeOffset + maxBinNode <= binNode->END);     //  Bounds checking for node count
        vec3 Xchild = X;                                        //  Initialize center coordinates of child node
        real_t r = R0 / (1 << (level + 1));                     //  Radius of cells for child's level
        for (int d=0; d<3; d++) {                               //  Loop over dimensions
          Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);       //   Shift center coordinates to that of child node
        }                                                       //  End loop over dimensions
        binNodeChild[i].BEGIN = binNodeOffset;                  //  Assign first memory address from offset
        binNodeChild[i].END = binNodeOffset + maxBinNode;       //  Keep track of last memory address
        buildNodes(octNode->CHILD[i], buffer, bodies,           //  Recursive call for children
                   octantOffset[i], octantOffset[i] + binNode->NBODY[i],
                   &binNodeChild[i], Xchild, R0, ncrit, level+1, !direction);
        binNodeOffset += maxBinNode;                            //  Increment offset for binNode memory address
      }                                                         // End loop over children
      for (int i=0; i<8; i++) {                                 // Loop over children
        if (octNode->CHILD[i]) octNode->NNODE += octNode->CHILD[i]->NNODE;// If child exists increment child node count
      }                                                         // End loop over chlidren
    }

    //! Recursive functor for creating cell data structure from nodes
    struct Nodes2cells {
      OctreeNode * octNode;                                     //!< Pointer to octree node
      B_iter B0;                                                //!< Iterator of first body
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      mutable C_iter CN;                                        //!< Iterator of cell counter
      vec3 X0;                                                  //!< Coordinate of root cell center
      real_t R0;                                                //!< Radius of root cell
      int nspawn;                                               //!< Threshold of NNODE for spawning new threads
      int & numLevels;                                          //!< Maximum tree level
      int level;                                                //!< Current tree level
      int iparent;                                              //!< Index of parent cell
      Nodes2cells(OctreeNode * _octNode, B_iter _B0, C_iter _C, // Constructor
		  C_iter _C0, C_iter _CN, vec3 _X0, real_t _R0,
		  int _nspawn, int & _numLevels, int _level=0, int _iparent=0) :
	octNode(_octNode), B0(_B0), C(_C), C0(_C0), CN(_CN),    // Initialize variables
	X0(_X0), R0(_R0), nspawn(_nspawn), numLevels(_numLevels), level(_level), iparent(_iparent) {}
      //! Get cell index
      uint64_t getKey(vec3 X, vec3 Xmin, real_t diameter) const {
	int iX[3] = {0, 0, 0};                                  // Initialize 3-D index
	for (int d=0; d<3; d++) iX[d] = int((X[d] - Xmin[d]) / diameter);// 3-D index
	uint64_t index = ((1 << 3 * level) - 1) / 7;            // Levelwise offset
	for (int l=0; l<level; l++) {                           // Loop over levels
	  for (int d=0; d<3; d++) index += (iX[d] & 1) << (3 * l + d); // Interleave bits into Morton key
	  for (int d=0; d<3; d++) iX[d] >>= 1;                  //  Bitshift 3-D index
	}                                                       // End loop over levels
	return index;                                           // Return Morton key
      }
      void operator() () const {                                // Overload operator()
	C->IPARENT = iparent;                                   //  Index of parent cell
	C->R       = R0 / (1 << level);                         //  Cell radius
	C->X       = octNode->X;                                //  Cell center
	C->NBODY   = octNode->NBODY;                            //  Number of decendant bodies
	C->IBODY   = octNode->IBODY;                            //  Index of first body in cell
	C->BODY    = B0 + C->IBODY;                             //  Iterator of first body in cell
	C->ICELL   = getKey(C->X, X0-R0, 2*C->R);               //  Get Morton key
	if (octNode->NNODE == 1) {                              //  If node has no children
	  C->ICHILD = 0;                                        //   Set index of first child cell to zero
	  C->NCHILD = 0;                                        //   Number of child cells
	  assert(C->NBODY > 0);                                 //   Check for empty leaf cells
	  numLevels = std::max(numLevels, level);               //   Update maximum level of tree
	} else {                                                //  Else if node has children
	  int nchild = 0;                                       //   Initialize number of child cells
	  int octants[8];                                       //   Map of child index to octants
	  for (int i=0; i<8; i++) {                             //   Loop over octants
	    if (octNode->CHILD[i]) {                            //    If child exists for that octant
	      octants[nchild] = i;                              //     Map octant to child index
	      nchild++;                                         //     Increment child cell counter
	    }                                                   //    End if for child existance
	  }                                                     //   End loop over octants
	  C_iter Ci = CN;                                       //   CN points to the next free memory address
	  C->ICHILD = Ci - C0;                                  //   Set Index of first child cell
	  C->NCHILD = nchild;                                   //   Number of child cells
	  assert(C->NCHILD > 0);                                //   Check for childless non-leaf cells
	  CN += nchild;                                         //   Increment next free memory address
	  mk_task_group;                                        //   Initialize tasks
	  for (int i=0; i<nchild; i++) {                        //   Loop over children
	    int octant = octants[i];                            //    Get octant from child index
	    Nodes2cells nodes2cells(octNode->CHILD[octant],     //    Instantiate recursive functor
				    B0, Ci, C0, CN, X0, R0, nspawn, numLevels, level+1, C-C0);
	    create_taskc_if(octNode->NNODE > nspawn,            //    Spawn task if number of sub-nodes is large
			    nodes2cells);                       //    Recursive call for each child
	    Ci++;                                               //    Increment cell iterator
	    CN += octNode->CHILD[octant]->NNODE - 1;            //    Increment next free memory address
	  }                                                     //   End loop over children
	  wait_tasks;                                           //   Synchronize tasks
	  for (int i=0; i<nchild; i++) {                        //   Loop over children
	    int octant = octants[i];                            //    Get octant from child index
	    delete octNode->CHILD[octant];                      //    Free child pointer to avoid memory leak
	  }                                                     //   End loop over children
	  numLevels = std::max(numLevels, level+1);             //   Update maximum level of tree
	}                                                       //  End if for child existance
      }                                                         // End overload operator()
    };

    //! Transform Xmin & Xmax to X (center) & R (radius)
    Box bounds2box(Bounds bounds) {
      vec3 Xmin = bounds.Xmin;                                  // Set local Xmin
      vec3 Xmax = bounds.Xmax;                                  // Set local Xmax
      Box box;                                                  // Bounding box
      for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2; // Calculate center of domain
      box.R = 0;                                                // Initialize localRadius
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	box.R = std::max(box.X[d] - Xmin[d], box.R);            //  Calculate min distance from center
	box.R = std::max(Xmax[d] - box.X[d], box.R);            //  Calculate max distance from center
      }                                                         // End loop over dimensions
      box.R *= 1.00001;                                         // Add some leeway to radius
      return box;                                               // Return box.X and box.R
    }

    //! Grow tree structure top down
    void growTree(Bodies & bodies, Bodies & buffer, Box box) {
      assert(box.R > 0);                                        // Check for bounds validity
      logger::startTimer("Grow tree");                          // Start timer
      B0 = bodies.begin();                                      // Bodies iterator
      BinaryTreeNode binNode[1];                                // Allocate root node of binary tree
      int maxBinNode = (4 * bodies.size()) / nspawn;            // Get maximum size of binary tree
      binNode->BEGIN = new BinaryTreeNode[maxBinNode];          // Allocate array for binary tree nodes
      binNode->END = binNode->BEGIN + maxBinNode;               // Set end pointer
      buildNodes(N0, bodies, buffer, 0, bodies.size(),          // Build octree nodes
                 binNode, box.X, box.R, ncrit);
      delete[] binNode->BEGIN;                                  // Deallocate binary tree array
      logger::stopTimer("Grow tree");                           // Stop timer
    }

    //! Link tree structure
    Cells linkTree(Box box) {
      logger::startTimer("Link tree");                          // Start timer
      Cells cells;                                              // Initialize cell array
      if (N0 != NULL) {                                         // If the node tree is not empty
	cells.resize(N0->NNODE);                                //  Allocate cells array
	C_iter C0 = cells.begin();                              //  Cell begin iterator
	Nodes2cells nodes2cells(N0, B0, C0, C0, C0+1, box.X, box.R, nspawn, numLevels);// Instantiate recursive functor
	nodes2cells();                                          //  Convert nodes to cells recursively
	delete N0;                                              //  Deallocate nodes
      }                                                         // End if for empty node tree
      logger::stopTimer("Link tree");                           // Stop timer
      return cells;                                             // Return cells array
    }

  public:
    BuildTree(int _ncrit, int _nspawn) : ncrit(_ncrit), nspawn(_nspawn), numLevels(0) {}

    //! Build tree structure top down
    Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
      Box box = bounds2box(bounds);                             // Get box from bounds
      if (bodies.empty()) {                                     // If bodies vector is empty
	N0 = NULL;                                              //  Reinitialize N0 with NULL
      } else {                                                  // If bodies vector is not empty
	if (bodies.size() > buffer.size()) buffer.resize(bodies.size());// Enlarge buffer if necessary
	growTree(bodies, buffer, box);                          //  Grow tree from root
      }                                                         // End if for empty root
      return linkTree(box);                                     // Form parent-child links in tree
    }

    //! Print tree structure statistics
    void printTreeData(Cells & cells) {
      if (logger::verbose && !cells.empty()) {                  // If verbose flag is true
	logger::printTitle("Tree stats");                       //  Print title
	std::cout  << std::setw(logger::stringLength) << std::left//  Set format
		   << "Bodies"     << " : " << cells.front().NBODY << std::endl// Print number of bodies
		   << std::setw(logger::stringLength) << std::left//  Set format
		   << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
		   << std::setw(logger::stringLength) << std::left//  Set format
		   << "Tree depth" << " : " << numLevels << std::endl;//  Print number of levels
      }                                                         // End if for verbose flag
    }
  };
}
#endif
