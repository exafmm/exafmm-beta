#ifndef build_tree_tbb_h
#define build_tree_tbb_h
#include "logger.h"
#include "namespace.h"
#include "thread.h"
#include "types.h"

namespace EXAFMM_NAMESPACE {
  class BuildTree {
    typedef vec<8,int> ivec8;                                   //!< Vector of 8 integer types

  private:
    //! Octree is used for building the FMM tree structure as "nodes", then transformed to "cells" data structure
    struct OctreeNode {
      int          IBODY;                                       //!< Index offset for first body in node
      int          NBODY;                                       //!< Number of descendant bodies
      int          NNODE;                                       //!< Number of descendant nodes
      OctreeNode * CHILD[8];                                    //!< Pointer to child node
      vec3         X;                                           //!< Coordinate at center
    };

    const int    ncrit;                                         //!< Number of bodies per leaf cell
    int          numLevels;                                     //!< Number of levels in tree
    B_iter       B0;                                            //!< Iterator of first body
    OctreeNode * N0;                                            //!< Pointer to octree root node

  private:
    //! Counting bodies in each octant
    void countBodies(Bodies & bodies, int begin, int end, vec3 X, ivec8 & NBODY) {
      for (int i=0; i<8; i++) NBODY[i] = 0;                     // Initialize number of bodies in octant
      for (int i=begin; i<end; i++) {                           // Loop over bodies in node
        vec3 x = bodies[i].X;                                   //  Coordinates of body
        if (bodies[i].ICELL < 0)                                //  If using residual index
          x = bodies[i+bodies[i].ICELL].X;                      //   Use coordinates of first body in residual group
        int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to
        NBODY[octant]++;                                        //  Increment body count in octant
      }                                                         // End loop over bodies in node
    }

    //! Sorting bodies according to octant (Morton order)
    void moveBodies(Bodies & bodies, Bodies & buffer, int begin, int end,
                           ivec8 octantOffset, vec3 X) {
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
    OctreeNode * makeOctNode(OctreeNode *& octNode, int begin, int end, vec3 X, bool nochild) {
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
    inline ivec8 exclusiveScan(ivec8 input, int offset) {
      ivec8 output;                                             // Output vector
      for (int i=0; i<8; i++) {                                 // Loop over elements
        output[i] = offset;                                     //  Set value
        offset += input[i];                                     //  Increment offset
      }                                                         // End loop over elements
      return output;                                            // Return output vector
    }

    //! Recursive functor for building nodes of an octree adaptively using a top-down approach
    void buildNodes(OctreeNode *& octNode, Bodies & bodies, Bodies & buffer,
                           int begin, int end, vec3 X, real_t R0,
                           int level=0, bool direction=false) {
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
      ivec8 NBODY;                                              // Number of bodies in node
      countBodies(bodies, begin, end, X, NBODY);                // Count bodies in each octant
      ivec8 octantOffset = exclusiveScan(NBODY, begin);         // Exclusive scan to obtain offset from octant count
      moveBodies(bodies, buffer, begin, end, octantOffset, X);  // Sort bodies according to octant
      for (int i=0; i<8; i++) {                                 // Loop over children
        vec3 Xchild = X;                                        //  Initialize center coordinates of child node
        real_t r = R0 / (1 << (level + 1));                     //  Radius of cells for child's level
        for (int d=0; d<3; d++) {                               //  Loop over dimensions
          Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);       //   Shift center coordinates to that of child node
        }                                                       //  End loop over dimensions
        buildNodes(octNode->CHILD[i], buffer, bodies,           //  Recursive call for children
                   octantOffset[i], octantOffset[i] + NBODY[i],
                   Xchild, R0, level+1, !direction);
      }                                                         // End loop over children
      for (int i=0; i<8; i++) {                                 // Loop over children
        if (octNode->CHILD[i]) octNode->NNODE += octNode->CHILD[i]->NNODE;// If child exists increment child node count
      }                                                         // End loop over chlidren
    }

    //! Get Morton key
    uint64_t getKey(vec3 X, vec3 Xmin, real_t diameter, int level) const {
      int iX[3] = {0, 0, 0};                                    // Initialize 3-D index
      for (int d=0; d<3; d++) iX[d] = int((X[d] - Xmin[d]) / diameter);// 3-D index
      uint64_t index = ((1 << 3 * level) - 1) / 7;              // Levelwise offset
      for (int l=0; l<level; l++) {                             // Loop over levels
        for (int d=0; d<3; d++) index += (iX[d] & 1) << (3 * l + d); // Interleave bits into Morton key
        for (int d=0; d<3; d++) iX[d] >>= 1;                    //  Bitshift 3-D index
      }                                                         // End loop over levels
      return index;                                             // Return Morton key
    }

    //! Creating cell data structure from nodes
    void nodes2cells(OctreeNode * octNode, C_iter C,
                     C_iter C0, C_iter CN, vec3 X0, real_t R0,
                     int & maxLevel, int level=0, int iparent=0) {
      C->IPARENT = iparent;                                     //  Index of parent cell
      C->R       = R0 / (1 << level);                           //  Cell radius
      C->X       = octNode->X;                                  //  Cell center
      C->NBODY   = octNode->NBODY;                              //  Number of decendant bodies
      C->IBODY   = octNode->IBODY;                              //  Index of first body in cell
      C->BODY    = B0 + C->IBODY;                               //  Iterator of first body in cell
      C->ICELL   = getKey(C->X, X0-R0, 2*C->R, level);          //  Get Morton key
      if (octNode->NNODE == 1) {                                //  If node has no children
        C->ICHILD = 0;                                          //   Set index of first child cell to zero
        C->NCHILD = 0;                                          //   Number of child cells
        assert(C->NBODY > 0);                                   //   Check for empty leaf cells
        maxLevel = std::max(maxLevel, level);                   //   Update maximum level of tree
      } else {                                                  //  Else if node has children
        int nchild = 0;                                         //   Initialize number of child cells
        int octants[8];                                         //   Map of child index to octants
        for (int i=0; i<8; i++) {                               //   Loop over octants
          if (octNode->CHILD[i]) {                              //    If child exists for that octant
            octants[nchild] = i;                                //     Map octant to child index
            nchild++;                                           //     Increment child cell counter
          }                                                     //    End if for child existance
        }                                                       //   End loop over octants
        C_iter Ci = CN;                                         //   CN points to the next free memory address
        C->ICHILD = Ci - C0;                                    //   Set Index of first child cell
        C->NCHILD = nchild;                                     //   Number of child cells
        assert(C->NCHILD > 0);                                  //   Check for childless non-leaf cells
        CN += nchild;                                           //   Increment next free memory address
        for (int i=0; i<nchild; i++) {                          //   Loop over children
          int octant = octants[i];                              //    Get octant from child index
          nodes2cells(octNode->CHILD[octant], Ci, C0, CN,       //    Recursive call for child cells
                      X0, R0, numLevels, level+1, C-C0);
          Ci++;                                                 //    Increment cell iterator
          CN += octNode->CHILD[octant]->NNODE - 1;              //    Increment next free memory address
        }                                                       //   End loop over children
        for (int i=0; i<nchild; i++) {                          //   Loop over children
          int octant = octants[i];                              //    Get octant from child index
          delete octNode->CHILD[octant];                        //    Free child pointer to avoid memory leak
        }                                                       //   End loop over children
        maxLevel = std::max(maxLevel, level+1);                 //   Update maximum level of tree
      }                                                         //  End if for child existance
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

  public:
    BuildTree(int _ncrit) : ncrit(_ncrit), numLevels(0) {}

    //! Build tree structure top down
    Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
      logger::startTimer("Grow tree");                          // Start timer
      Box box = bounds2box(bounds);                             // Get box from bounds
      if (bodies.empty()) {                                     // If bodies vector is empty
	N0 = NULL;                                              //  Reinitialize N0 with NULL
      } else {                                                  // If bodies vector is not empty
	if (bodies.size() > buffer.size()) buffer.resize(bodies.size());// Enlarge buffer if necessary
        assert(box.R > 0);                                      // Check for bounds validity
        B0 = bodies.begin();                                    // Bodies iterator
        buildNodes(N0, bodies, buffer, 0, bodies.size(),        // Build octree nodes
                   box.X, box.R);
      }                                                         // End if for empty root
      logger::stopTimer("Grow tree");                           // Stop timer
      logger::startTimer("Link tree");                          // Start timer
      Cells cells;                                              // Initialize cell array
      if (N0 != NULL) {                                         // If the node tree is not empty
	cells.resize(N0->NNODE);                                //  Allocate cells array
	C_iter C0 = cells.begin();                              //  Cell begin iterator
	nodes2cells(N0, C0, C0, C0+1, box.X, box.R, numLevels); // Instantiate recursive functor
	delete N0;                                              //  Deallocate nodes
      }                                                         // End if for empty node tree
      logger::stopTimer("Link tree");                           // Stop timer
      return cells;                                             // Return cells array
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
