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
  using Kernel<equation>::stringLength;                         //!< Max length of event name
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::Ci0;                                  //!< Begin iterator for target cells

private:
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

  int          MAXLEVEL;                                        //!< Maximum level of tree
  B_iter       B0;                                              //!< Iterator of first body
  OctreeNode * N0;                                              //!< Octree root node

private:
//! Get number of binary tree nodes for a given number of bodies
  inline int getNumBinNode(int n) const {
    if (n <= NSPAWN) return 1;                                  // If less then threshold, use only one node
    else return 4 * ((n - 1) / NSPAWN) - 1;                     // Else estimate number of binary tree nodes
  }

//! Get maximum number of binary tree nodes for a given number of bodies
  inline int getMaxBinNode(int n) const {
    return (4 * n) / NSPAWN;                                    // Conservative estimate of number of binary tree nodes
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
    if (end - begin <= NSPAWN) {                                // If number of bodies is less than threshold
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
      countBodies(bodies, begin, mid, X, binNode->LEFT);        //  Recursive call for left branch
      countBodies(bodies, mid, end, X, binNode->RIGHT);         //  Recursive call for right branch
      binNode->NBODY = binNode->LEFT->NBODY + binNode->RIGHT->NBODY;// Sum contribution from both branches
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
      moveBodies(bodies, buffer, begin, mid, binNode->LEFT, octantOffset, X);// Recursive call for left branch
      octantOffset += binNode->LEFT->NBODY;                     //  Increment the octant offset for right branch
      moveBodies(bodies, buffer, mid, end, binNode->RIGHT, octantOffset, X);// Recursive call for right branch
    }                                                           // End if for child existance
  }

//! Build nodes of octree adaptively using a top-down approach based on recursion (uses task based thread parallelism)
  OctreeNode * buildNodes(Bodies& bodies, Bodies& buffer, int begin, int end,
                          BinaryTreeNode * binNode, vec3 X, int level=0, bool direction=false) {
    assert(getMaxBinNode(end - begin) <= binNode->END - binNode->BEGIN);
    if (begin == end) return NULL;                              // If no bodies left, return null pointer
    if (end - begin <= NCRIT) {                                 // If number of bodies is less than threshold
      if (direction)                                            //  If direction of data is from bodies to buffer
        for (int i=begin; i<end; i++) buffer[i] = bodies[i];    //   Copy bodies to buffer
      return makeOctNode(begin, end, X, true);                  //  Create an octree node and return it's pointer
    }                                                           // End if for number of bodies
    OctreeNode * octNode = makeOctNode(begin, end, X, false);   // Create an octree node with child nodes
    countBodies(bodies, begin, end, X, binNode);                // Count bodies in each octant using binary recursion
    ivec8 octantOffset = exclusiveScan(binNode->NBODY, begin);  // Exclusive scan to obtain offset from octant count
    moveBodies(bodies, buffer, begin, end, binNode, octantOffset, X);// Sort bodies according to octant
    BinaryTreeNode * binNodeOffset = binNode->BEGIN;            // Initialize pointer offset for binary tree nodes
    for (int i=0; i<8; i++) {                                   // Loop over children
      int maxBinNode = getMaxBinNode(binNode->NBODY[i]);        //  Get maximum number of binary tree nodes
      assert(binNodeOffset + maxBinNode <= binNode->END);
      vec3 Xchild = X;                                          //  Initialize center position of child node
      real r = R0 / (1 << (level + 1));                         //  Radius of cells for child's level
      for (int d=0; d<3; d++) {                                 //  Loop over dimensions
        Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);         //   Shift center position to that of child node
      }                                                         //  End loop over dimensions
      BinaryTreeNode binNodeChild[1];                           //  Allocate new root for this branch
      binNodeChild->BEGIN = binNodeOffset;                      //  Assign first memory address from offset
      binNodeChild->END = binNodeOffset + maxBinNode;           //  Keep track of last memory address
      octNode->CHILD[i] = buildNodes(buffer, bodies,            //  Recursive call for each child
        octantOffset[i], octantOffset[i] + binNode->NBODY[i],   //  Range of bodies is calcuated from octant offset
        binNodeChild, Xchild, level+1, !direction);             //  Alternate copy direction bodies <-> buffer
      binNodeOffset += maxBinNode;                              //  Increment offset for binNode memory address
    }                                                           // End loop over children
    for (int i=0; i<8; i++) {                                   // Loop over children
      if (octNode->CHILD[i]) octNode->NNODE += octNode->CHILD[i]->NNODE;// If child exists increment child node counter
    }                                                           // End loop over chlidren
    return octNode;                                             // Return octree node
  }

//! Create cell data structure from nodes
  void nodes2cells(OctreeNode * octNode, C_iter C, C_iter CN, int level=0, int iparent=0) {
    C->PARENT = iparent;                                        // Index of parent cell
    C->R      = R0 / (1 << level);                     // Cell radius
    C->X      = octNode->X;                                     // Cell center
    C->NDBODY = octNode->NBODY;                                 // Number of decendant bodies
    C->BODY   = B0 + octNode->BODY;                             // Iterator of first body in cell
    if (octNode->NNODE == 1) {                                  // If node has no children
      C->CHILD  = 0;                                            //  Set index of first child cell to zero
      C->NCHILD = 0;                                            //  Number of child cells
      C->NCBODY = octNode->NBODY;                               //  Number of bodies in cell
      assert(C->NCBODY > 0);
      MAXLEVEL = std::max(MAXLEVEL, level);                     //  Update maximum level of tree
    } else {                                                    // Else if node has children
      C->NCBODY = 0;                                            //  Set number of bodies in cell to zero
      int nchild = 0;                                           //  Initialize number of child cells
      int octants[8];                                           //  Map of child index to octants (for when nchild < 8)
      for (int i=0; i<8; i++) {                                 //  Loop over octants
        if (octNode->CHILD[i]) {                                //   If child exists for that octant
          octants[nchild] = i;                                  //    Map octant to child index
          nchild++;                                             //    Increment child cell counter
        }                                                       //   End if for child existance
      }                                                         //  End loop over octants
      C_iter Ci = CN;                                           //  CN points to the next free memory address
      C->CHILD = Ci - Ci0;                                      //  Set Index of first child cell
      C->NCHILD = nchild;                                       //  Number of child cells
      assert(C->NCHILD > 0);
      CN += nchild;                                             //  Increment next free memory address
      for (int i=0; i<nchild; i++) {                            //  Loop over children
        int octant = octants[i];                                //   Get octant from child index
        nodes2cells(octNode->CHILD[octant], Ci, CN, level+1, C-Ci0);// Recursive call for each child
        Ci++;                                                   //   Increment cell iterator
        CN += octNode->CHILD[octant]->NNODE - 1;                //   Increment next free memory address
      }                                                         //  End loop over children
      for (int i=0; i<nchild; i++) {                            //  Loop over children
        int octant = octants[i];                                //   Get octant from child index
        delete octNode->CHILD[octant];                          //   Free child pointer to avoid memory leak
      }                                                         //  End loop over children
      MAXLEVEL = std::max(MAXLEVEL, level+1);                   //  Update maximum level of tree
    }                                                           // End if for child existance
  }

//! Get Xmin and Xmax of domain
  vec3Pair getBounds(B_iter BiBegin, B_iter BiEnd) {
    assert(BiEnd - BiBegin > 0);
    if (BiEnd - BiBegin < NSPAWN) {                             // If number of elements is small enough
      vec3 Xmin = BiBegin->X, Xmax = BiBegin->X;                //  Initialize Xmin and Xmax with first value
      for (B_iter B=BiBegin; B!=BiEnd; B++) {                   //  Loop over range of bodies
        Xmin = min(B->X, Xmin);                                 //   Update Xmin
        Xmax = max(B->X, Xmax);                                 //   Update Xmax
      }                                                         //  End loop over range of bodies
      return vec3Pair(Xmin, Xmax);                              //  Return Xmin and Xmax as pair
    } else {                                                    // Else if number of elements are large
      B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;           //  Middle iterator
      vec3Pair bounds0, bounds1;                                //  Pair : first Xmin : second Xmax
      bounds0 = getBounds(BiBegin, BiMid);                      //  Recursive call with new task
      bounds1 = getBounds(BiMid, BiEnd);                        //  Recursive call with old task
      bounds0.first = min(bounds0.first, bounds1.first);        //  Minimum of the two Xmins
      bounds0.second = max(bounds0.second, bounds1.second);     //  Maximum of the two Xmaxs
      return bounds0;                                           //  Return Xmin and Xmax
    }                                                           // End if for number fo elements
  }

protected:
//! Grow tree structure top down
  void growTree(Bodies &bodies) {
    Bodies buffer = bodies;                                     // Copy bodies to buffer
    startTimer("Grow tree");                                    // Start timer
    B0 = bodies.begin();                                        // Bodies iterator
    BinaryTreeNode binNode[1];                                  // Allocate root node of binary tree
    int maxBinNode = getMaxBinNode(bodies.size());              // Get maximum size of binary tree
    binNode->BEGIN = new BinaryTreeNode[maxBinNode];            // Allocate array for binary tree nodes
    binNode->END = binNode->BEGIN + maxBinNode;                 // Set end pointer
    N0 = buildNodes(bodies, buffer, 0, bodies.size(), binNode, X0);// Build tree recursively
    delete[] binNode->BEGIN;                                    // Deallocate binary tree array
    stopTimer("Grow tree",printNow);                            // Stop timer
  }

//! Link tree structure
  void linkTree(Cells &cells) {
    startTimer("Link tree");                                    // Start timer
    cells.resize(N0->NNODE);                                    // Allocate cells array
    Ci0 = cells.begin();                                        // Cell begin iterator
    nodes2cells(N0, Ci0, Ci0+1);                                // Convert nodes to cells recursively
    delete N0;                                                  // Deallocate nodes
    stopTimer("Link tree",printNow);                            // Stop timer
  }

  void printTreeData(Cells &cells) {
    std::cout << "--- FMM stats --------------------" << std::endl
    << std::setw(stringLength) << std::left                      // Set format
    << "Bodies"     << " : " << cells.front().NDBODY << std::endl// Print number of bodies
    << std::setw(stringLength) << std::left                      // Set format
    << "Cells"      << " : " << cells.size() << std::endl       // Print number of cells
    << std::setw(stringLength) << std::left                      // Set format
    << "Tree depth" << " : " << MAXLEVEL << std::endl           // Print number of levels
#if COUNT
    << std::setw(stringLength) << std::left                      // Set format
    << "P2P calls"  << " : " << NP2P << std::endl               // Print number of P2P calls
    << std::setw(stringLength) << std::left                      // Set format
    << "M2L calls"  << " : " << NM2L << std::endl               // Print number of M2l calls
#endif
    << "--- FMM stats --------------------" << std::endl;
  }
public:
  TopDown() : MAXLEVEL(0) {}

};

#endif
