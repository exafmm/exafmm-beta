#ifndef build_tree_omp2_h
#define build_tree_omp2_h
#include "logger.h"
#include "thread.h"
#include "types.h"

namespace exafmm {
  class BuildTree {
  private:
    const int ncrit;                                            //!< Number of bodies per leaf cell
    int numLevels;                                              //!< Number of levels

  private:
    //! Get permutation index for reordering bodies
    void reorder(Box box, int level, int * iX, vec3 * Xj,
		 int * permutation, int n, int * iwork, int * nbody) {
      int offset[9];                                            // Offset of bodies per octant
      vec3 X;                                                   // Declare temporary coordinates
      real_t R = box.R / (1 << level);                          // Current cell radius
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	X[d] = box.X[d] - box.R + iX[d] * R * 2 + R;            //  Coordinate of current cell center
      }                                                         // End loop over dimensions
      for (int i=0; i<8; i++) nbody[i] = 0;                     // Initialize nbody array
      for (int i=0; i<n; i++) {                                 // Loop over bodies
	int j = permutation[i];                                 //  Current body index
	int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);// Octant of current body
	nbody[octant]++;                                        //  Increment nbody counter for current octant
      }                                                         // End loop over bodies
      offset[0] = 0;                                            // Initialize offset array
      for (int i=0; i<8; i++) {                                 // Loop over octants
	offset[i+1] = offset[i] + nbody[i];                     //  Calculate offset from nbody
	nbody[i] = 0;                                           //  Initialize nbody again
      }                                                         // End loop over octants
      for (int i=0; i<n; i++) {                                 // Loop over bodies
	int j = permutation[i];                                 //  Current body index
	int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);// Octant of current body
	iwork[offset[octant]+nbody[octant]] = permutation[i];   //  Put permutation index into temporary buffer
	nbody[octant]++;                                        //  Increment nbody counter for current octant
      }                                                         // End loop over bodies
      for (int i=0; i<n; i++) {                                 // Loop over bodies
	permutation[i] = iwork[i];                              //  Copy back permutation array
      }                                                         // End loop over bodies
    }

    //! Get Morton key
    uint64_t getKey(ivec3 iX, int level) {
      uint64_t index = ((1 << 3 * level) - 1) / 7;              // Level offset
      for (int l=0; l<level; l++) {                             // Loop over levels
	for (int d=0; d<3; d++) {                               //  Loop over dimensions
	  index += (iX[d] & 1) << (3 * l + d);                  //   Interleave bits for each dimension
	  iX[d] >>= 1;                                          //   Shift bits for each dimension
	}                                                       //  End loop over dimensions
      }                                                         // End loop over levels
      return index;                                             // Return Morton key
    }

    //! Get level from Morton key
    int getLevel(uint64_t key) {
      int level = -1;                                           // Initialize level
      while( int(key) >= 0 ) {                                  // Loop while key has level offsets to subtract
	level++;                                                //  Increment level
	key -= 1 << 3*level;                                    //  Subtract level offset from key
      }                                                         // End while loop for level offsets to subtract
      return level;                                             // Return level
    }

    //! Transform Xmin & Xmax to X (center) & R (radius)
    Box bounds2box(Bounds bounds) {
      vec3 Xmin = bounds.Xmin;                                  // Set local Xmin
      vec3 Xmax = bounds.Xmax;                                  // Set local Xmax
      Box box;                                                  // Bounding box
      for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2;// Calculate center of domain 
      box.R = 0;                                                // Initialize cell radius
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	box.R = std::max(box.X[d] - Xmin[d], box.R);            //  Calculate min distance from center
	box.R = std::max(Xmax[d] - box.X[d], box.R);            //  Calculate max distance from center
      }                                                         // End loop over dimensions
      box.R *= 1.00001;                                         // Add some leeway to radius 
      return box;                                               // Return box.X and box.R 
    }

    //! Grow tree as link between node structures
    void growTree(Bodies & bodies, int (* nodes)[10], int & numCells,
		  int * permutation, int & numLevels, Box box) {
      logger::startTimer("Grow tree");                          // Start timer
      const int maxLevel = 30;                                  // Maximum levels in tree
      const int numBodies = bodies.size();                      // Number of bodies
      int nbody8[8];                                            // Number of bodies per octant
      int * iwork = new int [numBodies];                        // Allocate temporary work array of integers
      int * levelOffset = new int [maxLevel];                   // Allocate level offset array
      vec3 * Xj = new vec3 [numBodies];                         // Allocate temporary coordinate array
      nodes[0][0] = 0;                                          // Initialize level
      nodes[0][1] = 0;                                          // Initialize ix
      nodes[0][2] = 0;                                          // Initialize iy
      nodes[0][3] = 0;                                          // Initialize iz
      nodes[0][4] = 0;                                          // Initialize iparent
      nodes[0][5] = 0;                                          // Initialize ichild
      nodes[0][6] = 0;                                          // Initialize nchild
      nodes[0][7] = 0;                                          // Initialize ibody
      nodes[0][8] = numBodies;                                  // Initialize nbody
      levelOffset[0] = 0;                                       // Offset for level 0
      levelOffset[1] = 1;                                       // Offset for level 1
      for (int i=0; i<numBodies; i++) {                         // Loop over bodies
	permutation[i] = i;                                     //  Copy index
	Xj[i] = bodies[i].X;                                    //  Copy coordinates
      }                                                         // End loop over bodies
      numCells = 1;                                             // Initialize number of cells
      numLevels = 0;                                            // Initialize number of levels
      for (int level=0; level<maxLevel; level++) {              // Loop over levels
	for (int iparent=levelOffset[level]; iparent<levelOffset[level+1]; iparent++) {// Loop over cells in level
	  int nbody = nodes[iparent][8];                        //   Number of bodies in current cell
	  if (nbody > ncrit) {                                  //   If number of bodies is larger than threshold
	    int ibody = nodes[iparent][7];                      //    Index of first body in cell
	    reorder(box, level, &nodes[iparent][1], Xj, &permutation[ibody], nbody, iwork, nbody8);// Sort bodies
	    int nchild = 0;                                     //    Initialize number of child cells
	    int offset = ibody;                                 //    Initialize offset
	    nodes[iparent][5] = numCells;                       //    Store cell counter as ichild
	    for (int i=0; i<8; i++) {                           //    Loop over octants
	      nodes[numCells][0] = level + 1;                   //     Store level
	      nodes[numCells][1] = nodes[iparent][1] * 2 + i % 2;//    Store ix
	      nodes[numCells][2] = nodes[iparent][2] * 2 + (i / 2) % 2;// Store iy
	      nodes[numCells][3] = nodes[iparent][3] * 2 + i / 4;//    Store iz
	      nodes[numCells][4] = iparent;                     //     Store iparent
	      nodes[numCells][5] = 0;                           //     Initialize ichild
	      nodes[numCells][6] = 0;                           //     Initialize nchild
	      nodes[numCells][7] = offset;                      //     Store ibody
	      nodes[numCells][8] = nbody8[i];                   //     Store nbody
	      nchild++;                                         //     Increment number of child cells
	      offset += nbody8[i];                              //     Increment octant offset
	      numCells++;                                       //     Increment number of cells
	      numLevels=level+1;                                //     Update number of levels
	    }                                                   //    End loop over octants
	    nodes[iparent][6] = nchild;                         //    Store nchild
	  }                                                     //   End if for number of bodies threshold
	}                                                       //  End loop over cells in level
	levelOffset[level+2] = numCells;                        //  Update level offset
	if (levelOffset[level+1] == levelOffset[level+2]) break;//  If no cells were added then exit loop
      }                                                         // End loop over levels
      delete[] Xj;                                              // Deallocate temporary coordinate array
      delete[] levelOffset;                                     // Deallocate level offset array 
      delete[] iwork;                                           // Deallocate temporary work array of integers
      logger::stopTimer("Grow tree");                           // Stop timer
    }

    //! Convert nodes to cells
    Cells linkTree(Bodies & bodies, Bodies & buffer, int (* nodes)[10], int numCells,
		   int * permutation, Box box) {
      logger::startTimer("Link tree");                          // Start timer
      int numBodies = bodies.size();                            // Number of bodies
      Cells cells(numCells);                                    // Instantiate cells vector
      C_iter C = cells.begin();                                 // Iterator of first cell
      ivec3 iX;                                                 // 3-D index
      for (int i=0; i<numCells; i++,C++) {                      // Loop over cells
	int level = nodes[i][0];                                //  Current level
	iX[0]      = nodes[i][1];                               //  3-D index x component
	iX[1]      = nodes[i][2];                               //  3-D index y component
	iX[2]      = nodes[i][3];                               //  3-D index z component
	C->ICELL   = getKey(iX, level);                         //  Copy Morton key as icell
	C->IPARENT = nodes[i][4];                               //  Copy iparent
	C->ICHILD  = nodes[i][5];                               //  Copy ichild
	C->NCHILD  = nodes[i][6];                               //  Copy nchild
	C->IBODY   = nodes[i][7];                               //  Copy ibody
	C->NBODY   = nodes[i][8];                               //  Copy nbody
	real_t R = box.R / (1 << level);                        //  Cell radius
	C->R = R;                                               //  Store cell radius
	for (int d=0; d<3; d++) {                               //  Loop over dimensions
	  C->X[d] = box.X[d] - box.R + iX[d] * R * 2 + R;       //   Center of cell
	}                                                       //  End loop over dimensions
      }                                                         // End loop over cells
      buffer.resize(numBodies);                                 // Resize buffer
      for (int i=0; i<numBodies; i++) {                         // Loop over bodies
	buffer[i] = bodies[permutation[i]];                     //  Copy permuted bodies to buffer
      }                                                         // End loop over bodies
      bodies = buffer;                                          // Copy back to bodies
      B_iter B = bodies.begin();                                // Iterator of first body
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	C->BODY = B + C->IBODY;                                 //  Store iterator of first body in cell
      }                                                         // End loop over cells
      logger::stopTimer("Link tree");                           // Stop timer
      return cells;                                             // Return cells
    }

  public:
    BuildTree(int _ncrit, int ) : ncrit(_ncrit) {}              // Constructor

    //! Build tree structure
    Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
      int numCells;                                             // Number of cells
      int numBodies = bodies.size();                            // Number of bodies
      int (* nodes)[10] = new int [numBodies][10]();            // Allocate nodes array
      int * permutation = new int [numBodies];                  // Allocate permutation array
      Box box = bounds2box(bounds);                             // Bounding box
      growTree(bodies, nodes, numCells, permutation, numLevels, box);// Grow tree as link between node structures
      Cells cells = linkTree(bodies, buffer, nodes, numCells, permutation, box);// Convert nodes to cells
      delete[] permutation;                                     // Deallocate permutation array
      delete[] nodes;                                           // Deallocate nodes array
      return cells;                                             // Return cells
    }

    void printTreeData(Cells & cells) {
      if (logger::verbose && !cells.empty()) {                  // If verbose flag is true
	logger::printTitle("Tree stats");                       //  Print title
	std::cout  << std::setw(logger::stringLength) << std::left// Set format
		   << "Bodies"     << " : " << cells.front().NBODY << std::endl// Print number of bodies
		   << std::setw(logger::stringLength) << std::left// Set format
		   << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
		   << std::setw(logger::stringLength) << std::left// Set format
		   << "Tree depth" << " : " << numLevels << std::endl;//  Print number of levels
      }                                                         // End if for verbose flag
    }
  };
}
#endif
