#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

namespace exafmm {
  class BuildTreeFromCluster {
  public:
    typedef std::vector<int> ints;                              //!< Vector of integer types
    typedef std::vector<vec3> vec3s;                            //!< Vector of vec3 types
    ints iwrap;                                                 //!< Bit flag for wrap around periodic boundary

    BuildTreeFromCluster() {}                                   // Constructor

    //! Get number of cells (clusters) from body.ICELL
    int getNumCells(Bodies & bodies) {
      int icell = -1;                                           // Initialize cell index
      int numCells = 0;                                         // Initialize number of cells
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	int index = B->ICELL;                                   //  Current cell index
	if (index != icell) {                                   //  If cell index is different from the previous
	  numCells++;                                           //   Increment number of cells
	  icell = index;                                        //   Update current cell index
	}                                                       //  End if for different cell index
      }                                                         // End loop over bodies
      return numCells;                                          // Return number of cells
    }

    //! Get Xmin
    vec3s getXmin(Bodies & bodies, int numCells) {
      vec3s Xmin(numCells);                                     // Vector of Xmins
      int i = 0;                                                // Initialize cell counter
      int icell = bodies.begin()->ICELL;                        // First cell index
      Xmin[0] = bodies.begin()->X;                              // Initialize Xmin
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	int index = B->ICELL;                                   //  Current cell index
	if (index != icell) {                                   //  If cell index is different from the previous
	  i++;                                                  //   Increment cell counter
	  icell = index;                                        //   Update current cell index 
	  Xmin[i] = B->X;                                       //   Initialize Xmin
	}                                                       //  End if for different cell index
	Xmin[i] = min(Xmin[i], B->X);                           //  Update Xmin
      }                                                         // End loop over bodies
      return Xmin;                                              // Return Xmin
    }

    //! Set center coordinate of cluster
    Bodies setClusterCenter(Bodies & bodies, real_t cycle) {
      iwrap.resize(bodies.size());                              // Resize iwrap
      int numCells = getNumCells(bodies);                       // Get number of cells
      vec3s Xmin = getXmin(bodies, numCells);                   // Get Xmin
      Bodies cluster(numCells);                                 // Declare clusters as body vector type
      int icell = bodies.begin()->ICELL;                        // First cell index
      int numBodies = 0;                                        // Initialize number of bodies
      B_iter C=cluster.begin();                                 // Iterator of first clusters
      int b = 0, c = 0;                                         // Initialize body and cell counters
      C->X = 0;                                                 // Initialize center coordinate of cluster
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {// Loop over bodies
	int index = B->ICELL;                                   //  Current cell index
	if (index != icell) {                                   //  If cell index is different from the previous
	  assert(index > icell);                                //   Make sure cell index is accending
	  C->X /= numBodies;                                    //   Divide sum by number to get average
	  C->IBODY = B - bodies.begin() - numBodies;            //   Index of first body in current cell
	  C->ICELL = icell;                                     //   Store cell index
	  C++;                                                  //   Increment cell iterator
	  c++;                                                  //   Increment cell counter
	  C->X = 0;                                             //   Initialize center coordinate of cluster
	  numBodies = 0;                                        //   Initialize number of bodies
	  icell = index;                                        //   Update current cell index
	}                                                       //  End if for different cell index
	iwrap[b] = 0;                                           //  Initialize iwrap
	for (int d=0; d<3; d++) {                               //  Loop over dimensions
	  int flag = B->X[d] - Xmin[c][d] > cycle / 2;          //   ?
	  B->X[d] -= cycle * flag;                              //   ?
	  iwrap[b] |= flag << d;                                //   ?
	}                                                       //  End loop over dimensions
	C->X += B->X;                                           //  Accumulate body coordinates
	numBodies++;                                            //  Increment number of bodies
      }                                                         // End loop over bodies
      C->X /= numBodies;                                        // Divide sum by number to get average 
      C->IBODY = bodies.size() - numBodies;                     // Index of first body in last cell
      C->ICELL = icell;                                         // Store cell index
      return cluster;                                           // Return cluster
    }

    //! Upward pass (post order traversal) to determine NBODY, R, X
    void upwardPass(C_iter C, C_iter C0) {
      for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	upwardPass(CC, C0);                                     //  Recursive call to child cells
	C->NBODY += CC->NBODY;                                  //  Acculumate number of bodies
	C->R = std::max(C->R, 2 * CC->R);                       //  Update cell radius
	C->X += CC->X;                                          //  Accumulate cell coordinate
      }                                                         // End loop over child cells
      if (C->NCHILD != 0) C->X /= C->NCHILD;                    // Divide sum by number to get average
    }

    //! Turn clusters into cells
    void attachClusterBodies(Bodies & bodies, Cells & cells, real_t cycle) {
      B_iter B0 = bodies.begin();                               // Iterator of first body
      int numLeafs = 0;                                         // Initialize number of leafs
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	if (C->NCHILD == 0) numLeafs++;                         //  If cell has no child, increment number of leafs
      }                                                         // End loop over cells
      int dimLeafs = cbrt(numLeafs);                            // Number of leafs per dimension
      assert(dimLeafs*dimLeafs*dimLeafs == numLeafs);           // Check if numLeafs is a cube of an integer
      int numLevels = log2(dimLeafs) + 1;                       // Maximum number of levels
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	if (C->NCHILD == 0) {                                   //  If leaf cell
	  uint64_t icell = C->BODY->ICELL;                      //   Current cell index
	  B_iter B = B0 + C->BODY->IBODY;                       //   Iterator of first body in cell
	  int ix = icell / dimLeafs / dimLeafs;                 //   Index in x dimension
	  int iy = icell / dimLeafs % dimLeafs;                 //   Index in y dimension
	  int iz = icell % dimLeafs;                            //   Index in z dimension
	  int key = 0;                                          //   Initialize key
	  for (int l=0; l<numLevels; l++) {                     //   Loop over levels
	    key += (ix & 1) << (3 * l);                         //    Interleave x bits
	    key += (iy & 1) << (3 * l + 1);                     //    Interleave y bits
	    key += (iz & 1) << (3 * l + 2);                     //    Interleave z bits
	    ix >>= 1;                                           //    Shift x bits
	    iy >>= 1;                                           //    Shift y bits
	    iz >>= 1;                                           //    Shift z bits
	  }                                                     //   End loop over levels
	  C->X = C->BODY->X;                                    //   Store coordinates
	  C->BODY = B;                                          //   Store body iterator
	  C->IBODY = B - B0;                                    //   Store body index
	  C->NBODY = 0;                                         //   Store number of bodies
	  C->ICELL = key;                                       //   Store Morton key as cell index
	  C->R = cycle / dimLeafs / 2;                          //   Store cell radius
	  while (B->ICELL == icell) {                           //   Loop while icell is same (same cell)
	    C->NBODY++;                                         //    Increment number of bodies
	    if (B==bodies.end()-1) break;                       //    If end of bodies then exit loop
	    B++;                                                //    Increment body iterator
	  }                                                     //   End while loop for same icell
	} else {                                                //  If cell is not leaf
	  C->NBODY = 0;                                         //   Initialize number of bodies
	  C->ICELL = 0;                                         //   Initialize cell index
	  C->R = 0;                                             //   Initialize cell radius
	  C->X = 0;                                             //   Initialize cell center
	}                                                       //  End if for leaf cell
      }                                                         // End loop over cells
      upwardPass(cells.begin(), cells.begin());                 // Construct the non-leaf cells
    }

    //! Shift back bodies that were wrapped
    void shiftBackBodies(Bodies & bodies, real_t cycle) {
      int b = 0;                                                // Initialize body counter
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {// Loop over bodies
	unwrap(B->X,cycle,iwrap[b]);                            //  Unwrap body coordinate
      }                                                         // End loop over bodies
    }
  };
}
#endif
