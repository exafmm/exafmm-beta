#ifndef traversal_h
#define traversal_h
#include "kernel.h"
#include "logger.h"
#include <unordered_map>
#include "thread.h"

#if COUNT
#define count(N) N++
#else
#define count(N)
#endif

class Traversal : public Kernel, public Logger {
 private:
  typedef vec<3,int> ivec3;                                     //!< Vector of 3 integer types
  typedef std::unordered_map<long long,C_iter> CellMap;         //!< Map of cell index to cell iterator
  typedef std::unordered_map<long long,C_iter>::const_iterator M_iter;//!< Iterator for CellMap

  int nspawn;                                                   //!< Threshold of NDBODY for spawning new threads
  int images;                                                   //!< Number of periodic image sublevels
  real_t numP2P;                                                //!< Number of P2P kernel calls
  real_t numM2L;                                                //!< Number of M2L kernel calls
  real_t timeP2P;                                               //!< P2P execution time
  real_t timeM2L;                                               //!< M2L execution time
  C_iter Ci0;                                                   //!< Begin iterator for target cells
  C_iter Cj0;                                                   //!< Begin iterator for source cells

//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj, bool mutual) {
    M2L(Ci,Cj,mutual);                                          // M2L kernel
    count(numM2L);                                              // Increment M2L counter
  }

//! Tree traversal of periodic cells
  void traversePeriodic(real_t cycle) {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Periodic coordinate offset
    Cells pcells(27);                                           // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->CHILD = 0;                                              // Child cells for periodic center cell
    Ci->NCHILD = 26;                                            // Number of child cells for periodic center cell
    C_iter C0 = Cj0;                                            // Placeholder for Cj0
    for (int level=0; level<images-1; level++) {                // Loop over sublevels of tree
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            if (ix != 0 || iy != 0 || iz != 0) {                //     If periodic cell is not at center
              for (int cx=-1; cx<=1; cx++) {                    //      Loop over x periodic direction (child)
                for (int cy=-1; cy<=1; cy++) {                  //       Loop over y periodic direction (child)
                  for (int cz=-1; cz<=1; cz++) {                //        Loop over z periodic direction (child)
                    Xperiodic[0] = (ix * 3 + cx) * cycle;       //         Coordinate offset for x periodic direction
                    Xperiodic[1] = (iy * 3 + cy) * cycle;       //         Coordinate offset for y periodic direction
                    Xperiodic[2] = (iz * 3 + cz) * cycle;       //         Coordinate offset for z periodic direction
                    approximate(Ci0,Ci,false);                  //         Perform M2L kernel
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
#if Cartesian
      for (int i=1; i<NTERM; i++) Ci->M[i] *= Ci->M[0];         //  Normalize multipole expansion coefficients
#endif
      Cj0 = pcells.begin();                                     //  Redefine Cj0 for M2M
      C_iter Cj = Cj0;                                          //  Iterator of periodic neighbor cells
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              Cj->X[0] = Ci->X[0] + ix * cycle;                 //      Set new x coordinate for periodic image
              Cj->X[1] = Ci->X[1] + iy * cycle;                 //      Set new y cooridnate for periodic image
              Cj->X[2] = Ci->X[2] + iz * cycle;                 //      Set new z coordinate for periodic image
              Cj->M    = Ci->M;                                 //      Copy multipoles to new periodic image
              Cj++;                                             //      Increment periodic cell iterator
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Ci->RMAX = 0;                                             //  Initialize Rmax of periodic parent
      Ci->M = 0;                                                //  Reset multipoles of periodic parent
      M2M(Ci,Cj0);                                              //  Evaluate periodic M2M kernels for this sublevel
#if Cartesian
      for (int i=1; i<NTERM; i++) Ci->M[i] /= Ci->M[0];         //  Normalize multipole expansion coefficients
#endif
      cycle *= 3;                                               //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
#if Cartesian
    Ci0->L /= Ci0->M[0];                                        // Normalize local expansion coefficients
#endif
    stopTimer("Traverse periodic");                             // Stop timer
  }

//! Get level from Morton key
  int getLevel(long long key) {
    int level = 0;                                              // Initialize level counter
    while (key > 0) {                                           // While cell index is positive
      level++;                                                  //  Increment level
      key -= 1 << 3*level;                                      //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

//! Get 3-D index from Morton key
  ivec3 getIndex(long long key) {
    ivec3 iX = 0;                                               // Initialize 3-D index
    int d = 0, level = 0;                                       // Initialize dimension and level
    while (key > 0) {                                           // While key is positive
      iX[d] += (key & 1) * (1 << level);                        //  De-interleave bits into 3-D index
      key >>= 1;                                                //  Bitshift key
      d++;                                                      //  Increment dimension
      d = d > 2 ? 0 : d;                                        //  Wrap around dimension
      if (d == 0) level++;                                      //  Increment if dimension was wrapped around
    }                                                           // End while loop for key
    return iX;                                                  // Return 3-D index
  }

//! Get Morton key from 3-D index
  long long getKey(ivec3 iX, int level) {
    long long id = 0;                                           // Levelwise offset
    for (int l=0; l<level; l++) {                               // Loop over levels
      for (int d=0; d<3; d++) id += (iX[d] & 1) << (3 * l + d); //  Interleave bits into Morton key
      for (int d=0; d<3; d++) iX[d] >>= 1;                      //  Bitshift 3-D index
    }                                                           // End loop over levels
    return id;                                                  // Return Morton key
  }

 public:
  Traversal(int nspawn, int images) : nspawn(nspawn), images(images), numP2P(0), numM2L(0) {}

//! Evaluate P2P and M2L using list based traversal
  void dualTreeTraversal(Cells &icells, Cells &jcells, real_t cycle, bool mutual=false) {
    startTimer("Traverse");                                     // Start timer
    CellMap jcellMap;                                           // Map of cell index to cell iterator
    for (C_iter Cj=jcells.begin(); Cj!=jcells.end(); Cj++) {    // Loop over source cells
      jcellMap[Cj->ICELL] = Cj;                                 //  Assign cell iterator to cell index
    }                                                           // End loop over target cells
#pragma omp parallel for
    for (int ci=1; ci<icells.size(); ci++) {                   // Loop over target cells
      C_iter Ci=icells.begin()+ci;                              //  Target cell iterator
      int level = getLevel(Ci->ICELL);                          //  Get level from cell index
      long long targetOffset = ((1 << 3 * level) - 1) / 7;      //  Levelwise offset of target
      long long parentOffset = ((1 << 3 * (level - 1)) - 1) / 7;//  Levelwise offset of parent
      long long targetKey = Ci->ICELL - targetOffset;           //  Morton key from cell index
      long long parentKey = targetKey >> 3;                     //  Parent cell's Morton key
      ivec3 iX = getIndex(targetKey);                           //  3-D index of target cell
      ivec3 pX = getIndex(parentKey);                           //  3-D index of parent cell
      int ixmax = (1 << level) - 1;                             //  Maximum value of iX
      int pxmax = (1 << (level - 1)) - 1;                       //  Maximum value of iP
      ivec3 iXmin, iXmax, pXmin, pXmax;                         //  Declare loop bounds for M2L list
      for (int d=0; d<3; d++) {                                 //  Loop over dimensions
        iXmin[d] = std::max(iX[d]-1, 0);                        //   Minimum of target's neighbor loop
        iXmax[d] = std::min(iX[d]+1, ixmax);                    //   Maximum of target's neighbor loop
        pXmin[d] = std::max(pX[d]-1, 0);                        //   Minimum of parent's neighbor loop
        pXmax[d] = std::min(pX[d]+1, pxmax);                    //   Maximum of parent's neighbor loop
      }                                                         //  End loop over dimensions
      for (pX[0]=pXmin[0]; pX[0]<=pXmax[0]; pX[0]++) {          //  Loop over parent's neighbors x direction
        for (pX[1]=pXmin[1]; pX[1]<=pXmax[1]; pX[1]++) {        //   Loop over parent's neighbors y direction
          for (pX[2]=pXmin[2]; pX[2]<=pXmax[2]; pX[2]++) {      //    Loop over parent's neighbors z direction
            parentKey = getKey(pX, level-1);                    //     Parent's neighbor's Morton key
            M_iter M = jcellMap.find(parentKey+parentOffset);   //     Parent's neighbor's cell map iterator
            if (M != jcellMap.end()) {                          //     If the parent's neighbor exists
              C_iter Cj = M->second;                            //      Parent's neighbor's iterator
              if (Cj->NCHILD == 0) {                            //      If parent's neighbor has no children
                P2P(Ci,Cj,false);                               //       P2P kernel
                count(numP2P);                                  //       Increment P2P counter
              } else {                                          //      Else if parent's neighbor has children
                long long childKey = parentKey << 3;            //       Parent's neighbor's child's Morton key
                for (int i=0; i<8; i++) {                       //       Loop over parent's neighbor's children
                  iX = getIndex(childKey+i);                    //        Parent's neighbor's child's 3-D index
                  if (iX[0]<iXmin[0] || iXmax[0]<iX[0] ||       //        If parent's neighbor's child's 3-D index
                      iX[1]<iXmin[1] || iXmax[1]<iX[1] ||       //        is outside of the target's neighbor reigon
	              iX[2]<iXmin[2] || iXmax[2]<iX[2]) {       //        this cell belongs to the M2L interaction list
                    M = jcellMap.find(childKey+i+targetOffset); //         Cell map iterator
                    if (M != jcellMap.end()) {                  //         If the source cell exists
                      Cj = M->second;                           //          Source cell iterator
                      M2L(Ci,Cj,false);                         //          M2L kernel
                      count(numM2L);                            //          Increment M2L counter
                    }                                           //         End if for source cell existance
                  }                                             //        End if for M2L interaction list
                }                                               //       End loop over parent's neighbor's children
              }                                                 //      End if for parent's neighbor's child's existence
            }                                                   //     End if for parent's neighbor's existence
          }                                                     //    End loop over parent's neighbors z direction
        }                                                       //   End loop over parent's neighbors y direction
      }                                                         //  End loop over parent's neighbors x direction
      if (Ci->NCHILD == 0) {                                    //  If target cell is leaf cell
        for (iX[0]=iXmin[0]; iX[0]<=iXmax[0]; iX[0]++) {        //   Loop over target's neighbors x direction
          for (iX[1]=iXmin[1]; iX[1]<=iXmax[1]; iX[1]++) {      //    Loop over target's neighbors y direction
            for (iX[2]=iXmin[2]; iX[2]<=iXmax[2]; iX[2]++) {    //     Loop over target's neighbors z direction
              long long sourceKey = getKey(iX, level);          //      Target's neighbor's Morton key
              M_iter M = jcellMap.find(sourceKey+targetOffset); //      Cell map iterator
              if (M != jcellMap.end()) {                        //      If the source cell exists
                C_iter Cj = M->second;                          //       Source cell iterator
                P2P(Ci,Cj,false);                               //       P2P kernel
                count(numP2P);                                  //       Increment P2P counter
              }                                                 //      End if for source cell existance
            }                                                   //     End loop over target's neighbors z direction
          }                                                     //    End loop over target's neighbors y direction
        }                                                       //   End loop over target's neighbors x direction
      }                                                         //  End if for target leaf cell 
    }                                                           // End loop over target cells
    stopTimer("Traverse");                                      // Stop timer
    writeTrace();                                               // Write trace to file
  }

//! Direct summation
  void direct(Bodies &ibodies, Bodies &jbodies, real_t cycle) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;            // First cell is target, second cell is source
    Ci->BODY = ibodies.begin();                                 // Iterator of first target body
    Ci->NDBODY = ibodies.size();                                // Number of target bodies
    Cj->BODY = jbodies.begin();                                 // Iterator of first source body
    Cj->NDBODY = jbodies.size();                                // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<images; i++) {                              // Loop over periodic image sublevels
      prange += int(std::pow(3.,i));                            //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
    for (int ix=-prange; ix<=prange; ix++) {                    // Loop over x periodic direction
      for (int iy=-prange; iy<=prange; iy++) {                  //  Loop over y periodic direction
        for (int iz=-prange; iz<=prange; iz++) {                //   Loop over z periodic direction
          Xperiodic[0] = ix * cycle;                            //    Coordinate shift for x periodic direction
          Xperiodic[1] = iy * cycle;                            //    Coordinate shift for y periodic direction
          Xperiodic[2] = iz * cycle;                            //    Coordinate shift for z periodic direction
          P2P(Ci,Cj,false);                                     //    Evaluate P2P kernel
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
  }

//! Normalize bodies after direct summation
  void normalize(Bodies &bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG /= B->SRC;                                         //  Normalize by target charge
    }                                                           // End loop over bodies
  }

//! Time the kernel runtime for auto-tuning
  void timeKernels() {
    Bodies ibodies(1000), jbodies(1000);
    for (B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); Bi++, Bj++) {
      Bi->X = 0;
      Bj->X = 1;
    }
    Cells cells;
    cells.resize(2);
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;
    Ci->X = 0;
    Ci->NDBODY = 10;
    Ci->BODY = ibodies.begin();
    Ci->M = 0;
    Ci->L = 0;
    Cj->X = 1;
    Cj->NDBODY = 1000;
    Cj->BODY = jbodies.begin();
    Cj->M = 0;
    startTimer("P2P kernel test");
    P2P(Ci,Cj,false);
    timeP2P = stopTimer("P2P kernel test") / 10000;
    startTimer("M2L kernel test");
    for (int i=0; i<1000; i++) M2L(Ci,Cj,false);
    timeM2L = stopTimer("M2L kernel test") / 1000;
  }

//! Print traversal statistics
  void printTraversalData() {
#if COUNT
    if (verbose) {                                              // If verbose flag is true
      std::cout << "--- Traversal stats --------------" << std::endl// Print title
	      << std::setw(stringLength) << std::left           //  Set format
	      << "P2P calls"  << " : "                          //  Print title
	      << std::setprecision(0) << std::fixed             //  Set format
              << numP2P << std::endl                            //  Print number of P2P calls
	      << std::setw(stringLength) << std::left           //  Set format
	      << "M2L calls"  << " : "                          //  Print title
	      << std::setprecision(0) << std::fixed             //  Set format
              << numM2L << std::endl;                           //  Print number of M2L calls
    }                                                           // End if for verbose flag
#endif
  }
};
#endif
