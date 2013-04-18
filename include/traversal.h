#ifndef traversal_h
#define traversal_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

#if COUNT
#define count(N) N++
#else
#define count(N)
#endif

class Traversal : public Kernel, public Logger {
 private:
  real_t timeP2P;                                               //!< P2P execution time
  real_t timeM2L;                                               //!< M2L execution time
  real_t NP2P;                                                  //!< Number of P2P kernel calls
  real_t NM2L;                                                  //!< Number of M2L kernel calls
  C_iter Ci0;                                                   //!< Begin iterator for target cells
  C_iter Cj0;                                                   //!< Begin iterator for source cells

 public:
  int NCRIT;                                                    //!< Number of bodies per leaf cell
  int NSPAWN;                                                   //!< Threshold of NDBODY for spawning new threads
  int IMAGES;                                                   //!< Number of periodic image sublevels
  float THETA;                                                  //!< Multipole acceptance criteria 

  real_t localRadius;                                           //!< Radius of local root cell
  vec3   localCenter;                                           //!< Center of local root cell
  fvec3  localXmin;                                             //!< Local Xmin for a given rank
  fvec3  localXmax;                                             //!< Local Xmax for a given rank

  real_t periodicCycle;                                         //!< Diameter of global root cell

 private:
//! Calculate Bmax
  real_t getBmax(vec3 const &X, C_iter C) const {
    real_t rad = C->R;                                          // Radius of cell
    real_t dx = rad + std::abs(X[0]-C->X[0]);                   // Add x distance from center of mass
    real_t dy = rad + std::abs(X[1]-C->X[1]);                   // Add y distance from center of mass
    real_t dz = rad + std::abs(X[2]-C->X[2]);                   // Add z distance from center of mass
    return std::sqrt(dx * dx + dy * dy + dz * dz);              // Return scalar distance
  }

//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj, bool mutual) {
#if AUTO
    if (timeP2P*Ci->NDBODY*Cj->NDBODY > timeM2L) {              // If M2L is faster
      M2L(Ci, Cj, mutual);                                      //  M2L kernel
      count(NM2L);                                              //  Increment M2L counter
    } else {                                                    // Else if P2P is faster
      P2P(Ci, Cj, mutual);                                      //  P2P kernel
      count(NP2P);                                              //  Increment P2P counter
    }                                                           // End if for fastest kernel
#else
    M2L(Ci,Cj,mutual);                                          // M2L kernel
    count(NM2L);                                                // Increment M2L counter
#endif
  }

//! Dual tree traversal for a range of Ci and Cj
  void traverse(C_iter CiBegin, C_iter CiEnd, C_iter CjBegin, C_iter CjEnd, bool mutual) {
    if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {         // If only one cell in range
      if (CiBegin == CjBegin) {                                 //  If Ci == Cj
        assert(CiEnd == CjEnd);
        traverse(CiBegin, CjBegin, mutual);                     //   Call traverse for single pair
      } else {                                                  //  If Ci != Cj
        for (C_iter Ci=CiBegin; Ci!=CiEnd; Ci++) {              //   Loop over all Ci cells
          for (C_iter Cj=CjBegin; Cj!=CjEnd; Cj++) {            //    Loop over all Cj cells
            traverse(Ci, Cj, mutual);                           //     Call traverse for single pair
          }                                                     //    End loop over all Cj cells
        }                                                       //   End loop over all Ci cells
      }                                                         //  End if for Ci == Cj
    } else {                                                    // If many cells are in the range
      C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;           //  Split range of Ci cells in half
      C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;           //  Split range of Cj cells in half
      spawn_tasks {                                             //  Initialize task group
	spawn_task0(traverse(CiBegin, CiMid, CjBegin, CjMid, mutual));// Spawn Ci:former Cj:former
	traverse(CiMid, CiEnd, CjMid, CjEnd, mutual);             //  No spawn Ci:latter Cj:latter
	sync_tasks;                                           //  Synchronize task group
	spawn_task0(traverse(CiBegin, CiMid, CjMid, CjEnd, mutual));// Spawn Ci:former Cj:latter
	if (!mutual || CiBegin != CjBegin) {                      //  Exclude mutual & self interaction
	  traverse(CiMid, CiEnd, CjBegin, CjMid, mutual);         //   No spawn Ci:latter Cj:former
	} else {                                                  //  If mutual or self interaction
	  assert(CiEnd == CjEnd);                                 //   Check if mutual & self interaction
	}                                                         //  End if for mutual & self interaction
	sync_tasks;                                           //  Synchronize task group
      }
    }                                                           // End if for many cells in range
  }

//! Split cell and call traverse() recursively for child
  void splitCell(C_iter Ci, C_iter Cj, bool mutual) {
    if (Cj->NCHILD == 0) {                                      // If Cj is leaf
      assert(Ci->NCHILD > 0);                                   //  Make sure Ci is not leaf
      for (C_iter ci=Ci0+Ci->CHILD; ci!=Ci0+Ci->CHILD+Ci->NCHILD; ci++ ) {// Loop over Ci's children
        traverse(ci, Cj, mutual);                               //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else if (Ci->NCHILD == 0) {                               // Else if Ci is leaf
      assert(Cj->NCHILD > 0);                                   //  Make sure Cj is not leaf
      for (C_iter cj=Cj0+Cj->CHILD; cj!=Cj0+Cj->CHILD+Cj->NCHILD; cj++ ) {// Loop over Cj's children
        traverse(Ci, cj, mutual);                               //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    } else if (Ci->NDBODY + Cj->NDBODY >= NSPAWN || (mutual && Ci == Cj)) {// Else if cells are still large
      traverse(Ci0+Ci->CHILD, Ci0+Ci->CHILD+Ci->NCHILD,         //  Traverse for range of cell pairs
               Cj0+Cj->CHILD, Cj0+Cj->CHILD+Cj->NCHILD, mutual);
    } else if (Ci->RCRIT >= Cj->RCRIT) {                        // Else if Ci is larger than Cj
      for (C_iter ci=Ci0+Ci->CHILD; ci!=Ci0+Ci->CHILD+Ci->NCHILD; ci++ ) {// Loop over Ci's children
        traverse(ci, Cj, mutual);                               //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Cj is larger than Ci
      for (C_iter cj=Cj0+Cj->CHILD; cj!=Cj0+Cj->CHILD+Cj->NCHILD; cj++ ) {// Loop over Cj's children
        traverse(Ci, cj, mutual);                               //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

//! Dual tree traversal for a single pair of cells
  void traverse(C_iter Ci, C_iter Cj, bool mutual) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real_t R2 = norm(dX);                                       // Scalar distance squared
#if DUAL
    {                                                           // Dummy bracket
#else
    if (Ci->RCRIT != Cj->RCRIT) {                               // If cell is not at the same level
      splitCell(Ci, Cj, mutual);                                //  Split cell and call function recursively for child
    } else {                                                    // If we don't care if cell is not at the same level
#endif
      if (R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {   //  If distance is far enough
        approximate(Ci, Cj, mutual);                            //   Use approximate kernels
      } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {          //  Else if both cells are bodies
        if (Cj->NCBODY == 0) {                                  //   If the bodies weren't sent from remote node
          approximate(Ci, Cj, mutual);                          //    Use approximate kernels
        } else {                                                //   Else if the bodies were sent
          if (Ci == Cj) {                                       //    If source and target are same
            P2P(Ci);                                            //     P2P kernel for single cell
          } else {                                              //    Else if source and target are different
            P2P(Ci, Cj, mutual);                                //     P2P kernel for pair of cells
          }                                                     //    End if for same source and target
          count(NP2P);                                          //    Increment P2P counter
        }                                                       //   End if for bodies
      } else {                                                  //  Else if cells are close but not bodies
        splitCell(Ci, Cj, mutual);                              //   Split cell and call function recursively for child
      }                                                         //  End if for multipole acceptance
    }                                                           // End if for same level cells
  }


//! Tree traversal of periodic cells
  void traversePeriodic(real_t Length) {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Periodic coordinate offset
    Cells pcells(27);                                           // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->CHILD = 0;                                              // Child cells for periodic center cell
    Ci->NCHILD = 26;                                            // Number of child cells for periodic center cell
    C_iter C0 = Cj0;                                            // Placeholder for Cj0
    for (int level=0; level<IMAGES-1; level++) {                // Loop over sublevels of tree
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            if (ix != 0 || iy != 0 || iz != 0) {                //     If periodic cell is not at center
              for (int cx=-1; cx<=1; cx++) {                    //      Loop over x periodic direction (child)
                for (int cy=-1; cy<=1; cy++) {                  //       Loop over y periodic direction (child)
                  for (int cz=-1; cz<=1; cz++) {                //        Loop over z periodic direction (child)
                    Xperiodic[0] = (ix * 3 + cx) * Length;      //         Coordinate offset for x periodic direction
                    Xperiodic[1] = (iy * 3 + cy) * Length;      //         Coordinate offset for y periodic direction
                    Xperiodic[2] = (iz * 3 + cz) * Length;      //         Coordinate offset for z periodic direction
                    M2L(Ci0,Ci,false);                          //         Perform M2L kernel
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cj0 = pcells.begin();                                     //  Redefine Cj0 for M2M
      C_iter Cj = Cj0;                                          //  Iterator of periodic neighbor cells
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++, Cj++) {                  //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              Cj->X[0] = Ci->X[0] + ix * Length;                //      Set new x coordinate for periodic image
              Cj->X[1] = Ci->X[0] + iy * Length;                //      Set new y cooridnate for periodic image
              Cj->X[2] = Ci->X[0] + iz * Length;                //      Set new z coordinate for periodic image
              Cj->M    = Ci->M;                                 //      Copy multipoles to new periodic image
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Ci->M = 0;                                                //  Reset multipoles of periodic parent
      setCenter(Ci,Cj0);                                        //  Set center of mass for periodic parent
      M2M(Ci,Cj0);                                              //  Evaluate periodic M2M kernels for this sublevel
      Length *= 3;                                              //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
    stopTimer("Traverse periodic",printNow);                    // Stop timer
  }

 protected:
//! Set center of expansion to center of mass
  void setCenter(C_iter C, C_iter C0) const {
    real_t m = 0;                                               // Initialize mass
    vec3 X = 0;                                                 // Initialize coordinates
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {         // Loop over bodies
      m += B->SRC;                                              //  Accumulate mass
      X += B->X * B->SRC;                                       //  Accumulate dipole
    }                                                           // End loop over bodies
    for (C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; c++) { // Loop over child cells
      m += std::abs(c->M[0]);                                   //  Accumulate mass
      X += c->X * std::abs(c->M[0]);                            //  Accumulate dipole
    }                                                           // End loop over child cells
    X /= m;                                                     // Center of mass
#if USE_BMAX
    C->R = getBmax(X,C);                                        // Use Bmax as cell radius
#endif
#if COMcenter
    C->X = X;                                                   // Use center of mass as center of expansion
#endif
    C->RMAX = 0;                                                // Initialize Rmax
  }

 public:
  Traversal() : NP2P(0), NM2L(0) {}
  ~Traversal() {}

//! Evaluate P2P and M2L using dual tree traversal
  void dualTreeTraversal(Cells &icells, Cells &jcells, bool mutual=false) {
    Ci0 = icells.begin();                                       // Set iterator of target root cell
    Cj0 = jcells.begin();                                       // Set iterator of source root cell
    startTimer("Traverse");                                     // Start timer
    if (IMAGES == 0) {                                          // If non-periodic boundary condition
      Xperiodic = 0;                                            //  No periodic shift
      traverse(Ci0,Cj0,mutual);                                 //  Traverse the tree
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            Xperiodic[0] = ix * periodicCycle;                  //     Coordinate shift for x periodic direction
            Xperiodic[1] = iy * periodicCycle;                  //     Coordinate shift for y periodic direction
            Xperiodic[2] = iz * periodicCycle;                  //     Coordinate shift for z periodic direction
            traverse(Ci0,Cj0,false);                            //     Traverse the tree for this periodic image
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic(periodicCycle);                          //  Traverse tree for periodic images
    }                                                           // End if for periodic boundary condition
    stopTimer("Traverse",printNow);                             // Stop timer
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
    startTimer("P2P kernel");
    P2P(Ci,Cj,false);
    timeP2P = stopTimer("P2P kernel") / 10000;
    startTimer("M2L kernel");
    for (int i=0; i<1000; i++) M2L(Ci,Cj,false);
    timeM2L = stopTimer("M2L kernel") / 1000;
  }

//! Print traversal statistics
  void printTraversalData() {
#if COUNT
    std::cout << "--- Traversal stats --------------" << std::endl
	      << std::setw(stringLength) << std::left           // Set format
	      << "P2P calls"  << " : " << NP2P << std::endl     // Print number of P2P calls
	      << std::setw(stringLength) << std::left           // Set format
	      << "M2L calls"  << " : " << NM2L << std::endl;    // Print number of M2l calls
#endif
  }
};
#endif
