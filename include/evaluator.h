#ifndef evaluator_h
#define evaluator_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"
#if COUNT
#define count(N) N++
#else
#define count(N)
#endif

class Evaluator : public Kernel, public Logger {
private:
  real_t timeP2P;                                               //!< P2P execution time
  real_t timeM2L;                                               //!< M2L execution time

protected:
  real_t NP2P;                                                  //!< Number of P2P kernel calls
  real_t NM2L;                                                  //!< Number of M2L kernel calls

public:
  int NSPAWN;                                                   //!< Threshold of NDBODY for spawning new threads
  int IMAGES;                                                   //!< Number of periodic image sublevels
  float THETA;                                                  //!< Multipole acceptance criteria

private:
//! Calculate Bmax
  real_t getBmax(vec3 const &X, C_iter C) const {
    real_t rad = C->R;                                          // Radius of cell
    real_t dx = rad+std::abs(X[0]-C->X[0]);                     // Add x distance from center of mass
    real_t dy = rad+std::abs(X[1]-C->X[1]);                     // Add y distance from center of mass
    real_t dz = rad+std::abs(X[2]-C->X[2]);                     // Add z distance from center of mass
    return std::sqrt( dx*dx + dy*dy + dz*dz );                  // Return scalar distance
  }

//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj, bool mutual) {
#if AUTO
    if (timeP2P*Ci->NDBODY*Cj->NDBODY > timeM2L) {              // If M2L is faster
      M2L(Ci,Cj,mutual);                                        //  M2L kernel
      count(NM2L);                                              //  Increment M2L counter
    } else {                                                    // Else if P2P is faster
      P2P(Ci,Cj,mutual);                                        //  P2P kernel
      count(NP2P);                                              //  Increment P2P counter
    }                                                           // End if for fastest kernel
#else
    M2L(Ci,Cj,mutual);                                          // M2L kernel
    count(NM2L);                                                // Increment M2L counter
#endif
  }

  void traverse(C_iter CiBegin, C_iter CiEnd, C_iter CjBegin, C_iter CjEnd, bool mutual) {
    if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {
      if (CiBegin == CjBegin) {                                 // TODO : unecessary?
        assert(CiEnd == CjEnd);
        traverse(CiBegin, CjBegin, mutual);
      } else {
        for (C_iter Ci=CiBegin; Ci!=CiEnd; Ci++) {
          for (C_iter Cj=CjBegin; Cj!=CjEnd; Cj++) {
            traverse(Ci, Cj, mutual);
          }
        }
      }
    } else {
      C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;
      C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;
      __init_tasks__;
      spawn_task0(traverse(CiBegin, CiMid, CjBegin, CjMid, mutual));
      traverse(CiMid, CiEnd, CjMid, CjEnd, mutual);
      __sync_tasks__;
      spawn_task0(traverse(CiBegin, CiMid, CjMid, CjEnd, mutual));
      if (!mutual || CiBegin != CjBegin) {
        traverse(CiMid, CiEnd, CjBegin, CjMid, mutual);
      } else {
        assert(CiEnd == CjEnd);
      }
      __sync_tasks__;
    }
  }

  void splitCell(C_iter Ci, C_iter Cj, bool mutual) {
    if (Cj->NCHILD == 0) {
      assert(Ci->NCHILD > 0);
      for (C_iter ci=Ci0+Ci->CHILD; ci!=Ci0+Ci->CHILD+Ci->NCHILD; ci++ ) {
        traverse(ci,Cj,mutual);
      }
    } else if (Ci->NCHILD == 0) {
      assert(Cj->NCHILD > 0);
      for (C_iter cj=Cj0+Cj->CHILD; cj!=Cj0+Cj->CHILD+Cj->NCHILD; cj++ ) {
        traverse(Ci,cj,mutual);
      }
    } else if (Ci->NDBODY + Cj->NDBODY >= NSPAWN || (mutual && Ci == Cj)) {
      traverse(Ci0+Ci->CHILD, Ci0+Ci->CHILD+Ci->NCHILD,
               Cj0+Cj->CHILD, Cj0+Cj->CHILD+Cj->NCHILD, mutual);
    } else if (Ci->RCRIT >= Cj->RCRIT) {
      for (C_iter ci=Ci0+Ci->CHILD; ci!=Ci0+Ci->CHILD+Ci->NCHILD; ci++ ) {
        traverse(ci,Cj,mutual);
      } 
    } else {
      for (C_iter cj=Cj0+Cj->CHILD; cj!=Cj0+Cj->CHILD+Cj->NCHILD; cj++ ) {
        traverse(Ci,cj,mutual);
      }
    }
  }

protected:
//! Set center of expansion to center of mass
  void setCenter(C_iter C) const {
    real_t m = 0;                                               // Initialize mass
    vec3 X = 0;                                                 // Initialize coordinates
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {         // Loop over bodies
      m += B->SRC;                                              //  Accumulate mass
      X += B->X * B->SRC;                                       //  Accumulate dipole
    }                                                           // End loop over bodies
    for (C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; c++) {// Loop over child cells
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
  }

//! Dual tree traversal
  void traverse(C_iter Ci, C_iter Cj, bool mutual) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real_t R2 = norm(dX);                                       // Scalar distance squared
#if DUAL
    {                                                           // Dummy bracket
#else
    if (Ci->RCRIT != Cj->RCRIT) {                               // If cell is not at the same level
      splitCell(Ci,Cj,mutual);                                  //  Split cell and call function recursively for child
    } else {                                                    // If we don't care if cell is not at the same level
#endif
      if (R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {   //  If distance is far enough
        approximate(Ci,Cj,mutual);                              //   Use approximate kernels
      } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {          //  Else if both cells are bodies
        if (Cj->NCBODY == 0) {                                  //   If the bodies weren't sent from remote node
          approximate(Ci,Cj,mutual);                            //    Use approximate kernels
        } else {                                                //   Else if the bodies were sent
          if (Ci == Cj) {                                       //    If source and target are same
            P2P(Ci);                                            //     P2P kernel for single cell
          } else {                                              //    Else if source and target are different
            P2P(Ci,Cj,mutual);                                  //     P2P kernel for pair of cells
          }                                                     //    End if for same source and target
          count(NP2P);                                          //    Increment P2P counter
        }                                                       //   End if for bodies
      } else {                                                  //  Else if cells are close but not bodies
        splitCell(Ci,Cj,mutual);                                //   Split cell and call function recursively for child
      }                                                         //  End if for multipole acceptance
    }                                                           // End if for same level cells
  }


//! Tree traversal of periodic cells
  void traversePeriodic(real_t R) {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Periodic coordinate offset
    Cells pcells(28);                                           // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->CHILD = 0;                                              // Child cells for periodic center cell
    Ci->NCHILD = 27;                                            // Number of child cells for periodic center cell
    C_iter C0 = Cj0;                                            // Placeholder for Cj0
    for (int level=0; level<IMAGES-1; level++) {                // Loop over sublevels of tree
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            if (ix != 0 || iy != 0 || iz != 0) {                //     If periodic cell is not at center
              for (int cx=-1; cx<=1; cx++) {                    //      Loop over x periodic direction (child)
                for (int cy=-1; cy<=1; cy++) {                  //       Loop over y periodic direction (child)
                  for (int cz=-1; cz<=1; cz++) {                //        Loop over z periodic direction (child)
                    Xperiodic[0] = (ix * 3 + cx) * 2 * R;       //         Coordinate offset for x periodic direction
                    Xperiodic[1] = (iy * 3 + cy) * 2 * R;       //         Coordinate offset for y periodic direction
                    Xperiodic[2] = (iz * 3 + cz) * 2 * R;       //         Coordinate offset for z periodic direction
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
            Cj->X[0] = Ci->X[0] + ix * 2 * R;                   //     Set new x coordinate for periodic image
            Cj->X[1] = Ci->X[0] + iy * 2 * R;                   //     Set new y cooridnate for periodic image
            Cj->X[2] = Ci->X[0] + iz * 2 * R;                   //     Set new z coordinate for periodic image
            Cj->M    = Ci->M;                                   //     Copy multipoles to new periodic image
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Ci->M = 0;                                                //  Reset multipoles of periodic parent
      real_t Rmax = 0;                                          //  Dummy parameter for calling M2M
      setCenter(Ci);                                            //  Set center of mass for periodic parent
      M2M(Ci,Rmax);                                             //  Evaluate periodic M2M kernels for this sublevel
      R *= 3;                                                   //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
    stopTimer("Traverse periodic",printNow);                    // Stop timer
  }

public:
  Evaluator() : NP2P(0), NM2L(0) {}
  ~Evaluator() {}

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
};

#endif
