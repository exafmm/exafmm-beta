#ifndef updownpass_h
#define updownpass_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

class UpDownPass : public Kernel, public Logger {
 public:
  int IMAGES;                                                   //!< Number of periodic image sublevels
  real_t THETA;                                                 //!< Multipole acceptance criteria

 private:
//! Error optimization of Rcrit
  void setRcrit(C_iter C, C_iter C0, real_t c) {
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(setRcrit(CC, C0, c));                       //  Recursive call with new task
      }                                                         // End loop over child cells
      sync_tasks;                                               // Synchronize tasks
    }
#if Cartesian
    for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];            // Normalize multipole expansion coefficients
#endif
    real_t x = 1.0 / THETA;                                     // Inverse of theta
#if ERROR_OPT
    assert(THETA != 1.0);
    real_t a = c * powf(std::abs(C->M[0]),1.0/3);               // Cell coefficient
    for (int i=0; i<5; i++) {                                   // Newton-Rhapson iteration
      real_t f = x * x - 2 * x + 1 - a * pow(x,-P);             //  Function value
      real_t df = (P + 2) * x - 2 * (P + 1) + P / x;            //  Function derivative value
      x -= f / df;                                              //  Increment x
    }                                                           // End Newton-Rhapson iteration
#endif
    C->RCRIT *= x;                                              // Multiply Rcrit by error optimized parameter x
  }

//! Recursive call for upward pass
  void upwardRecursion(C_iter C, C_iter C0) {
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(upwardRecursion(CC, C0));                   //  Recursive call with new task
      }                                                         // End loop over child cells
      sync_tasks;                                               // Synchronize tasks
    }
    setCenter(C,C0);                                            // Set center of cell to center of mass
    C->M = 0;                                                   // Initialize multipole expansion coefficients
    C->L = 0;                                                   // Initialize local expansion coefficients
    P2M(C);                                                     // P2M kernel
    M2M(C,C0);                                                  // M2M kernel
  }

//! Recursive call for downward pass
  void downwardRecursion(C_iter C, C_iter C0) const {
    L2L(C,C0);                                                  // L2L kernel
    L2P(C);                                                     // L2P kernel
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(downwardRecursion(CC, C0));                 //  Recursive call with new task
      }                                                         // End loop over chlid cells
      sync_tasks;                                               // Synchronize tasks
    }
  }

 public:
  UpDownPass(int images, real_t theta) : IMAGES(images), THETA(theta) {}
  ~UpDownPass() {}

//! Upward pass (P2M, M2M)
  void upwardPass(Cells &cells) {
    startTimer("Upward pass");                                  // Start timer
    C_iter C0 = cells.begin();                                  // Set iterator of target root cell
    upwardRecursion(C0, C0);                                    // Recursive call for upward pass
    real_t c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / powf(std::abs(C0->M[0]),1.0/3); // Root coefficient
    setRcrit(C0, C0, c);                                        // Error optimization of Rcrit
    if( cells.size() > 9 ) {                                    // If tree has more than 2 levels
      for (C_iter C=C0; C!=C0+9; C++) {                         //  Loop over top 2 levels of cells
        C->RCRIT *= 10;                                         //   Prevent approximation
      }                                                         //  End loop over top 2 levels of cells
    }                                                           // End if for tree levels
    stopTimer("Upward pass",printNow);                          // Stop timer
  }

//! Downward pass (L2L, L2P)
  void downwardPass(Cells &cells) { 
    startTimer("Downward pass");                                // Start timer
    C_iter C0 = cells.begin();                                  // Root cell
    L2P(C0);                                                    // If root is the only cell do L2P
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C0->CHILD; CC!=C0+C0->CHILD+C0->NCHILD; CC++) {// Loop over child cells
	spawn_task0(downwardRecursion(CC, C0));                 //  Recursive call for downward pass
      }                                                         // End loop over child cells
      sync_tasks;                                               // Synchronize tasks
    }
    stopTimer("Downward pass",printNow);                        // Stop timer
  }

  void direct(Bodies &ibodies, Bodies &jbodies, real_t cycle) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;            // First cell is target, second cell is source
    Ci->BODY = ibodies.begin();                                 // Iterator of first target body
    Ci->NDBODY = ibodies.size();                                // Number of target bodies
    Cj->BODY = jbodies.begin();                                 // Iterator of first source body
    Cj->NDBODY = jbodies.size();                                // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<IMAGES; i++) {                              // Loop over periodic image sublevels
      prange += int(powf(3,i));                                 //  Accumulate range of periodic images
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
};
#endif
