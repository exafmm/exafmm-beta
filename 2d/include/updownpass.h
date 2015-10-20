#ifndef updownpass_h
#define updownpass_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

class UpDownPass : public Kernel, public Logger {
 public:
  real_t theta;                                                 //!< Multipole acceptance criteria

 private:
//! Error optimization of Rcrit
  void setRcrit(C_iter C, C_iter C0, real_t c) {
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(setRcrit(CC, C0, c));                       //  Recursive call with new task
      }                                                         // End loop over child cells
      sync_tasks;                                               // Synchronize tasks
    }
    C->RCRIT = C->R / theta;                                    // Multiply Rcrit by error optimized parameter x
  }

//! Recursive call for upward pass
  void postOrderTraversal(C_iter C, C_iter C0) {
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(postOrderTraversal(CC, C0));                //   Recursive call with new task
      }                                                         //  End loop over child cells
      sync_tasks;                                               //  Synchronize tasks
    }                                                           // Finalize tasks
    C->M = 0;                                                   // Initialize multipole expansion coefficients
    C->L = 0;                                                   // Initialize local expansion coefficients
    P2M(C);                                                     // P2M kernel
    M2M(C,C0);                                                  // M2M kernel
  }

//! Recursive call for downward pass
  void preOrderTraversal(C_iter C, C_iter C0) const {
    L2L(C,C0);                                                  // L2L kernel
    L2P(C);                                                     // L2P kernel
    spawn_tasks {                                               // Initialize tasks
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
	spawn_task0(preOrderTraversal(CC, C0));                 //   Recursive call with new task
      }                                                         //  End loop over chlid cells
      sync_tasks;                                               //  Synchronize tasks
    }                                                           // Finalize tasks
  }

 public:
  UpDownPass(real_t _theta, real_t eps2) : Kernel(eps2), theta(_theta) {}

//! Upward pass (P2M, M2M)
  void upwardPass(Cells &cells) {
    startTimer("Upward pass");                                  // Start timer
    if (!cells.empty()) {                                       // If cell vector is not empty
      C_iter C0 = cells.begin();                                //  Set iterator of target root cell
      postOrderTraversal(C0, C0);                               //  Recursive call for upward pass
      real_t c = (1 - theta) * (1 - theta) / std::pow(theta,P+2) / powf(std::abs(C0->M[0]),1.0/3); // Root coefficient
      setRcrit(C0, C0, c);                                      //  Error optimization of Rcrit
      if( cells.size() > 9 ) {                                  //  If tree has more than 2 levels
        for (C_iter C=C0; C!=C0+9; C++) {                       //   Loop over top 2 levels of cells
          C->RCRIT = 1e12;                                      //    Prevent approximation
        }                                                       //   End loop over top 2 levels of cells
      }                                                         //  End if for tree levels
    }                                                           // End if for empty cell vector
    stopTimer("Upward pass");                                   // Stop timer
  }

//! Downward pass (L2L, L2P)
  void downwardPass(Cells &cells) {
    startTimer("Downward pass");                                // Start timer
    if (!cells.empty()) {                                       // If cell vector is not empty
      C_iter C0 = cells.begin();                                //  Root cell
      L2P(C0);                                                  //  If root is the only cell do L2P
      spawn_tasks {                                             //  Initialize tasks
        for (C_iter CC=C0+C0->CHILD; CC!=C0+C0->CHILD+C0->NCHILD; CC++) {// Loop over child cells
          spawn_task0(preOrderTraversal(CC, C0));               //    Recursive call for downward pass
        }                                                       //   End loop over child cells
        sync_tasks;                                             //   Synchronize tasks
      }                                                         //  Finalize tasks
    }                                                           // End if for empty cell vector
    stopTimer("Downward pass");                                 // Stop timer
  }
};
#endif
