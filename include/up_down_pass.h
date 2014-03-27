#ifndef up_down_pass_h
#define up_down_pass_h
#include "kernel.h"
#include "logger.h"
#include <tpswitch/tpswitch.h>

class UpDownPass : public Kernel {
public:
  real_t theta;                                                 //!< Multipole acceptance criteria

private:
  //! Recursive functor for error optimization of Rcrit
  struct SetRcrit {
    C_iter C;                                                   //!< Iterator of current cell
    C_iter C0;                                                  //!< Iterator of first cell
    real_t c;                                                   //!< Root coefficient
    real_t theta;                                               //!< Multipole acceptance criteria
    SetRcrit(C_iter _C, C_iter _C0, real_t _c, real_t _theta) : // Constructor
      C(_C), C0(_C0), c(_c), theta(_theta) {}                   // Initialize variables
    void operator() () {                                        // Overload operator()
      mk_task_group;                                            //  Initialize tasks
      for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	SetRcrit setRcrit(CC, C0, c, theta);                    //   Initialize recusive functor
	create_taskc(setRcrit);                                 //   Create new task for recursive call
      }                                                         //  End loop over child cells
      wait_tasks;                                               //  Synchronize tasks
#if Cartesian
      for (int i=1; i<NTERM; i++) C->M[i] /= C->M[0];           //  Normalize multipole expansion coefficients
#endif
      real_t x = 1.0 / theta;                                   //  Inverse of theta
#if ERROR_OPT
      assert(theta != 1.0);                                     //  Newton-Raphson won't work for theta==1
      real_t a = c * powf(std::abs(C->M[0]),1.0/3);             //  Cell coefficient
      for (int i=0; i<5; i++) {                                 //  Newton-Raphson iteration
	real_t f = x * x - 2 * x + 1 - a * std::pow(x,-P);      //   Function value
	real_t df = (P + 2) * x - 2 * (P + 1) + P / x;          //   Function derivative value
	x -= f / df;                                            //   Increment x
      }                                                         //  End Newton-Raphson iteration
#endif
      C->RCRIT *= x;                                            //  Multiply Rcrit by error optimized parameter x
    }                                                           // End overload operator()
  };

  //! Recursive functor for the post-order traversal during upward pass
  struct PostOrderTraversal : public Kernel {
    C_iter C;                                                   //!< Iterator of current cell
    C_iter C0;                                                  //!< Iterator of first cell
    PostOrderTraversal(C_iter _C, C_iter _C0) :                 // Constructor
      C(_C), C0(_C0) {}                                         // Initialize variables
    void operator() () {                                        // Overload operator()
      mk_task_group;                                            //  Initialize tasks
      for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	PostOrderTraversal postOrderTraversal(CC, C0);          //    Instantiate recursive functor
	create_taskc(postOrderTraversal);                       //    Create new task for recursive call
      }                                                         //   End loop over child cells
      wait_tasks;                                               //   Synchronize tasks
      C->RMAX = 0;                                              //  Initialzie Rmax
      C->M = 0;                                                 //  Initialize multipole expansion coefficients
      C->L = 0;                                                 //  Initialize local expansion coefficients
      if(C->NCHILD==0) P2M(C);                                  //  P2M kernel
      else M2M(C,C0);                                           //  M2M kernel
    }                                                           // End overload operator()
  };

  //! Recursive functor for the pre-order traversal during downward pass
  struct PreOrderTraversal : public Kernel {
    C_iter C;                                                   //!< Iterator of current cell
    C_iter C0;                                                  //!< Iterator of first cell
    PreOrderTraversal(C_iter _C, C_iter _C0) :                  // Constructor
      C(_C), C0(_C0) {}                                         // Initialize variables
    void operator() () {                                        // Overload operator()
      L2L(C,C0);                                                //  L2L kernel
      if (C->NCHILD==0) L2P(C);                                 //  L2P kernel
      mk_task_group;                                            //  Initialize tasks
      for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	PreOrderTraversal preOrderTraversal(CC, C0);            //   Instantiate recursive functor
	create_taskc(preOrderTraversal);                        //   Create new task for recursive call
      }                                                         //  End loop over chlid cells
      wait_tasks;                                               //  Synchronize tasks
    }                                                           // End overload operator()
  };

public:
  //! Constructor
  UpDownPass(real_t _theta) : theta(_theta) {}                  // Initialize variables

  //! Upward pass (P2M, M2M)
  void upwardPass(Cells & cells) {
    logger::startTimer("Upward pass");                          // Start timer
    if (!cells.empty()) {                                       // If cell vector is not empty
      C_iter C0 = cells.begin();                                //  Set iterator of target root cell
      PostOrderTraversal postOrderTraversal(C0, C0);            //  Instantiate recursive functor
      postOrderTraversal();                                     //  Recursive call for upward pass
      real_t c = (1 - theta) * (1 - theta) / std::pow(theta,P+2) / powf(std::abs(C0->M[0]),1.0/3); // Root coefficient
      SetRcrit setRcrit(C0, C0, c, theta);                      //  Instantiate recursive functor
      setRcrit();                                               //  Error optimization of Rcrit
      if( cells.size() > 9 ) {                                  //  If tree has more than 2 levels
        for (C_iter C=C0; C!=C0+9; C++) {                       //   Loop over top 2 levels of cells
          C->RCRIT = 1e12;                                      //    Prevent approximation
        }                                                       //   End loop over top 2 levels of cells
      }                                                         //  End if for tree levels
    }                                                           // End if for empty cell vector
    logger::stopTimer("Upward pass");                           // Stop timer
  }

  //! Downward pass (L2L, L2P)
  void downwardPass(Cells & cells) {
    logger::startTimer("Downward pass");                        // Start timer
    if (!cells.empty()) {                                       // If cell vector is not empty
      C_iter C0 = cells.begin();                                //  Root cell
      if (C0->NCHILD == 0) L2P(C0);                             //  If root is the only cell do L2P
      mk_task_group;                                            //  Initialize tasks
      for (C_iter CC=C0+C0->ICHILD; CC!=C0+C0->ICHILD+C0->NCHILD; CC++) {// Loop over child cells
	PreOrderTraversal preOrderTraversal(CC, C0);            //    Instantiate recursive functor
	create_taskc(preOrderTraversal);                        //    Recursive call for downward pass
      }                                                         //   End loop over child cells
      wait_tasks;                                               //   Synchronize tasks
    }                                                           // End if for empty cell vector
    logger::stopTimer("Downward pass");                         // Stop timer
  }

  //! Get dipole of entire system
  vec3 getDipole(Bodies & bodies, vec3 X0) {
    vec3 dipole = 0;                                            // Initialize dipole correction
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      dipole += (B->X - X0) * B->SRC;                           //  Calcuate dipole of the whole system
    }                                                           // End loop over bodies
    return dipole;                                              // Return dipole
  }

  //! Dipole correction
  void dipoleCorrection(Bodies & bodies, vec3 dipole, int numBodies, real_t cycle) {
    real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle);       // Precalcualte constant
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG[0] -= coef * norm(dipole) / numBodies / B->SRC;    //  Dipole correction for potential
      for (int d=0; d!=3; d++) {                                //  Loop over dimensions
        B->TRG[d+1] -= coef * dipole[d];                        //   Dipole correction for forces
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
  }
};
#endif
