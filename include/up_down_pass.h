#ifndef up_down_pass_h
#define up_down_pass_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

namespace exafmm {
  class UpDownPass {
  private:
    const real_t theta;                                         //!< Multipole acceptance criteria
    const bool useRmax;                                         //!< Use maximum distance for MAC
    const bool useRopt;                                         //!< Use error optimized theta for MAC

  private:
    //! Recursive functor for error optimization of R
    struct SetRopt {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      real_t c;                                                 //!< Root coefficient
      real_t theta;                                             //!< Multipole acceptance criteria
      SetRopt(C_iter _C, C_iter _C0, real_t _c, real_t _theta) :// Constructor
	C(_C), C0(_C0), c(_c), theta(_theta) {}                 // Initialize variables
      void operator() () {                                      // Overload operator()
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	  SetRopt setRopt(CC, C0, c, theta);                    //   Initialize recusive functor
	  create_taskc(setRopt);                                //   Create new task for recursive call
	}                                                       //  End loop over child cells
	wait_tasks;                                             //  Synchronize tasks
#if EXAFMM_MASS
	for (int i=1; i<NTERM; i++) C->M[i] /= C->M[0];         //  Normalize multipole expansion coefficients
#endif
	real_t x = 1.0 / theta;                                 //  Inverse of theta
	assert(theta != 1.0);                                   //  Newton-Raphson won't work for theta==1
	real_t a = c * powf(std::abs(C->M[0]),1.0/3);           //  Cell coefficient
	for (int i=0; i<5; i++) {                               //  Loop for Newton-Raphson iteration
	  real_t f = x * x - 2 * x + 1 - a * std::pow(x,-P);    //   Function value
	  real_t df = (P + 2) * x - 2 * (P + 1) + P / x;        //   Function derivative value
	  x -= f / df;                                          //   Increment x
	}                                                       //  End loop for Newton-Raphson iteration
	C->R *= x * theta;                                      //  Multiply R by error optimized parameter x
      }                                                         // End overload operator()
    };

    //! Recursive functor for setting cell scale
    struct SetScaleFromRadius {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      SetScaleFromRadius(C_iter _C, C_iter _C0) :               // Constructor
	C(_C), C0(_C0) {}                                       // Initialize variables
      void operator() () {                                      // Overload operator()
	C->SCALE = 2 * C->R;                                    //  Set cell scale
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	  SetScaleFromRadius setScaleFromRadius(CC, C0);        //   Instantiate recursive functor
	  create_taskc(setScaleFromRadius);                     //   Create new task for recursive call
	}                                                       //  End loop over chlid cells
	wait_tasks;                                             //  Synchronize tasks
      }                                                         // End overload operator()
    };

    //! Recursive functor for the post-order traversal during upward pass
    struct PostOrderTraversal {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      real_t theta;                                             //!< Multipole acceptance criteria
      bool useRmax;                                             //!< Use maximum distance for MAC
      //! Redefine cell radius R based on maximum distance
      void setRmax() {
	real_t Rmax = 0;                                        // Initialize Rmax
	if (C->NCHILD == 0) {                                   // If leaf cell
	  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {    //  Loop over bodies in cell
	    vec3 dX = C->X - B->X;                              //   Distance vector from particles to center of expansion
	    real_t R = std::sqrt(norm(dX));                     //   Scalar distance
	    if (R > Rmax) Rmax = R;                             //   Maximum distance
	  }                                                     //  End loop over bodies in cell
	} else {                                                // If not leaf cell
	  for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	    vec3 dX = C->X - CC->X;                             //   Distance vector from particles to center of expansion
	    real_t R = std::sqrt(norm(dX)) + CC->R;             //   Scalar distance
	    if (R > Rmax) Rmax = R;                             //   Maximum distance
	  }                                                     //  End loop over child cells
	}                                                       // End if for leaf cell
	C->R = std::min(C->R,Rmax);                             // Redefine R based on maximum distance
      }
      PostOrderTraversal(C_iter _C, C_iter _C0, real_t _theta, bool _useRmax) : // Constructor
	C(_C), C0(_C0), theta(_theta), useRmax(_useRmax) {}     // Initialize variables
      void operator() () {                                      // Overload operator()
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	  PostOrderTraversal postOrderTraversal(CC, C0, theta, useRmax); // Instantiate recursive functor
	  create_taskc(postOrderTraversal);                     //    Create new task for recursive call
	}                                                       //   End loop over child cells
	wait_tasks;                                             //   Synchronize tasks
	C->M = 0;                                               //  Initialize multipole expansion coefficients
	C->L = 0;                                               //  Initialize local expansion coefficients
	if(C->NCHILD==0) kernel::P2M(C);                        //  P2M kernel
	else kernel::M2M(C, C0);                                //  M2M kernel
	if (useRmax) setRmax();                                 //  Redefine cell radius R based on maximum distance
	C->R /= theta;                                          //  Divide R by theta
      }                                                         // End overload operator()
    };

    //! Recursive functor for the pre-order traversal during downward pass
    struct PreOrderTraversal {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      PreOrderTraversal(C_iter _C, C_iter _C0) :                // Constructor
	C(_C), C0(_C0) {}                                       // Initialize variables
      void operator() () {                                      // Overload operator()
	kernel::L2L(C, C0);                                     //  L2L kernel
	if (C->NCHILD==0) kernel::L2P(C);                       //  L2P kernel
#if EXAFMM_USE_WEIGHT
	C_iter CP = C0 + C->IPARENT;                            // Parent cell
	C->WEIGHT += CP->WEIGHT;                                // Add parent's weight
	if (C->NCHILD==0) {                                     // If leaf cell
	  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {    //  Loop over bodies in cell
	    B->WEIGHT += C->WEIGHT;                             //   Add cell weights to bodies
	  }                                                     //  End loop over bodies in cell
	}                                                       // End if for leaf cell
#endif
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
	  PreOrderTraversal preOrderTraversal(CC, C0);          //   Instantiate recursive functor
	  create_taskc(preOrderTraversal);                      //   Create new task for recursive call
	}                                                       //  End loop over chlid cells
	wait_tasks;                                             //  Synchronize tasks
      }                                                         // End overload operator()
    };

  public:
    //! Constructor
    UpDownPass(real_t _theta, bool _useRmax, bool _useRopt) :
      theta(_theta), useRmax(_useRmax), useRopt(_useRopt) {     // Initialize variables
    }

    //! Upward pass (P2M, M2M)
    void upwardPass(Cells & cells) {
      logger::startTimer("Upward pass");                        // Start timer
      if (!cells.empty()) {                                     // If cell vector is not empty
	C_iter C0 = cells.begin();                              //  Set iterator of target root cell
	SetScaleFromRadius setScaleFromRadius(C0, C0);          //  Instantiate recursive functor
	setScaleFromRadius();                                   //  Recursive call for setting cell scale
	PostOrderTraversal postOrderTraversal(C0, C0, theta, useRmax); // Instantiate recursive functor
	postOrderTraversal();                                   //  Recursive call for upward pass
	real_t c = (1 - theta) * (1 - theta) / std::pow(theta,P+2) / powf(std::abs(C0->M[0]),1.0/3); // Root coefficient
	if (useRopt) {                                          //  If using error optimized theta
	  SetRopt setRopt(C0, C0, c, theta);                    //   Instantiate recursive functor
	  setRopt();                                            //   Error optimization of R
	}                                                       //  End if for using error optimized theta
      }                                                         // End if for empty cell vector
      logger::stopTimer("Upward pass");                         // Stop timer
    }

    //! Downward pass (L2L, L2P)
    void downwardPass(Cells & cells) {
      logger::startTimer("Downward pass");                      // Start timer
      if (!cells.empty()) {                                     // If cell vector is not empty
	C_iter C0 = cells.begin();                              //  Root cell
	if (C0->NCHILD == 0) kernel::L2P(C0);                   //  If root is the only cell do L2P
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C0->ICHILD; CC!=C0+C0->ICHILD+C0->NCHILD; CC++) {// Loop over child cells
	  PreOrderTraversal preOrderTraversal(CC, C0);          //    Instantiate recursive functor
	  create_taskc(preOrderTraversal);                      //    Recursive call for downward pass
	}                                                       //   End loop over child cells
	wait_tasks;                                             //   Synchronize tasks
      }                                                         // End if for empty cell vector
      logger::stopTimer("Downward pass");                       // Stop timer
    }

#if EXAFMM_LAPLACE
    //! Get dipole of entire system
    vec3 getDipole(Bodies & bodies, vec3 X0) {
      vec3 dipole = 0;                                          // Initialize dipole correction
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	dipole += (B->X - X0) * std::real(complex_t(B->SRC));   //  Calcuate dipole of the whole system
      }                                                         // End loop over bodies
      return dipole;                                            // Return dipole
    }

    //! Dipole correction
    void dipoleCorrection(Bodies & bodies, vec3 dipole, int numBodies, real_t cycle) {
      real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle);     // Precalcualte constant
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	B->TRG[0] -= coef * norm(dipole) / numBodies / B->SRC;  //  Dipole correction for potential
	for (int d=0; d!=3; d++) {                              //  Loop over dimensions
	  B->TRG[d+1] -= coef * dipole[d];                      //   Dipole correction for forces
	}                                                       //  End loop over dimensions
      }                                                         // End loop over bodies
    }
#endif
  };
}
#endif
