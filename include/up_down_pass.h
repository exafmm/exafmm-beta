#ifndef up_down_pass_h
#define up_down_pass_h
#include "logger.h"
#include "thread.h"
#include "types.h"

namespace exafmm {
  template<typename Kernel>
  class UpDownPass {
    typedef typename Kernel::Bodies Bodies;                     //!< Vector of bodies
    typedef typename Kernel::Cells Cells;                       //!< Vector of cells
    typedef typename Kernel::B_iter B_iter;                     //!< Iterator of body vector
    typedef typename Kernel::C_iter C_iter;                     //!< Iterator of cell vecto
    static const int P = Kernel::P;                             //!< Set order of expansion

  private:
    const real_t theta;                                         //!< Multipole acceptance criteria

  private:
    //! Recursive functor for setting cell scale
    struct SetScaleFromRadius {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      SetScaleFromRadius(C_iter _C, C_iter _C0) :               // Constructor
	C(_C), C0(_C0) {}                                       // Initialize variables
      void operator() () const {                                // Overload operator()
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
      PostOrderTraversal(C_iter _C, C_iter _C0, real_t _theta) : // Constructor
	C(_C), C0(_C0), theta(_theta) {}     // Initialize variables
      void operator() () const {                                // Overload operator()
	mk_task_group;                                          //  Initialize tasks
        for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) { // Loop over child cells
	  PostOrderTraversal postOrderTraversal(CC, C0, theta); // Instantiate recursive functor
	  create_taskc(postOrderTraversal);                     //    Create new task for recursive call
	}                                                       //   End loop over child cells
	wait_tasks;                                             //   Synchronize tasks
	C->M = 0;                                               //  Initialize multipole expansion coefficients
	C->L = 0;                                               //  Initialize local expansion coefficients
	if(C->NCHILD==0) Kernel::P2M(C);                        //  P2M kernel
	else {                                                  //  If not leaf cell
          Kernel::M2M(C, C0);                                   //   M2M kernel
        }                                                       //  End if for non leaf cell
	C->R /= theta;                                          //  Divide R by theta
      }                                                         // End overload operator()
    };

    //! Recursive functor for the pre-order traversal during downward pass
    struct PreOrderTraversal {
      C_iter C;                                                 //!< Iterator of current cell
      C_iter C0;                                                //!< Iterator of first cell
      PreOrderTraversal(C_iter _C, C_iter _C0) :                // Constructor
	C(_C), C0(_C0) {}                                       // Initialize variables
      void operator() () const {                                // Overload operator()
	Kernel::L2L(C, C0);                                     //  L2L kernel
	if (C->NCHILD==0) {                                     //  If leaf cell
          Kernel::L2P(C);                                       //  L2P kernel
        }                                                       // End if for leaf cell
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
    UpDownPass(real_t _theta) :
      theta(_theta) {     // Initialize variables
    }

    //! Upward pass (P2M, M2M)
    void upwardPass(Cells & cells) {
      logger::startTimer("Upward pass");                        // Start timer
      if (!cells.empty()) {                                     // If cell vector is not empty
	C_iter C0 = cells.begin();                              //  Set iterator of target root cell
	SetScaleFromRadius setScaleFromRadius(C0, C0);          //  Instantiate recursive functor
	setScaleFromRadius();                                   //  Recursive call for setting cell scale
	PostOrderTraversal postOrderTraversal(C0, C0, theta);   // Instantiate recursive functor
	postOrderTraversal();                                   //  Recursive call for upward pass
	real_t c = (1 - theta) * (1 - theta) / std::pow(theta,P+2) / powf(std::abs(C0->M[0]),1.0/3); // Root coefficient
      }                                                         // End if for empty cell vector
      logger::stopTimer("Upward pass");                         // Stop timer
    }

    //! Downward pass (L2L, L2P)
    void downwardPass(Cells & cells) {
      logger::startTimer("Downward pass");                      // Start timer
      if (!cells.empty()) {                                     // If cell vector is not empty
	C_iter C0 = cells.begin();                              //  Root cell
	if (C0->NCHILD == 0 ) {                                 //  If root is the only cell
          Kernel::L2P(C0);                                      //   L2P kernel
        }                                                       //  End if root is the only cell
	mk_task_group;                                          //  Initialize tasks
	for (C_iter CC=C0+C0->ICHILD; CC!=C0+C0->ICHILD+C0->NCHILD; CC++) {// Loop over child cells
	  PreOrderTraversal preOrderTraversal(CC, C0);          //    Instantiate recursive functor
	  create_taskc(preOrderTraversal);                      //    Recursive call for downward pass
	}                                                       //   End loop over child cells
	wait_tasks;                                             //   Synchronize tasks
      }                                                         // End if for empty cell vector
      logger::stopTimer("Downward pass");                       // Stop timer
    }

    //! Get dipole of entire system
    vec3 getDipole(Bodies & bodies, vec3 X0) {
      vec3 dipole = 0;                                          // Initialize dipole correction
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	dipole += (B->X - X0) * std::real(complex_t(B->SRC));   //  Calcuate dipole of the whole system
      }                                                         // End loop over bodies
      return dipole;                                            // Return dipole
    }

    //! Dipole correction
    void dipoleCorrection(Bodies & bodies, vec3 dipole, int numBodies, vec3 cycle) {
      real_t coef = 4 * M_PI / (3 * cycle[0] * cycle[1] * cycle[2]);// Precalcualte constant
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	B->TRG[0] -= coef * norm(dipole) / numBodies / B->SRC;  //  Dipole correction for potential
	for (int d=0; d!=3; d++) {                              //  Loop over dimensions
	  B->TRG[d+1] -= coef * dipole[d];                      //   Dipole correction for forces
	}                                                       //  End loop over dimensions
      }                                                         // End loop over bodies
    }
  };
}
#endif
