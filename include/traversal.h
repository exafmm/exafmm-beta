#ifndef traversal_h
#define traversal_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

#if COUNT
#define countKernel(N) N++
#else
#define countKernel(N)
#endif

class Traversal {
private:
  int nspawn;                                                   //!< Threshold of NBODY for spawning new threads
  int images;                                                   //!< Number of periodic image sublevels
#if COUNT
  real_t numP2P;                                                //!< Number of P2P kernel calls
  real_t numM2L;                                                //!< Number of M2L kernel calls
#endif
  C_iter Ci0;                                                   //!< Iterator of first target cell
  C_iter Cj0;                                                   //!< Iterator of first source cell

private:
#if USE_WEIGHT
  //! Accumulate interaction weights of cells
  void countWeight(C_iter Ci, C_iter Cj, bool mutual, real_t weight) {
    Ci->WEIGHT += weight;                                       // Increment weight of target cell
    if (mutual) Cj->WEIGHT += weight;                           // Increment weight of source cell
  }
#else
  void countWeight(C_iter, C_iter, bool, real_t) {}
#endif

  //! Dual tree traversal for a single pair of cells
  void traverse(C_iter Ci, C_iter Cj, vec3 Xperiodic, bool mutual, real_t remote) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real_t R2 = norm(dX);                                       // Scalar distance squared
    if (R2 > (Ci->R+Cj->R) * (Ci->R+Cj->R)) {                   // If distance is far enough
      kernel::M2L(Ci, Cj, Xperiodic, mutual);                   //  M2L kernel
      countKernel(numM2L);                                      //  Increment M2L counter
      countWeight(Ci, Cj, mutual, remote);                      //  Increment M2L weight
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are bodies
      if (Cj->NBODY == 0) {                                     //  If the bodies weren't sent from remote node
	std::cout << Ci->ICELL << " " << Cj->ICELL << std::endl;
	kernel::M2L(Ci, Cj, Xperiodic, mutual);                 //   M2L kernel
        countKernel(numM2L);                                    //   Increment M2L counter
	countWeight(Ci, Cj, mutual, remote);                    //   Increment M2L weight
      } else {                                                  //  Else if the bodies were sent
	if (R2 == 0 && Ci == Cj) {                              //   If source and target are same
	  kernel::P2P(Ci);                                      //    P2P kernel for single cell
	} else {                                                //   Else if source and target are different
	  kernel::P2P(Ci, Cj, Xperiodic, mutual);               //    P2P kernel for pair of cells
	}                                                       //   End if for same source and target
	countKernel(numP2P);                                    //   Increment P2P counter
	countWeight(Ci, Cj, mutual, remote);                    //   Increment P2P weight
      }                                                         //  End if for bodies
    } else {                                                    // Else if cells are close but not bodies
      splitCell(Ci, Cj, Xperiodic, mutual, remote);             //  Split cell and call function recursively for child
    }                                                           // End if for multipole acceptance
  }

  //! Recursive functor for dual tree traversal of a range of Ci and Cj
  struct TraverseRange {
    Traversal * traversal;                                      //!< Traversal object
    C_iter CiBegin;                                             //!< Begin iterator of target cells
    C_iter CiEnd;                                               //!< End iterator of target cells
    C_iter CjBegin;                                             //!< Begin Iterator of source cells
    C_iter CjEnd;                                               //!< End iterator of source cells
    vec3 Xperiodic;                                             //!< Periodic coordinate offset
    bool mutual;                                                //!< Flag for mutual interaction
    real_t remote;                                              //!< Weight for remote work load
    TraverseRange(Traversal * _traversal, C_iter _CiBegin, C_iter _CiEnd,// Constructor
		  C_iter _CjBegin, C_iter _CjEnd, vec3 _Xperiodic, bool _mutual, real_t _remote) :
      traversal(_traversal), CiBegin(_CiBegin), CiEnd(_CiEnd),  // Initialize variables
      CjBegin(_CjBegin), CjEnd(_CjEnd), Xperiodic(_Xperiodic),
      mutual(_mutual), remote(_remote) {}
    void operator() () {                                        // Overload operator()
      Tracer tracer;                                            //  Instantiate tracer
      logger::startTracer(tracer);                              //  Start tracer
      if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {       //  If only one cell in range
	if (CiBegin == CjBegin) {                               //   If Ci == Cj
	  assert(CiEnd == CjEnd);                               //    Check if mutual & self interaction
	  traversal->traverse(CiBegin, CjBegin, Xperiodic, mutual, remote);// Call traverse for single pair
	} else {                                                //   If Ci != Cj
	  for (C_iter Ci=CiBegin; Ci!=CiEnd; Ci++) {            //    Loop over all Ci cells
	    for (C_iter Cj=CjBegin; Cj!=CjEnd; Cj++) {          //     Loop over all Cj cells
	      traversal->traverse(Ci, Cj, Xperiodic, mutual, remote); // Call traverse for single pair
	    }                                                   //     End loop over all Cj cells
	  }                                                     //    End loop over all Ci cells
	}                                                       //   End if for Ci == Cj
      } else {                                                  //  If many cells are in the range
	C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;         //   Split range of Ci cells in half
	C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;         //   Split range of Cj cells in half
	mk_task_group;                                          //   Initialize task group
	{
	  TraverseRange leftBranch(traversal, CiBegin, CiMid,   //    Instantiate recursive functor
				   CjBegin, CjMid, Xperiodic, mutual, remote);
	  create_taskc(leftBranch);                             //    Ci:former Cj:former
	  TraverseRange rightBranch(traversal, CiMid, CiEnd,    //    Instantiate recursive functor
				    CjMid, CjEnd, Xperiodic, mutual, remote);
	  rightBranch();                                        //    Ci:latter Cj:latter
	  wait_tasks;                                           //    Synchronize task group
	}
	{
	  TraverseRange leftBranch(traversal, CiBegin, CiMid,   //    Instantiate recursive functor
				   CjMid, CjEnd, Xperiodic, mutual, remote);
	  create_taskc(leftBranch);                             //    Ci:former Cj:latter
	  if (!mutual || CiBegin != CjBegin) {                  //    Exclude mutual & self interaction
            TraverseRange rightBranch(traversal, CiMid, CiEnd,  //    Instantiate recursive functor
				      CjBegin, CjMid, Xperiodic, mutual, remote);
	    rightBranch();                                      //    Ci:latter Cj:former
	  } else {                                              //    If mutual or self interaction
	    assert(CiEnd == CjEnd);                             //     Check if mutual & self interaction
	  }                                                     //    End if for mutual & self interaction
	  wait_tasks;                                           //    Synchronize task group
	}
      }                                                         //  End if for many cells in range
      logger::stopTracer(tracer);                               //  Stop tracer
    }                                                           // End overload operator()
  };

  //! Tree traversal of periodic cells
  void traversePeriodic(real_t cycle) {
    logger::startTimer("Traverse periodic");                    // Start timer
    vec3 Xperiodic = 0;                                         // Periodic coordinate offset
    Cells pcells; pcells.resize(27);                            // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->ICHILD = 0;                                             // Child cells for periodic center cell
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
		    kernel::M2L(Ci0, Ci, Xperiodic, false);     //         M2L kernel
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
#if MASS
      for (int i=1; i<NTERM; i++) Ci->M[i] *= Ci->M[0];         //  Normalize multipole expansion coefficients
#endif
      Cj0 = pcells.begin();                                     //  Redefine Cj0 for M2M
      C_iter Cj = Cj0;                                          //  Iterator of periodic neighbor cells
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            if (ix != 0 || iy != 0 || iz != 0) {                //     If periodic cell is not at center
              Cj->X[0] = Ci->X[0] + ix * cycle;                 //      Set new x coordinate for periodic image
              Cj->X[1] = Ci->X[1] + iy * cycle;                 //      Set new y cooridnate for periodic image
              Cj->X[2] = Ci->X[2] + iz * cycle;                 //      Set new z coordinate for periodic image
              Cj->M    = Ci->M;                                 //      Copy multipoles to new periodic image
              Cj++;                                             //      Increment periodic cell iterator
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Ci->M = 0;                                                //  Reset multipoles of periodic parent
      kernel::M2M(Ci,Cj0);                                      //  Evaluate periodic M2M kernels for this sublevel
#if MASS
      for (int i=1; i<NTERM; i++) Ci->M[i] /= Ci->M[0];         //  Normalize multipole expansion coefficients
#endif
      cycle *= 3;                                               //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
#if MASS
    Ci0->L /= Ci0->M[0];                                        // Normalize local expansion coefficients
#endif
    logger::stopTimer("Traverse periodic");                     // Stop timer
  }

  //! Split cell and call traverse() recursively for child
  void splitCell(C_iter Ci, C_iter Cj, vec3 Xperiodic, bool mutual, real_t remote) {
    if (Cj->NCHILD == 0) {                                      // If Cj is leaf
      assert(Ci->NCHILD > 0);                                   //  Make sure Ci is not leaf
      for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        traverse(ci, Cj, Xperiodic, mutual, remote);            //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else if (Ci->NCHILD == 0) {                               // Else if Ci is leaf
      assert(Cj->NCHILD > 0);                                   //  Make sure Cj is not leaf
      for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
        traverse(Ci, cj, Xperiodic, mutual, remote);            //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    } else if (Ci->NBODY + Cj->NBODY >= nspawn || (mutual && Ci == Cj)) {// Else if cells are still large
      TraverseRange traverseRange(this, Ci0+Ci->ICHILD, Ci0+Ci->ICHILD+Ci->NCHILD,// Instantiate recursive functor
				  Cj0+Cj->ICHILD, Cj0+Cj->ICHILD+Cj->NCHILD, Xperiodic, mutual, remote);
      traverseRange();                                          //  Traverse for range of cell pairs
    } else if (Ci->R >= Cj->R) {                                // Else if Ci is larger than Cj
      for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        traverse(ci, Cj, Xperiodic, mutual, remote);            //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Cj is larger than Ci
      for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
        traverse(Ci, cj, Xperiodic, mutual, remote);            //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

public:
  //! Constructor
  Traversal(int _nspawn, int _images) : nspawn(_nspawn), images(_images) // Initialize variables
#if COUNT
				      , numP2P(0), numM2L(0)
#endif
  {}

#if USE_WEIGHT
  //! Initialize interaction weights of bodies and cells
  void initWeight(Cells & cells) {
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {         // Loop over cells
      C->WEIGHT = 0;                                            //  Initialize cell weights
      if (C->NCHILD==0) {                                       //  If leaf cell
	for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {      //   Loop over bodies in cell
	  B->WEIGHT = 0;                                        //    Initialize body weights
	}                                                       //   End loop over bodies in cell
      }                                                         //  End if for leaf cell
    }                                                           // End loop over cells
  }
#else
  void initWeight(Cells) {}
#endif

  //! Evaluate P2P and M2L using dual tree traversal
  void dualTreeTraversal(Cells & icells, Cells & jcells, real_t cycle, bool mutual, real_t remote=1) {
    if (icells.empty() || jcells.empty()) return;               // Quit if either of the cell vectors are empty
    logger::startTimer("Traverse");                             // Start timer
    logger::initTracer();                                       // Initialize tracer
    Ci0 = icells.begin();                                       // Set iterator of target root cell
    Cj0 = jcells.begin();                                       // Set iterator of source root cell
    vec3 Xperiodic = 0;                                         // Periodic coordinate offset
    if (images == 0) {                                          // If non-periodic boundary condition
      traverse(Ci0, Cj0, Xperiodic, mutual, remote);            //  Traverse the tree
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
	for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
	  for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
	    Xperiodic[0] = ix * cycle;                          //     Coordinate shift for x periodic direction
	    Xperiodic[1] = iy * cycle;                          //     Coordinate shift for y periodic direction
	    Xperiodic[2] = iz * cycle;                          //     Coordinate shift for z periodic direction
	    traverse(Ci0, Cj0, Xperiodic, false, remote);       //     Traverse the tree for this periodic image
	  }                                                     //    End loop over z periodic direction
	}                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic(cycle);                                  //  Traverse tree for periodic images
    }                                                           // End if for periodic boundary condition
    logger::stopTimer("Traverse");                              // Stop timer
    logger::writeTracer();                                      // Write tracer to file
  }

  struct DirectRecursion {
    C_iter Ci;                                                  //!< Iterator of target cell
    C_iter Cj;                                                  //!< Iterator of source cell
    int prange;                                                 //!< Range of periodic images
    real_t cycle;                                               //!< Periodic cycle
    DirectRecursion(C_iter _Ci, C_iter _Cj, int _prange, real_t _cycle) :// Constructor
      Ci(_Ci), Cj(_Cj), prange(_prange), cycle(_cycle) {}       // Initialize variables
    void operator() () {                                        // Overload operator
      if (Ci->NBODY < 25) {                                     // If number of target bodies is less than threshold
	vec3 Xperiodic = 0;                                     //   Periodic coordinate offset
	for (int ix=-prange; ix<=prange; ix++) {                //   Loop over x periodic direction
	  for (int iy=-prange; iy<=prange; iy++) {              //    Loop over y periodic direction
	    for (int iz=-prange; iz<=prange; iz++) {            //     Loop over z periodic direction
	      Xperiodic[0] = ix * cycle;                        //      Coordinate shift for x periodic direction
	      Xperiodic[1] = iy * cycle;                        //      Coordinate shift for y periodic direction
	      Xperiodic[2] = iz * cycle;                        //      Coordinate shift for z periodic direction
	      kernel::P2P(Ci, Cj, Xperiodic, false);            //      Evaluate P2P kernel
	    }                                                   //     End loop over z periodic direction
	  }                                                     //    End loop over y periodic direction
	}                                                       //   End loop over x periodic direction
      } else {                                                  // If number of target bodies is more than threshold
        Cells cells; cells.resize(1);                           //  Initialize new cell vector
	C_iter Ci2 = cells.begin();                             //  New cell iterator for right branch
	Ci2->BODY = Ci->BODY + Ci->NBODY / 2;                   //  Set begin iterator to handle latter half
	Ci2->NBODY = Ci->NBODY - Ci->NBODY / 2;                 //  Set range to handle latter half
	Ci->NBODY = Ci->NBODY / 2;                              //  Set range to handle first half
	mk_task_group;                                          //  Initialize task group
        DirectRecursion leftBranch(Ci, Cj, prange, cycle);      //  Instantiate recursive functor
	create_taskc(leftBranch);                               //  Create new task for left branch
	DirectRecursion rightBranch(Ci2, Cj, prange, cycle);    //  Instantiate recursive functor
	rightBranch();                                          //  Use old task for right branch
	wait_tasks;                                             //  Synchronize task group
      }                                                         // End if for NBODY threshold
    }                                                           // End operator
  };

  //! Direct summation
  void direct(Bodies & ibodies, Bodies & jbodies, real_t cycle) {
    Cells cells; cells.resize(2);                               // Define a pair of cells to pass to P2P kernel
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;            // First cell is target, second cell is source
    Ci->BODY = ibodies.begin();                                 // Iterator of first target body
    Ci->NBODY = ibodies.size();                                 // Number of target bodies
    Cj->BODY = jbodies.begin();                                 // Iterator of first source body
    Cj->NBODY = jbodies.size();                                 // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<images; i++) {                              // Loop over periodic image sublevels
      prange += int(std::pow(3.,i));                            //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
#if 1
    DirectRecursion directRecursion(Ci, Cj, prange, cycle);     // Instantiate recursive functor
    directRecursion();                                          // Recursive call for direct summation
#else
    for (int ix=-prange; ix<=prange; ix++) {                    // Loop over x periodic direction
      for (int iy=-prange; iy<=prange; iy++) {                  //  Loop over y periodic direction
	for (int iz=-prange; iz<=prange; iz++) {                //   Loop over z periodic direction
	  Xperiodic[0] = ix * cycle;                            //    Coordinate shift for x periodic direction
	  Xperiodic[1] = iy * cycle;                            //    Coordinate shift for y periodic direction
	  Xperiodic[2] = iz * cycle;                            //    Coordinate shift for z periodic direction
	  kernel::P2P(Ci, Cj, Xperiodic, false);                //    Evaluate P2P kernel
	}                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
#endif
  }

  //! Normalize bodies after direct summation
  void normalize(Bodies & bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG /= B->SRC;                                         //  Normalize by target charge
    }                                                           // End loop over bodies
  }

  //! Print traversal statistics
  void printTraversalData() {
#if COUNT
    if (logger::verbose) {                                      // If verbose flag is true
      std::cout << "--- Traversal stats --------------" << std::endl// Print title
		<< std::setw(stringLength) << std::left         //  Set format
		<< "P2P calls"  << " : "                        //  Print title
		<< std::setprecision(0) << std::fixed           //  Set format
		<< numP2P << std::endl                          //  Print number of P2P calls
		<< std::setw(stringLength) << std::left         //  Set format
		<< "M2L calls"  << " : "                        //  Print title
		<< std::setprecision(0) << std::fixed           //  Set format
		<< numM2L << std::endl;                         //  Print number of M2L calls
    }                                                           // End if for verbose flag
#endif
  }
};
#endif
