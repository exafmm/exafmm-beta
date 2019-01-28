#ifndef traversal_h
#define traversal_h
#include "logger.h"
#include "namespace.h"
#include "thread.h"
#include "types.h"

#if EXAFMM_COUNT_KERNEL
#define countKernel(N) N++
#else
#define countKernel(N)
#endif

namespace EXAFMM_NAMESPACE {
  class Traversal {
  private:
    Kernel & kernel;                                            //!< Kernel class
    const real_t theta;                                         //!< Multipole acceptance criteria
    const int nspawn;                                           //!< Threshold of NBODY for spawning new threads
    const int images;                                           //!< Number of periodic image sublevels
    const char * path;                                          //!< Path to save files
#if EXAFMM_COUNT_KERNEL
    real_t numP2P;                                              //!< Number of P2P kernel calls
    real_t numM2L;                                              //!< Number of M2L kernel calls
#endif
    C_iter Ci0;                                                 //!< Iterator of first target cell
    C_iter Cj0;                                                 //!< Iterator of first source cell

  private:
#if EXAFMM_COUNT_LIST
    //! Accumulate interaction list size of cells
    void countList(C_iter Ci, C_iter Cj, bool isP2P) {
      if (isP2P) Ci->numP2P++;                                  // If P2P, increment P2P counter of target cell
      else Ci->numM2L++;                                        // Else, increment M2L counter of target cell
    }
#else
    void countList(C_iter, C_iter, bool) {}
#endif

#if EXAFMM_USE_WEIGHT
    //! Accumulate interaction weights of cells
    void countWeight(C_iter Ci, C_iter Cj, real_t weight) {
      Ci->WEIGHT += weight;                                     // Increment weight of target cell
    }
#else
    void countWeight(C_iter, C_iter, real_t) {}
#endif

    //! Get level from key
    int getLevel(uint64_t key) {
      int level = -1;                                           // Initialize level
      while( int(key) >= 0 ) {                                  // While key has level offsets to subtract
	level++;                                                //  Increment level
	key -= 1 << 3*level;                                    //  Subtract level offset
      }                                                         // End while loop for level offsets
      return level;                                             // Return level
    }

    //! Get 3-D index from key
    ivec3 getIndex(uint64_t key) {
      int level = -1;                                           // Initialize level
      while( int(key) >= 0 ) {                                  // While key has level offsets to subtract
	level++;                                                //  Increment level
	key -= 1 << 3*level;                                    //  Subtract level offset
      }                                                         // End while loop for level offsets
      key += 1 << 3*level;                                      // Compensate for over-subtraction
      level = 0;                                                // Initialize level
      ivec3 iX = 0;                                             // Initialize 3-D index
      int d = 0;                                                // Initialize dimension
      while( key > 0 ) {                                        // While key has bits to shift
	iX[d] += (key % 2) * (1 << level);                      //  Deinterleave key bits to 3-D bits
	key >>= 1;                                              //  Shift bits in key
	d = (d+1) % 3;                                          //  Increment dimension
	if( d == 0 ) level++;                                   //  Increment level
      }                                                         // End while loop for key bits to shift
      return iX;                                                // Return 3-D index
    }

    //! Get 3-D index from periodic key
    ivec3 getPeriodicIndex(int key) {
      ivec3 iX;                                                 // Initialize 3-D periodic index
      iX[0] = key % 3;                                          // x periodic index
      iX[1] = (key / 3) % 3;                                    // y periodic index
      iX[2] = key / 9;                                          // z periodic index
      iX -= 1;                                                  // {0,1,2} -> {-1,0,1}
      return iX;                                                // Return 3-D periodic index
    }

    //! Split cell and call traverse() recursively for child
    void splitCell(C_iter Ci, C_iter Cj, real_t remote) {
      if (Cj->NCHILD == 0) {                                    // If Cj is leaf
	assert(Ci->NCHILD > 0);                                 //  Make sure Ci is not leaf
	for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
	  dualTreeTraversal(ci, Cj, remote);                    //   Traverse a single pair of cells
	}                                                       //  End loop over Ci's children
      } else if (Ci->NCHILD == 0) {                             // Else if Ci is leaf
	assert(Cj->NCHILD > 0);                                 //  Make sure Cj is not leaf
	for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
	  dualTreeTraversal(Ci, cj, remote);                    //   Traverse a single pair of cells
	}                                                       //  End loop over Cj's children
      } else if (Ci->NBODY + Cj->NBODY >= nspawn || (Ci == Cj)) {// Else if cells are still large
	TraverseRange traverseRange(this, Ci0+Ci->ICHILD, Ci0+Ci->ICHILD+Ci->NCHILD,// Instantiate recursive functor
				    Cj0+Cj->ICHILD, Cj0+Cj->ICHILD+Cj->NCHILD, remote);
	traverseRange();                                        //  Traverse for range of cell pairs
      } else if (Ci->R >= Cj->R) {                              // Else if Ci is larger than Cj
	for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
	  dualTreeTraversal(ci, Cj, remote);                    //   Traverse a single pair of cells
	}                                                       //  End loop over Ci's children
      } else {                                                  // Else if Cj is larger than Ci
	for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
	  dualTreeTraversal(Ci, cj, remote);                    //   Traverse a single pair of cells
	}                                                       //  End loop over Cj's children
      }                                                         // End if for leafs and Ci Cj size
    }

    //! Dual tree traversal for a single pair of cells
    void dualTreeTraversal(C_iter Ci, C_iter Cj, real_t remote) {
      vec3 dX = Ci->X - Cj->X - kernel.Xperiodic;               // Distance vector from source to target
      real_t RT2 = norm(dX) * theta * theta;                    // Scalar distance squared
      if (RT2 > (Ci->R+Cj->R) * (Ci->R+Cj->R) * (1 - 1e-3)) {   // If distance is far enough
	kernel.M2L(Ci, Cj);                                     //  M2L kernel
	countKernel(numM2L);                                    //  Increment M2L counter
	countList(Ci, Cj, false);                               //  Increment M2L list
	countWeight(Ci, Cj, remote);                            //  Increment M2L weight
      } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {          // Else if both cells are bodies
#if EXAFMM_NO_P2P
	int index = Ci->ICELL;
	int iX[3] = {0, 0, 0};
	int d = 0, level = 0;
	while( index != 0 ) {
	  iX[d] += (index % 2) * (1 << level);
	  index >>= 1;
	  d = (d+1) % 3;
	  if( d == 0 ) level++;
	}
	index = Cj->ICELL;
	int jX[3] = {0, 0, 0};
	d = 0; level = 0;
	while( index != 0 ) {
	  jX[d] += (index % 2) * (1 << level);
	  index >>= 1;
	  d = (d+1) % 3;
	  if( d == 0 ) level++;
	}
	int isNeighbor = 1;
	for (d=0; d<3; d++) {
	  if (kernel.Xperiodic[d] > 1e-3) jX[d] += 5;
	  if (kernel.Xperiodic[d] < -1e-3) jX[d] -= 5;
	  isNeighbor &= abs(iX[d] - jX[d]) <= 1;
	}
#endif
	if (Cj->NBODY == 0) {                                   //  If the bodies weren't sent from remote node
	  //std::cout << "Warning: icell " << Ci->ICELL << " needs bodies from jcell" << Cj->ICELL << std::endl;
	  kernel.M2L(Ci, Cj);                                   //   M2L kernel
	  countKernel(numM2L);                                  //   Increment M2L counter
	  countList(Ci, Cj, false);                             //   Increment M2L list
	  countWeight(Ci, Cj, remote);                          //   Increment M2L weight
#if EXAFMM_NO_P2P
	} else if (!isNeighbor) {                               //  If GROAMCS handles neighbors
	  kernel.M2L(Ci, Cj);                                   //   M2L kernel
	  countKernel(numM2L);                                  //   Increment M2L counter
	  countList(Ci, Cj, false);                             //   Increment M2L list
	  countWeight(Ci, Cj, remote);                          //   Increment M2L weight
	} else {
	  countList(Ci, Cj, true);                              //   Increment P2P list
#else
	} else {
          kernel.P2P(Ci, Cj);                                   //   P2P kernel for pair of cells
	  countKernel(numP2P);                                  //   Increment P2P counter
	  countList(Ci, Cj, true);                              //   Increment P2P list
	  countWeight(Ci, Cj, remote);                          //   Increment P2P weight
#endif
	}                                                       //  End if for bodies
      } else {                                                  // Else if cells are close but not bodies
	splitCell(Ci, Cj, remote);                              //  Split cell and call function recursively for child
      }                                                         // End if for multipole acceptance
    }

    //! Recursive functor for dual tree traversal of a range of Ci and Cj
    struct TraverseRange {
      Traversal * traversal;                                    //!< Traversal object
      C_iter CiBegin;                                           //!< Begin iterator of target cells
      C_iter CiEnd;                                             //!< End iterator of target cells
      C_iter CjBegin;                                           //!< Begin Iterator of source cells
      C_iter CjEnd;                                             //!< End iterator of source cells
      real_t remote;                                            //!< Weight for remote work load
      TraverseRange(Traversal * _traversal, C_iter _CiBegin, C_iter _CiEnd,// Constructor
		    C_iter _CjBegin, C_iter _CjEnd, real_t _remote) :
	traversal(_traversal), CiBegin(_CiBegin), CiEnd(_CiEnd),// Initialize variables
	CjBegin(_CjBegin), CjEnd(_CjEnd), remote(_remote) {}
      void operator() () const {                                // Overload operator()
	Tracer tracer;                                          //  Instantiate tracer
	logger::startTracer(tracer);                            //  Start tracer
	if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {     //  If only one cell in range
	  if (CiBegin == CjBegin) {                             //   If Ci == Cj
	    assert(CiEnd == CjEnd);                             //    Check if self interaction
	    traversal->dualTreeTraversal(CiBegin, CjBegin, remote);//   Call traverse for single pair
	  } else {                                              //   If Ci != Cj
	    for (C_iter Ci=CiBegin; Ci!=CiEnd; Ci++) {          //    Loop over all Ci cells
	      for (C_iter Cj=CjBegin; Cj!=CjEnd; Cj++) {        //     Loop over all Cj cells
		traversal->dualTreeTraversal(Ci, Cj, remote);   //      Call traverse for single pair
	      }                                                 //     End loop over all Cj cells
	    }                                                   //    End loop over all Ci cells
	  }                                                     //   End if for Ci == Cj
	} else {                                                //  If many cells are in the range
	  C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;       //   Split range of Ci cells in half
	  C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;       //   Split range of Cj cells in half
	  mk_task_group;                                        //   Initialize task group
	  {
	    TraverseRange leftBranch(traversal, CiBegin, CiMid, //    Instantiate recursive functor
				     CjBegin, CjMid, remote);
	    create_taskc(leftBranch);                           //    Ci:former Cj:former
	    TraverseRange rightBranch(traversal, CiMid, CiEnd,  //    Instantiate recursive functor
				      CjMid, CjEnd, remote);
	    rightBranch();                                      //    Ci:latter Cj:latter
	    wait_tasks;                                         //    Synchronize task group
	  }
	  {
	    TraverseRange leftBranch(traversal, CiBegin, CiMid, //    Instantiate recursive functor
				     CjMid, CjEnd, remote);
	    create_taskc(leftBranch);                           //    Ci:former Cj:latter
            TraverseRange rightBranch(traversal, CiMid, CiEnd,  //    Instantiate recursive functor
                                      CjBegin, CjMid, remote);
            rightBranch();                                      //    Ci:latter Cj:former
	    wait_tasks;                                         //    Synchronize task group
	  }
	}                                                       //  End if for many cells in range
	logger::stopTracer(tracer);                             //  Stop tracer
      }                                                         // End overload operator()
    };

    //! Tree traversal of periodic cells
    void traversePeriodic(vec3 cycle) {
      logger::startTimer("Traverse periodic");                  // Start timer
      int prange = .5 / theta + 1;                              // Periodic range
      int neighbor = 2 * prange + 1;                            // Neighbor region size
      int listSize = neighbor * neighbor * neighbor;            // Neighbor list size
      Cells pcells; pcells.resize(listSize);                    // Create cells
      for (C_iter C=pcells.begin(); C!=pcells.end(); C++) {     // Loop over periodic cells
        C->M.resize(kernel.NTERM, 0.0);                         //  Allocate & initialize M coefs
        C->L.resize(kernel.NTERM, 0.0);                         //  Allocate & initialize L coefs
      }                                                         // End loop over periodic cells
      C_iter Ci = pcells.end()-1;                               // Last cell is periodic parent cell
      *Ci = *Cj0;                                               // Copy values from source root
      Ci->ICHILD = 0;                                           // Child cells for periodic center cell
      Ci->NCHILD = listSize - 1;                                // Number of child cells for periodic center cell
      C_iter C0 = Cj0;                                          // Placeholder for Cj0
      for (int level=0; level<images-1; level++) {              // Loop over sublevels of tree
	for (int ix=-prange; ix<=prange; ix++) {                //  Loop over x periodic direction
	  for (int iy=-prange; iy<=prange; iy++) {              //   Loop over y periodic direction
	    for (int iz=-prange; iz<=prange; iz++) {            //    Loop over z periodic direction
	      if (ix != 0 || iy != 0 || iz != 0) {              //     If periodic cell is not at center
		for (int cx=-prange; cx<=prange; cx++) {        //      Loop over x periodic direction (child)
		  for (int cy=-prange; cy<=prange; cy++) {      //       Loop over y periodic direction (child)
		    for (int cz=-prange; cz<=prange; cz++) {    //        Loop over z periodic direction (child)
		      kernel.Xperiodic[0] = (ix * neighbor + cx) * cycle[0];//   Coordinate offset for x periodic direction
		      kernel.Xperiodic[1] = (iy * neighbor + cy) * cycle[1];//   Coordinate offset for y periodic direction
		      kernel.Xperiodic[2] = (iz * neighbor + cz) * cycle[2];//   Coordinate offset for z periodic direction
		      kernel.M2L(Ci0, Ci);                      //         M2L kernel
		    }                                           //        End loop over z periodic direction (child)
		  }                                             //       End loop over y periodic direction (child)
		}                                               //      End loop over x periodic direction (child)
	      }                                                 //     Endif for periodic center cell
	    }                                                   //    End loop over z periodic direction
	  }                                                     //   End loop over y periodic direction
	}                                                       //  End loop over x periodic direction
	Cj0 = pcells.begin();                                   //  Redefine Cj0 for M2M
	C_iter Cj = Cj0;                                        //  Iterator of periodic neighbor cells
	for (int ix=-prange; ix<=prange; ix++) {                //  Loop over x periodic direction
	  for (int iy=-prange; iy<=prange; iy++) {              //   Loop over y periodic direction
	    for (int iz=-prange; iz<=prange; iz++) {            //    Loop over z periodic direction
	      if (ix != 0 || iy != 0 || iz != 0) {              //     If periodic cell is not at center
		Cj->X[0] = Ci->X[0] + ix * cycle[0];            //      Set new x coordinate for periodic image
		Cj->X[1] = Ci->X[1] + iy * cycle[1];            //      Set new y cooridnate for periodic image
		Cj->X[2] = Ci->X[2] + iz * cycle[2];            //      Set new z coordinate for periodic image
		Cj->M = Ci->M;                                  //      Copy multipoles to new periodic image
		Cj++;                                           //      Increment periodic cell iterator
	      }                                                 //     Endif for periodic center cell
	    }                                                   //    End loop over z periodic direction
	  }                                                     //   End loop over y periodic direction
	}                                                       //  End loop over x periodic direction
	kernel.M2M(Ci,Cj0);                                     //  Evaluate periodic M2M kernels for this sublevel
	cycle *= neighbor;                                      //  Increase periodic cycle by number of neighbors
	Cj0 = C0;                                               //  Reset Cj0 back
      }                                                         // End loop over sublevels of tree
      logger::stopTimer("Traverse periodic");                   // Stop timer
    }

  public:
    //! Constructor
    Traversal(Kernel & _kernel, real_t _theta, int _nspawn, int _images, const char * _path) : // Constructor
      kernel(_kernel), theta(_theta), nspawn(_nspawn), images(_images), path(_path) // Initialize variables
#if EXAFMM_COUNT_KERNEL
      , numP2P(0), numM2L(0)
#endif
    {}

#if EXAFMM_COUNT_LIST
    //! Initialize size of P2P and M2L interaction lists per cell
    void initListCount(Cells & cells) {
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	C->numP2P = C->numM2L = 0;                              //  Initialize size of interaction list
      }                                                         // End loop over cells
    }
#else
    void initListCount(Cells) {}
#endif

#if EXAFMM_USE_WEIGHT
    //! Initialize interaction weights of bodies and cells
    void initWeight(Cells & cells) {
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over cells
	C->WEIGHT = 0;                                          //  Initialize cell weights
	if (C->NCHILD==0) {                                     //  If leaf cell
	  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {    //   Loop over bodies in cell
	    B->WEIGHT = 0;                                      //    Initialize body weights
	  }                                                     //   End loop over bodies in cell
	}                                                       //  End if for leaf cell
      }                                                         // End loop over cells
    }
#else
    void initWeight(Cells) {}
#endif

    //! Evaluate P2P and M2L using list based traversal
    void traverse(Cells & icells, Cells & jcells, vec3 cycle, bool dual, real_t remote=1) {
      if (icells.empty() || jcells.empty()) return;             // Quit if either of the cell vectors are empty
      logger::startTimer("Traverse");                           // Start timer
      logger::initTracer();                                     // Initialize tracer
      int prange = .5 / theta + 1;                              // Periodic range
      prange = 1;
      Ci0 = icells.begin();                                     // Iterator of first target cell
      Cj0 = jcells.begin();                                     // Iterator of first source cell
      kernel.Xperiodic = 0;                                     // Set periodic coordinate offset to 0
      if (images == 0) {                                        //  If non-periodic boundary condition
        dualTreeTraversal(Ci0, Cj0, remote);                    //   Traverse the tree
      } else {                                                  //  If periodic boundary condition
        for (int ix=-prange; ix<=prange; ix++) {                //   Loop over x periodic direction
          for (int iy=-prange; iy<=prange; iy++) {              //    Loop over y periodic direction
            for (int iz=-prange; iz<=prange; iz++) {            //     Loop over z periodic direction
              kernel.Xperiodic[0] = ix * cycle[0];              //      Coordinate shift for x periodic direction
              kernel.Xperiodic[1] = iy * cycle[1];              //      Coordinate shift for y periodic direction
              kernel.Xperiodic[2] = iz * cycle[2];              //      Coordinate shift for z periodic direction
              dualTreeTraversal(Ci0, Cj0, remote);              //      Traverse the tree for this periodic image
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
        traversePeriodic(cycle);                                //   Traverse tree for periodic images
      }                                                         //  End if for periodic boundary condition
      logger::stopTimer("Traverse");                            // Stop timer
      logger::writeTracer();                                    // Write tracer to file
    }

    //! Direct summation
    void direct(Bodies & ibodies, Bodies & jbodies, vec3 cycle) {
      Cells cells; cells.resize(2);                             // Define a pair of cells to pass to P2P kernel
      C_iter Ci = cells.begin(), Cj = cells.begin()+1;          // First cell is target, second cell is source
      int neighbor = 2 * (.5 / theta + 1) + 1;                  // Neighbor reigon
      int prange = 0;                                           // Range of periodic images
      for (int i=0; i<images; i++) {                            // Loop over periodic image sublevels
	prange += int(pow(neighbor,i));                         //  Accumulate range of periodic images
      }                                                         // End loop over perioidc image sublevels
      for (int ix=-prange; ix<=prange; ix++) {                  //  Loop over x periodic direction
	for (int iy=-prange; iy<=prange; iy++) {                //   Loop over y periodic direction
	  for (int iz=-prange; iz<=prange; iz++) {              //    Loop over z periodic direction
	    kernel.Xperiodic[0] = ix * cycle[0];                //     Coordinate shift for x periodic direction
	    kernel.Xperiodic[1] = iy * cycle[1];                //     Coordinate shift for y periodic direction
	    kernel.Xperiodic[2] = iz * cycle[2];                //     Coordinate shift for z periodic direction
	    Ci->BODY = ibodies.begin();                         //     Iterator of first target body
	    Ci->NBODY = ibodies.size();                         //     Number of target bodies
	    Cj->BODY = jbodies.begin();                         //     Iterator of first source body
	    Cj->NBODY = jbodies.size();                         //     Number of source bodies
            kernel.P2P(Ci, Cj);                                 //     Evaluate P2P kenrel
	  }                                                     //    End loop over z periodic direction
	}                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
    }

    //! Write G matrix to file
    void writeMatrix(Bodies & bodies, Bodies & jbodies) {
      std::stringstream name;                                   // File name
      name << path << "matrix.dat";                             // Create file name for matrix
      std::ofstream matrixFile(name.str().c_str());             // Open matrix log file
      for (B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {  // Loop over target bodies
	for (B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {//  Loop over source bodies
	  vec3 dX = Bi->X - Bj->X;                              //   Distance vector
	  real_t R2 = norm(dX) + kernel.eps2;                   //   Distance squared
	  real_t G = R2 == 0 ? 0.0 : 1.0 / sqrt(R2);            //   Green's function
	  matrixFile << G << " ";                               //   Write to matrix data file
	}                                                       //  End loop over source bodies
	matrixFile << std::endl;                                //  Line break
      }                                                         // End loop over target bodies
    }

    //! Print traversal statistics
    void printTraversalData() {
#if EXAFMM_COUNT_KERNEL
      if (logger::verbose) {                                    // If verbose flag is true
	std::cout << "--- Traversal stats --------------" << std::endl// Print title
		  << std::setw(logger::stringLength) << std::left //  Set format
		  << "P2P calls"  << " : "                      //  Print title
		  << std::setprecision(0) << std::fixed         //  Set format
		  << numP2P << std::endl                        //  Print number of P2P calls
		  << std::setw(logger::stringLength) << std::left //  Set format
		  << "M2L calls"  << " : "                      //  Print title
		  << std::setprecision(0) << std::fixed         //  Set format
		  << numM2L << std::endl;                       //  Print number of M2L calls
      }                                                         // End if for verbose flag
#endif
    }
#if EXAFMM_COUNT_LIST
    void writeList(Cells cells, int mpirank) {
      std::stringstream name;                                   // File name
      name << path << "list" << std::setfill('0') << std::setw(6) // Set format
	   << mpirank << ".dat";                                // Create file name for list
      std::ofstream listFile(name.str().c_str());               // Open list log file
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {       // Loop over all lists
	listFile << std::setw(logger::stringLength) << std::left//  Set format
		 << C->ICELL << " " << C->numP2P << " " << C->numM2L << std::endl; // Print list size
      }                                                         // End loop over all lists
    }
#else
    void writeList(Cells, int) {}
#endif
  };
}
#endif
