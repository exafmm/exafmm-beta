#ifndef traversal_h
#define traversal_h
#include "kernel.h"
#include "logger.h"
#include "thread.h"

#if COUNT_KERNEL
#define countKernel(N) N++
#else
#define countKernel(N)
#endif

class Traversal {
private:
  const int nspawn;                                             //!< Threshold of NBODY for spawning new threads
  const int images;                                             //!< Number of periodic image sublevels
  int (* listOffset)[3];                                        //!< Offset in interaction lists
  int (* lists)[2];                                             //!< Interaction lists
#if COUNT_KERNEL
  real_t numP2P;                                                //!< Number of P2P kernel calls
  real_t numM2L;                                                //!< Number of M2L kernel calls
#endif
  C_iter Ci0;                                                   //!< Iterator of first target cell
  C_iter Cj0;                                                   //!< Iterator of first source cell

private:
#if COUNT_LIST
  //! Accumulate interaction list size of cells
  void countList(C_iter Ci, C_iter Cj, bool mutual, bool isP2P) {
    if (isP2P) Ci->numP2P++;                                    // If P2P, increment P2P counter of target cell
    else Ci->numM2L++;                                          // Else, increment M2L counter of target cell
    if (mutual) {                                               // If mutual interaction in on
      if (isP2P) Cj->numP2P++;                                  //  If P2P, increment P2P counter of source cell
      else Cj->numM2L++;                                        //  Else, increment M2L counter of source cell
    }                                                           // End if for mutual interaction
  }
#else
  void countList(C_iter, C_iter, bool, bool) {}
#endif

#if USE_WEIGHT
  //! Accumulate interaction weights of cells
  void countWeight(C_iter Ci, C_iter Cj, bool mutual, real_t weight) {
    Ci->WEIGHT += weight;                                       // Increment weight of target cell
    if (mutual) Cj->WEIGHT += weight;                           // Increment weight of source cell
  }
#else
  void countWeight(C_iter, C_iter, bool, real_t) {}
#endif

  //! Get 3-D index from key
  ivec3 getIndex(uint64_t key) {
    int level = -1;                                             // Initialize level
    while( int(key) >= 0 ) {                                    // While key has level offsets to subtract
      level++;                                                  //  Increment level
      key -= 1 << 3*level;                                      //  Subtract level offset
    }                                                           // End while loop for level offsets
    key += 1 << 3*level;                                        // Compensate for over-subtraction
    level = 0;                                                  // Initialize level
    ivec3 iX = 0;                                               // Initialize 3-D index
    int d = 0;                                                  // Initialize dimension
    while( key > 0 ) {                                          // While key has bits to shift
      iX[d] += (key % 2) * (1 << level);                        //  Deinterleave key bits to 3-D bits
      key >>= 1;                                                //  Shift bits in key
      d = (d+1) % 3;                                            //  Increment dimension
      if( d == 0 ) level++;                                     //  Increment level
    }                                                           // End while loop for key bits to shift
    return iX;                                                  // Return 3-D index
  }

  //! Reset cell radius
  void resetCellRadius(C_iter C, C_iter C0, real_t R0, int level) {
    C->R = R0 / (1 << level);                                   // Set cell radius
    for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {// Loop over child cells
      resetCellRadius(CC, C0, R0, level+1);                     //  Recursive call for child cells
    }                                                           // End loop over child cells
  }

  //! Get interaction list
  void getList(int itype, int icell, int * list, int & numList) {
    int ilast = listOffset[icell][itype];                       // Initialize list pointer
    numList = 0;                                                // Initialize list size
    while (ilast >= 0) {                                        // While pointer exists
      if (lists[ilast][1] > 0) {                                //  If pointer is valid
	list[numList] = lists[ilast][1];                        //   Store interaction list in list
	numList++;                                              //   Increment list size
      }                                                         //  End if for valid pointer
      ilast = lists[ilast][0];                                  //  Increment pointer
    }                                                           // End while loop for pointer existence
  }

  //! Set one interaction list
  void setList(int itype, int icell, int list, int & numLists) {
    lists[numLists][0] = listOffset[icell][itype];              // Store list pointer
    lists[numLists][1] = list;                                  // Store list
    listOffset[icell][itype] = numLists;                        // Store list size
    numLists++;                                                 // Increment list size
  }

  //! Set all interaction lists
  void setLists(Cells & cells) {
    int numCells = cells.size();                                // Number of cells
    int childs[216], neighbors[27];                             // Array of parents' neighbors' children and neighbors
    C_iter C0 = cells.begin();                                  // Iterator of first cell
    for (int i=0; i<numCells; i++) {                            // Loop over number of cells
      for (int j=0; j<3; j++) {                                 //  Loop over list types
	listOffset[i][j] = -1;                                  //   Set initial value to -1
      }                                                         //  End loop over list types
    }                                                           // End loop over number of cells
    int numLists = 0;                                           // Initialize number of lists
    for (int icell=1; icell<numCells; icell++) {                // Loop over target cells
      C_iter Ci = C0 + icell;                                   //  Iterator of current target cell
      int iparent = Ci->IPARENT;                                //  Index of parent target cell
      neighbors[0] = iparent;                                   //  Include self to neighbor list
      int numNeighbors;                                         //  Number of neighbor parents
      getList(2, iparent, &neighbors[1], numNeighbors);         //  Get list of parents' neighbors
      numNeighbors++;                                           //  Increment number of parents' neighbors
      ivec3 iX = getIndex(Ci->ICELL);                           //  Get 3-D index from key
      int nchilds = 0;                                          //  Initialize number of parents' neighbors' children
      for (int i=0; i<numNeighbors; i++) {                      //  Loop over parents' neighbors
	int jparent = neighbors[i];                             //   Index of parent source cell
	C_iter Cj = C0 + jparent;                               //   Iterator of parent source cell
	for (int j=0; j<Cj->NCHILD; j++) {                      //   Loop over children of parents' neighbors
	  int jcell = Cj->ICHILD+j;                             //    Index of source cell
	  if (jcell != icell) {                                 //    If target != source
	    childs[nchilds] = jcell;                            //     Store index of source cell
	    nchilds++;                                          //     Increment number of parents' neighbors' children
	  }                                                     //    End if for target != source
	}                                                       //   End loop over children of parents' neighbors
      }                                                         //  End loop over parents' neighbors
      for (int i=0; i<nchilds; i++) {                           //  Loop over children of parents' neighbors
	int jcell = childs[i];                                  //   Index of source cell
	C_iter Cj = C0 + jcell;                                 //   Iterator of source cell
	ivec3 jX = getIndex(Cj->ICELL);                         //   3-D index of source cell
	if (iX[0]-1 <= jX[0] && jX[0] <= iX[0]+1 &&             //   If neighbor in x dimension and
	    iX[1]-1 <= jX[1] && jX[1] <= iX[1]+1 &&             //               in y dimension and
	    iX[2]-1 <= jX[2] && jX[2] <= iX[2]+1) {             //               in z dimension
	  setList(2, icell, jcell, numLists);                   //    Store neighbor list (not P2P unless leaf)
	} else {                                                //   If non-neighbor
	  setList(1, icell, jcell, numLists);                   //    Store M2L list
	}                                                       //   End if for non-neighbor
      }                                                         //  End loop over children of parents' neighbors
    }                                                           // End loop over target cells
    for (int icell=0; icell<numCells; icell++) {                // Loop over target cells
      C_iter Ci = C0 + icell;                                   //  Iterator of target cell
      if (Ci->NCHILD == 0) {                                    //  If taget cell is leaf
	int numNeighbors;                                       //   Number of neighbors
	getList(2, icell, neighbors, numNeighbors);             //   Get list of neighbor cells
	for (int j=0; j<numNeighbors; j++) {                    //   Loop over neighbor cells
	  int jcell = neighbors[j];                             //    Index of source cell
	  C_iter Cj = C0 + jcell;                               //    Iterator of source cell
	  if (Cj->NCHILD == 0) {                                //    If source cell is leaf
	    setList(0, icell, jcell, numLists);                 //     Store P2P list
	  }                                                     //    End if for source cell leaf
	}                                                       //   End loop over neighbor cells
      }                                                         //  End if for target cell leaf
    }                                                           // End loop over target cells
  }

  //! Dual tree traversal for a single pair of cells
  void traverse(C_iter Ci, C_iter Cj, bool mutual, real_t remote) {
    vec3 dX = Ci->X - Cj->X;                                    // Distance vector from source to target
    real_t R2 = norm(dX);                                       // Scalar distance squared
    if (R2 > (Ci->R+Cj->R) * (Ci->R+Cj->R) * (1 - 1e-3)) {      // If distance is far enough
      kernel::M2L(Ci, Cj, mutual);                              //  M2L kernel
      countKernel(numM2L);                                      //  Increment M2L counter
      countList(Ci, Cj, mutual, false);                         //  Increment M2L list
      countWeight(Ci, Cj, mutual, remote);                      //  Increment M2L weight
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are bodies
#if NO_P2P
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
	if (Xperiodic[d] > 1e-3) jX[d] += 5;
	if (Xperiodic[d] < -1e-3) jX[d] -= 5;
	isNeighbor &= abs(iX[d] - jX[d]) <= 1;
      }
#endif
      if (Cj->NBODY == 0) {                                     //  If the bodies weren't sent from remote node
	//std::cout << "Warning: icell " << Ci->ICELL << " needs bodies from jcell" << Cj->ICELL << std::endl;
	kernel::M2L(Ci, Cj, mutual);                            //   M2L kernel
        countKernel(numM2L);                                    //   Increment M2L counter
	countList(Ci, Cj, mutual, false);                       //   Increment M2L list
	countWeight(Ci, Cj, mutual, remote);                    //   Increment M2L weight
#if NO_P2P
      } else if (!isNeighbor) {                                 //  If GROAMCS handles neighbors
	kernel::M2L(Ci, Cj, mutual);                            //   M2L kernel
        countKernel(numM2L);                                    //   Increment M2L counter
        countList(Ci, Cj, mutual, false);                       //   Increment M2L list
        countWeight(Ci, Cj, mutual, remote);                    //   Increment M2L weight
      } else {
	countList(Ci, Cj, mutual, true);                        //   Increment P2P list
#else
      } else {
	if (R2 == 0 && Ci == Cj) {                              //   If source and target are same
	  kernel::P2P(Ci);                                      //    P2P kernel for single cell
	} else {                                                //   Else if source and target are different
	  kernel::P2P(Ci, Cj, mutual);                          //    P2P kernel for pair of cells
	}                                                       //   End if for same source and target
	countKernel(numP2P);                                    //   Increment P2P counter
	countList(Ci, Cj, mutual, true);                        //   Increment P2P list
	countWeight(Ci, Cj, mutual, remote);                    //   Increment P2P weight
#endif
      }                                                         //  End if for bodies
    } else {                                                    // Else if cells are close but not bodies
      splitCell(Ci, Cj, mutual, remote);                        //  Split cell and call function recursively for child
    }                                                           // End if for multipole acceptance
  }

  //! Recursive functor for dual tree traversal of a range of Ci and Cj
  struct TraverseRange {
    Traversal * traversal;                                      //!< Traversal object
    C_iter CiBegin;                                             //!< Begin iterator of target cells
    C_iter CiEnd;                                               //!< End iterator of target cells
    C_iter CjBegin;                                             //!< Begin Iterator of source cells
    C_iter CjEnd;                                               //!< End iterator of source cells
    bool mutual;                                                //!< Flag for mutual interaction
    real_t remote;                                              //!< Weight for remote work load
    TraverseRange(Traversal * _traversal, C_iter _CiBegin, C_iter _CiEnd,// Constructor
		  C_iter _CjBegin, C_iter _CjEnd,
		  bool _mutual, real_t _remote) :
      traversal(_traversal), CiBegin(_CiBegin), CiEnd(_CiEnd),  // Initialize variables
      CjBegin(_CjBegin), CjEnd(_CjEnd), mutual(_mutual), remote(_remote) {}
    void operator() () {                                        // Overload operator()
      Tracer tracer;                                            //  Instantiate tracer
      logger::startTracer(tracer);                              //  Start tracer
      if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {       //  If only one cell in range
	if (CiBegin == CjBegin) {                               //   If Ci == Cj
	  assert(CiEnd == CjEnd);                               //    Check if mutual & self interaction
	  traversal->traverse(CiBegin, CjBegin, mutual, remote);//   Call traverse for single pair
	} else {                                                //   If Ci != Cj
	  for (C_iter Ci=CiBegin; Ci!=CiEnd; Ci++) {            //    Loop over all Ci cells
	    for (C_iter Cj=CjBegin; Cj!=CjEnd; Cj++) {          //     Loop over all Cj cells
	      traversal->traverse(Ci, Cj, mutual, remote);      //      Call traverse for single pair
	    }                                                   //     End loop over all Cj cells
	  }                                                     //    End loop over all Ci cells
	}                                                       //   End if for Ci == Cj
      } else {                                                  //  If many cells are in the range
	C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;         //   Split range of Ci cells in half
	C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;         //   Split range of Cj cells in half
	mk_task_group;                                          //   Initialize task group
	{
	  TraverseRange leftBranch(traversal, CiBegin, CiMid,   //    Instantiate recursive functor
				   CjBegin, CjMid, mutual, remote);
	  create_taskc(leftBranch);                             //    Ci:former Cj:former
	  TraverseRange rightBranch(traversal, CiMid, CiEnd,    //    Instantiate recursive functor
				    CjMid, CjEnd, mutual, remote);
	  rightBranch();                                        //    Ci:latter Cj:latter
	  wait_tasks;                                           //    Synchronize task group
	}
	{
	  TraverseRange leftBranch(traversal, CiBegin, CiMid,   //    Instantiate recursive functor
				   CjMid, CjEnd, mutual, remote);
	  create_taskc(leftBranch);                             //    Ci:former Cj:latter
	  if (!mutual || CiBegin != CjBegin) {                  //    Exclude mutual & self interaction
            TraverseRange rightBranch(traversal, CiMid, CiEnd,  //    Instantiate recursive functor
				      CjBegin, CjMid, mutual, remote);
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
		    kernel::Xperiodic[0] = (ix * 3 + cx) * cycle;//        Coordinate offset for x periodic direction
		    kernel::Xperiodic[1] = (iy * 3 + cy) * cycle;//        Coordinate offset for y periodic direction
		    kernel::Xperiodic[2] = (iz * 3 + cz) * cycle;//        Coordinate offset for z periodic direction
		    kernel::M2L(Ci0, Ci, false);                //         M2L kernel
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
  void splitCell(C_iter Ci, C_iter Cj, bool mutual, real_t remote) {
    if (Cj->NCHILD == 0) {                                      // If Cj is leaf
      assert(Ci->NCHILD > 0);                                   //  Make sure Ci is not leaf
      for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        traverse(ci, Cj, mutual, remote);                       //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else if (Ci->NCHILD == 0) {                               // Else if Ci is leaf
      assert(Cj->NCHILD > 0);                                   //  Make sure Cj is not leaf
      for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
        traverse(Ci, cj, mutual, remote);                       //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    } else if (Ci->NBODY + Cj->NBODY >= nspawn || (mutual && Ci == Cj)) {// Else if cells are still large
      TraverseRange traverseRange(this, Ci0+Ci->ICHILD, Ci0+Ci->ICHILD+Ci->NCHILD,// Instantiate recursive functor
				  Cj0+Cj->ICHILD, Cj0+Cj->ICHILD+Cj->NCHILD, mutual, remote);
      traverseRange();                                          //  Traverse for range of cell pairs
    } else if (Ci->R >= Cj->R) {                                // Else if Ci is larger than Cj
      for (C_iter ci=Ci0+Ci->ICHILD; ci!=Ci0+Ci->ICHILD+Ci->NCHILD; ci++) {// Loop over Ci's children
        traverse(ci, Cj, mutual, remote);                       //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Cj is larger than Ci
      for (C_iter cj=Cj0+Cj->ICHILD; cj!=Cj0+Cj->ICHILD+Cj->NCHILD; cj++) {// Loop over Cj's children
        traverse(Ci, cj, mutual, remote);                       //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

public:
  //! Constructor
  Traversal(int _nspawn, int _images) :                         // Constructor
    nspawn(_nspawn), images(_images)                            // Initialize variables
#if COUNT_KERNEL
    , numP2P(0), numM2L(0)
#endif
  {}

#if COUNT_LIST
  //! Initialize size of P2P and M2L interaction lists per cell
  void initListCount(Cells & cells) {
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {         // Loop over cells
      C->numP2P = C->numM2L = 0;                                //  Initialize size of interaction list
    }                                                           // End loop over cells
  }
#else
  void initListCount(Cells) {}
#endif

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

  //! Evaluate P2P and M2L using list based traversal
  void listBasedTraversal(Cells & cells, real_t remote=1) {
    int numCells = cells.size();                                // Number of cells
    C_iter C0 = cells.begin();                                  // Iterator of first target cell
    real_t R0 = C0->R;                                          // Radius of root cell
    kernel::Xperiodic = 0;                                      // Set periodic coordinate offset to 0
    bool mutual = false;                                        // Set mutual interaction flag to false
    int list[189];                                              // Current interaction list
    listOffset = new int [numCells][3]();                       // Offset of interaction lists
    lists = new int [189*numCells][2]();                        // All interaction lists
    resetCellRadius(C0, C0, R0, 0);                             // Reset cell radius
    setLists(cells);                                            // Set P2P and M2L interaction lists

    logger::startTimer("M2L");                                  // Start timer
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=0; icell<numCells; icell++) {                // Loop over target cells
      C_iter Ci = C0 + icell;                                   //  Iterator of target cell
      int nlist;                                                //  Interaction list size
      getList(1, icell, list, nlist);                           //   Get M2L interaction list
      for (int ilist=0; ilist<nlist; ilist++) {                 //   Loop over M2L interaction list
	int jcell = list[ilist];                                //    Index of source cell
	C_iter Cj = C0 + jcell;                                 //    Iterator of source cell
	kernel::M2L(Ci, Cj, mutual);                            //    M2L kernel
        countKernel(numM2L);                                    //  Increment M2L counter
        countList(Ci, Cj, mutual, false);                       //  Increment M2L list
        countWeight(Ci, Cj, mutual, remote);                    //  Increment M2L weight
      }                                                         //   End loop over M2L interaction list
    }                                                           // End loop over target cells
    logger::stopTimer("M2L");                                   // Stop timer

#ifndef NO_P2P
    logger::startTimer("P2P");                                  // Start timer
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=0; icell<numCells; icell++) {                // Loop over target cells
      C_iter Ci = C0 + icell;                                   //  Iterator of target cell
      if (Ci->NCHILD == 0) {                                    //  If target cell is leaf
	kernel::P2P(Ci, Ci, mutual);                            //   P2P kernel for self
	countKernel(numP2P);                                    //   Increment P2P counter
	countList(Ci, Ci, mutual, true);                        //   Increment P2P list
	countWeight(Ci, Ci, mutual, remote);                    //   Increment P2P weight
	int nlist;                                              //   Interaction list size
	getList(0, icell, list, nlist);                         //   Get P2P interaction list
	for (int ilist=0; ilist<nlist; ilist++) {               //   Loop over P2P interaction list
	  int jcell = list[ilist];                              //    Index of source cell
	  C_iter Cj = C0 + jcell;                               //    Iterator of source cell
	  kernel::P2P(Ci, Cj, mutual);                          //    P2P kernel
  	  countKernel(numP2P);                                  //   Increment P2P counter
	  countList(Ci, Cj, mutual, true);                      //   Increment P2P list
	  countWeight(Ci, Cj, mutual, remote);                  //   Increment P2P weight
	}                                                       //   End loop over P2P interaction list
      }                                                         //  End if for target cell leaf
    }                                                           // End loop over target cells
    logger::stopTimer("P2P");                                   // Stop timer
#endif
    delete[] listOffset;                                        // Deallocate offset of lists
    delete[] lists;                                             // Deallocate lists
  }

  //! Evaluate P2P and M2L using dual tree traversal
  void dualTreeTraversal(Cells & icells, Cells & jcells, real_t cycle, bool mutual, real_t remote=1) {
    if (icells.empty() || jcells.empty()) return;               // Quit if either of the cell vectors are empty
    logger::startTimer("Traverse");                             // Start timer
    logger::initTracer();                                       // Initialize tracer
    Ci0 = icells.begin();                                       // Set iterator of target root cell
    Cj0 = jcells.begin();                                       // Set iterator of source root cell
    kernel::Xperiodic = 0;                                      // Periodic coordinate offset
    if (images == 0) {                                          // If non-periodic boundary condition
      traverse(Ci0, Cj0, mutual, remote);                       //  Traverse the tree
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
	for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
	  for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
	    kernel::Xperiodic[0] = ix * cycle;                  //     Coordinate shift for x periodic direction
	    kernel::Xperiodic[1] = iy * cycle;                  //     Coordinate shift for y periodic direction
	    kernel::Xperiodic[2] = iz * cycle;                  //     Coordinate shift for z periodic direction
	    traverse(Ci0, Cj0, false, remote);                  //     Traverse the tree for this periodic image
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
	for (int ix=-prange; ix<=prange; ix++) {                //   Loop over x periodic direction
	  for (int iy=-prange; iy<=prange; iy++) {              //    Loop over y periodic direction
	    for (int iz=-prange; iz<=prange; iz++) {            //     Loop over z periodic direction
	      kernel::Xperiodic[0] = ix * cycle;                //      Coordinate shift for x periodic direction
	      kernel::Xperiodic[1] = iy * cycle;                //      Coordinate shift for y periodic direction
	      kernel::Xperiodic[2] = iz * cycle;                //      Coordinate shift for z periodic direction
	      kernel::P2P(Ci, Cj);                              //      Evaluate P2P kernel
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
    DirectRecursion directRecursion(Ci, Cj, prange, cycle);     // Instantiate recursive functor
    directRecursion();                                          // Recursive call for direct summation
  }

  //! Normalize bodies after direct summation
  void normalize(Bodies & bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG /= B->SRC;                                         //  Normalize by target charge
    }                                                           // End loop over bodies
  }

  //! Print traversal statistics
  void printTraversalData() {
#if COUNT_KERNEL
    if (logger::verbose) {                                      // If verbose flag is true
      std::cout << "--- Traversal stats --------------" << std::endl// Print title
		<< std::setw(logger::stringLength) << std::left //  Set format
		<< "P2P calls"  << " : "                        //  Print title
		<< std::setprecision(0) << std::fixed           //  Set format
		<< numP2P << std::endl                          //  Print number of P2P calls
		<< std::setw(logger::stringLength) << std::left //  Set format
		<< "M2L calls"  << " : "                        //  Print title
		<< std::setprecision(0) << std::fixed           //  Set format
		<< numM2L << std::endl;                         //  Print number of M2L calls
    }                                                           // End if for verbose flag
#endif
  }
#if COUNT_LIST
  void writeList(Cells cells, int mpirank) {
    std::stringstream name;                                     // File name
    name << "list" << std::setfill('0') << std::setw(6)         // Set format
	 << mpirank << ".dat";                                  // Create file name for list
    std::ofstream listFile(name.str().c_str());                 // Open list log file
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {         // Loop over all lists
      listFile << std::setw(logger::stringLength) << std::left  //  Set format
	       << C->ICELL << " " << C->numP2P << " " << C->numM2L << std::endl; // Print list size
    }                                                           // End loop over all lists
  }
#else
  void writeList(Cells, int) {}
#endif
};
#endif
