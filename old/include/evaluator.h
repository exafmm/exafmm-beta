/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef evaluator_h
#define evaluator_h
#include "dataset.h"
#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R)

//! Interface between tree and kernel
template<Equation equation>
class Evaluator : public Dataset<equation> {
private:
  real        timeM2L;                                          //!< M2L execution time
  real        timeM2P;                                          //!< M2P execution time
  real        timeP2P;                                          //!< P2P execution time

protected:
  int         NCRIT;                                            //!< Number of particles per leaf
  C_iter      CiB;                                              //!< icells begin per call
  C_iter      CiE;                                              //!< icells end per call
  Lists       listM2L;                                          //!< M2L interaction list
  Lists       listM2P;                                          //!< M2P interaction list
  Lists       listP2P;                                          //!< P2P interaction list

  int         Iperiodic;                                        //!< Periodic image flag (using each bit for images)
  const int   Icenter;                                          //!< Periodic image flag at center
  Maps        flagM2L;                                          //!< Existance of periodic image for M2L
  Maps        flagM2P;                                          //!< Existance of periodic image for M2P
  Maps        flagP2P;                                          //!< Existance of periodic image for P2P

  real        NP2P;                                             //!< Number of P2P kernel calls
  real        NM2P;                                             //!< Number of M2P kernel calls
  real        NM2L;                                             //!< Number of M2L kernel calls

public:
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::writeTrace;                           //!< Write traces of all events
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::Ci0;                                  //!< Begin iterator for target cells
  using Kernel<equation>::Cj0;                                  //!< Begin iterator for source cells
  using Kernel<equation>::ALPHA;                                //!< Scaling parameter for Ewald summation
  using Kernel<equation>::keysHost;                             //!< Offsets for rangeHost
  using Kernel<equation>::rangeHost;                            //!< Offsets for sourceHost
  using Kernel<equation>::constHost;                            //!< Constants on host
  using Kernel<equation>::sourceHost;                           //!< Sources on host
  using Kernel<equation>::targetHost;                           //!< Targets on host
  using Kernel<equation>::sourceBegin;                          //!< Define map for offset of source cells
  using Kernel<equation>::sourceSize;                           //!< Define map for size of source cells
  using Kernel<equation>::targetBegin;                          //!< Define map for offset of target cells
  using Kernel<equation>::allocate;                             //!< Allocate GPU kernels
  using Kernel<equation>::hostToDevice;                         //!< Copy from host to device
  using Kernel<equation>::deviceToHost;                         //!< Copy from device to host
  using Kernel<equation>::P2M;                                  //!< Evaluate P2M kernel
  using Kernel<equation>::M2M;                                  //!< Evaluate M2M kernel
  using Kernel<equation>::M2L;                                  //!< Evaluate M2L kernel
  using Kernel<equation>::M2P;                                  //!< Evaluate M2P kernel
  using Kernel<equation>::P2P;                                  //!< Evaluate P2P kernel
  using Kernel<equation>::L2L;                                  //!< Evaluate L2L kernel
  using Kernel<equation>::L2P;                                  //!< Evaluate L2P kernel
  using Kernel<equation>::EwaldReal;                            //!< Evaluate Ewald real part
  using Dataset<equation>::initSource;                          //!< Initialize source values
  using Dataset<equation>::initTarget;                          //!< Initialize target values

private:
//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {// If P2P is fastest
      evalP2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {// If M2P is fastest
      evalM2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else {                                                    // If M2L is fastest
      evalM2L(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    }                                                           // End if for kernel selection
#elif TREECODE
    evalM2P(Ci,Cj);                                             // Evaluate on CPU, queue on GPU
#else
    evalM2L(Ci,Cj);                                             // Evalaute on CPU, queue on GPU
#endif
  }

#if QUARK
  inline void interact(C_iter Ci, C_iter Cj, Quark *quark);     //!< interact() function using QUARK
#endif

//! Traverse a single tree using a stack
  void traverseStack(C_iter Ci, C_iter C) {
    CellStack cellStack;                                        // Stack of cells
    cellStack.push(C);                                          // Push pair to queue
    while( !cellStack.empty() ) {                               // While traversal stack is not empty
      C = cellStack.top();                                      //  Get cell from top of stack
      cellStack.pop();                                          //  Pop traversal stack
      for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {// Loop over cell's children
        vect dX = Ci->X - Cj->X - Xperiodic;                    //   Distance vector from source to target
        real Rq = std::sqrt(norm(dX));                          //   Scalar distance
        if( Rq * THETA < Ci->R + Cj->R && Cj->NCHILD == 0 ) {   //   If twigs are close
          evalEwaldReal(Ci,Cj);                                 //    Ewald real part
        } else if( Cj->NCHILD != 0 ) {                          //   If cells are not twigs
          cellStack.push(Cj);                                   //    Push source cell to stack
        }                                                       //   End if for twig cells
      }                                                         //  End loop over cell's children
    }                                                           // End while loop for traversal stack
  }

//! Traverse a pair of trees using a queue
  void traverseQueue(Pair pair) {
    PairQueue pairQueue;                                        // Queue of interacting cell pairs
#if QUARK
    Quark *quark = QUARK_New(1);                                // Initialize QUARK object
#endif
    pairQueue.push_back(pair);                                  // Push pair to queue
    while( !pairQueue.empty() ) {                               // While dual traversal queue is not empty
      pair = pairQueue.front();                                 //  Get interaction pair from front of queue
      pairQueue.pop_front();                                    //  Pop dual traversal queue
      if(splitFirst(pair.first,pair.second)) {                  //  If first cell is larger
        C_iter C = pair.first;                                  //   Split the first cell
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {// Loop over first cell's children
          interact(Ci,pair.second,pairQueue);                   //    Calculate interaction between cells
        }                                                       //   End loop over fist cell's children
      } else {                                                  //  Else if second cell is larger
        C_iter C = pair.second;                                 //   Split the second cell
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {// Loop over second cell's children
          interact(pair.first,Cj,pairQueue);                    //    Calculate interaction betwen cells
        }                                                       //   End loop over second cell's children
      }                                                         //  End if for which cell to split
#if QUARK
      if( int(pairQueue.size()) > 100 ) {                       //  When queue size reaches threshold
        while( !pairQueue.empty() ) {                           //   While dual traversal queue is not empty
          pair = pairQueue.front();                             //    Get interaction pair from front of queue
          pairQueue.pop_front();                                //    Pop dual traversal queue
          interact(pair.first,pair.second,quark);               //    Schedule interact() task on QUARK
        }                                                       //   End while loop for dual traversal queue
      }                                                         //  End if for queue size
#endif
    }                                                           // End while loop for dual traversal queue
#if QUARK
    QUARK_Delete(quark);                                        // Delete QUARK object 
    writeTrace();                                               // Write event trace to file
#endif
  }

//! Get range of periodic images
  int getPeriodicRange() {
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

protected:
//! Get level from cell index
  int getLevel(bigint index) {
    int i = index;                                              // Copy to dummy index
    int level = -1;                                             // Initialize level counter
    while( i >= 0 ) {                                           // While cell index is non-negative
      level++;                                                  //  Increment level
      i -= 1 << 3*level;                                        //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void timeKernels();                                           //!< Time all kernels for auto-tuning

//! Upward phase for periodic cells
  void upwardPeriodic(Cells &jcells) {
    Cells pccells, pjcells;                                     // Periodic center cell and jcell
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      Cell cell;                                                //  New periodic cell at next sublevel
      C_iter C = pccells.end() - 1;                             //  Set previous periodic center cell as source
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              for( int cx=-1; cx<=1; ++cx ) {                   //      Loop over x periodic direction (child)
                for( int cy=-1; cy<=1; ++cy ) {                 //       Loop over y periodic direction (child)
                  for( int cz=-1; cz<=1; ++cz ) {               //        Loop over z periodic direction (child)
                    cell.X[0]  = C->X[0] + (ix * 6 + cx * 2) * C->R;//     Set new x coordinate for periodic image
                    cell.X[1]  = C->X[1] + (iy * 6 + cy * 2) * C->R;//     Set new y cooridnate for periodic image
                    cell.X[2]  = C->X[2] + (iz * 6 + cz * 2) * C->R;//     Set new z coordinate for periodic image
                    cell.M     = C->M;                          //         Copy multipoles to new periodic image
                    cell.NCLEAF = cell.NDLEAF = cell.NCHILD = 0;//         Initialize NCLEAF, NDLEAF, & NCHILD
                    jcells.push_back(cell);                     //         Push cell into periodic jcell vector
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              cell.X[0] = C->X[0] + ix * 2 * C->R;              //      Set new x coordinate for periodic image
              cell.X[1] = C->X[1] + iy * 2 * C->R;              //      Set new y cooridnate for periodic image
              cell.X[2] = C->X[2] + iz * 2 * C->R;              //      Set new z coordinate for periodic image
              cell.M = C->M;                                    //      Copy multipoles to new periodic image
              pjcells.push_back(cell);                          //      Push cell into periodic jcell vector
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increases three times
      pccells.pop_back();                                       //  Pop periodic center cell from vector
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      C_iter Ci = pccells.end() - 1;                            //  Set current cell as target for M2M
      Ci->CHILD = 0;                                            //  Set child cells for periodic M2M
      Ci->NCHILD = 26;                                          //  Set number of child cells for periodic M2M
      evalM2M(pccells,pjcells);                                 // Evaluate periodic M2M kernels for this sublevel
      pjcells.clear();                                          // Clear periodic jcell vector
    }                                                           // End loop over sublevels of tree
  }

//! Traverse tree for periodic cells
  void traversePeriodic(Cells &cells, Cells &jcells) {          
    Xperiodic = 0;                                              // Set periodic coordinate offset
    Iperiodic = Icenter;                                        // Set periodic flag to center
    C_iter Cj = jcells.end()-1;                                 // Initialize iterator for periodic source cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      for( int I=0; I!=26*27; ++I, --Cj ) {                     //  Loop over periodic images (exclude center)
#if TREECODE
        for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) { //   Loop over cells
          if( Ci->NCHILD == 0 ) {                               //    If cell is twig
            evalM2P(Ci,Cj);                                     //     Perform M2P kernel
          }                                                     //    Endif for twig
        }                                                       //   End loop over cells
#else
        C_iter Ci = cells.end() - 1;                            //   Set root cell as target iterator
        evalM2L(Ci,Cj);                                         //   Perform M2P kernel
#endif
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }

public:
//! Constructor
  Evaluator() : Icenter(1 << 13), NP2P(0), NM2P(0), NM2L(0) {}
//! Destructor
  ~Evaluator() {}

//! Random distribution in [-1,1]^3 cube
  void cube(Bodies &bodies, int seed=0, int numSplit=1) {
    srand48(seed);                                              // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = drand48() * 2 * M_PI - M_PI;                  //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Random distribution on r = 1 sphere
  void sphere(Bodies &bodies, int seed=0, int numSplit=1) {
    srand48(seed);                                              // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = drand48() * 2 - 1;                            //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->X));                           //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Uniform distribution on [-pi,pi]^3 lattice (for debugging)
  void lattice(Bodies &bodies) {
    int nx = int(powf(bodies.size()*MPISIZE,1./3));             // Size of lattice per dimension
    assert( int(bodies.size()*MPISIZE) == nx*nx*nx );           // Number of bodies must be cube of integer
    int sign = 1;                                               // Initialize alternating sign
    B_iter B=bodies.begin();                                    // Initialize body iterator
    for( int ix=0; ix!=nx; ++ix ) {                             // Loop over x direction
      for( int iy=0; iy!=nx; ++iy ) {                           //  Loop over y direction
        for( int iz=0; iz!=nx; ++iz, ++B ) {                    //   Loop over z direction
          B->X[0] = -M_PI + (2 * M_PI * ix + M_PI) / nx;        //    Set x coordinate
          B->X[1] = -M_PI + (2 * M_PI * iy + M_PI) / nx;        //    Set y coordinate
          B->X[2] = -M_PI + (2 * M_PI * iz + M_PI) / nx;        //    Set z coordinate
          B->SRC = sign;                                        //    Source alternates between -1 and 1
          B->IBODY = B-bodies.begin();                          //    Tag body with initial index
          B->IPROC = MPIRANK;                                   //    Tag body with initial MPI rank
          sign = -sign;                                         //    Alternate sign in x direction
        }                                                       //   End loop over z direction
        sign = -sign;                                           //   Alternate sign in y direction
      }                                                         //  End loop over y direction
      sign = -sign;                                             //  Alternate sign in z direction
    }                                                           // End loop over x direction
    initTarget(bodies);                                         // Initialize target values
  }

//! Add single list for kernel unit test
  void addM2L(C_iter Cj) {
    listM2L.resize(1);                                          // Resize vector of M2L interation lists
    flagM2L.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2L[0].push_back(Cj);                                   // Push single cell into list
    flagM2L[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

//! Add single list for kernel unit test
  void addM2P(C_iter Cj) {
    listM2P.resize(1);                                          // Resize vector of M2P interation lists
    flagM2P.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2P[0].push_back(Cj);                                   // Push single cell into list
    flagM2P[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

//! Use multipole acceptance criteria to determine whether to approximate, do P2P, or subdivide
  void interact(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    vect dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real Rq = std::sqrt(norm(dX));                              // Scalar distance
    if( Rq * THETA > Ci->R + Cj->R ) {                          // If distance if far enough
      approximate(Ci,Cj);                                       //  Use approximate kernels, e.g. M2L, M2P
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {             // Else if both cells are leafs
      evalP2P(Ci,Cj);                                           //  Use P2P
    } else {                                                    // If cells are close but not leafs
      Pair pair(Ci,Cj);                                         //  Form a pair of cell iterators
      pairQueue.push_back(pair);                                //  Push pair to queue
    }                                                           // End if for multipole acceptance
  }

//! Dual tree traversal
  void traverse(Cells &cells, Cells &jcells) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * 27 * (IMAGES - 1);        //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
#if QUEUE
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagM2L.resize(cells.size());                               // Resize M2L periodic image flag
    flagM2P.resize(cells.size());                               // Resize M2P periodic image flag
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
#endif
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      traverseQueue(pair);                                      //  Traverse a pair of trees
    } else {                                                    // If periodic boundary condition
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
            Iperiodic = 1 << I;                                 //     Set periodic image flag
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            Pair pair(root,jroot);                              //     Form pair of root cells
            traverseQueue(pair);                                //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
#if QUEUE
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2L[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listM2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listP2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listP2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
#endif
    }                                                           // Endif for periodic boundary condition
  }

//! Traverse neighbor cells only (for cutoff based methods)
  void neighbor(Cells &cells, Cells &jcells) {
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
    for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {     // Loop over target cells
      if( Ci->NCHILD == 0 ) {                                   //  If cell is a twig
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              Iperiodic = 1 << I;                               //      Set periodic image flag
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              traverseStack(Ci,jroot);                          //      Traverse the source tree
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
        for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B) {  //   Loop over all leafs in cell
          B->TRG[0] -= M_2_SQRTPI * B->SRC * ALPHA;             //    Self term of Ewald real part
        }                                                       //   End loop over all leafs in cell
      }                                                         //  End if for twig cells
      listP2P[Ci-Ci0].sort();                                   //  Sort interaction list
      listP2P[Ci-Ci0].unique();                                 //  Eliminate duplicate periodic entries
    }                                                           // End loop over target cells
  }

  void setSourceBody();                                         //!< Set source buffer for bodies (for GPU)
  void setSourceCell(bool isM);                                 //!< Set source buffer for cells (for GPU)
  void setTargetBody(Lists lists, Maps flags);                  //!< Set target buffer for bodies (for GPU)
  void setTargetCell(Lists lists, Maps flags);                  //!< Set target buffer for cells (for GPU)
  void getTargetBody(Lists &lists);                             //!< Get body values from target buffer (for GPU)
  void getTargetCell(Lists &lists, bool isM);                   //!< Get cell values from target buffer (for GPU)
  void clearBuffers();                                          //!< Clear GPU buffers

  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate all P2P kernels (all pairs)
  void periodicP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate all P2P kernels (all pairs) for periodic
  void evalP2M(Cells &cells);                                   //!< Evaluate all P2M kernels
  void evalM2M(Cells &cells, Cells &jcells);                    //!< Evaluate all M2M kernels
  void evalM2L(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2L(Cells &cells);                                   //!< Evaluate queued M2L kernels
  void evalM2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2P(Cells &cells);                                   //!< Evaluate queued M2P kernels
  void evalP2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalP2P(Cells &cells);                                   //!< Evaluate queued P2P kernels (near field)
  void evalL2L(Cells &cells);                                   //!< Evaluate all L2L kernels
  void evalL2P(Cells &cells);                                   //!< Evaluate all L2P kernels
  void evalEwaldReal(C_iter Ci, C_iter Cj);                     //!< Evaluate on CPU, queue on GPU
  void evalEwaldReal(Cells &cells);                             //!< Evaluate queued Ewald real kernels
};

#if CPU
#include "../kernel/CPUEvaluator.cxx"
#elif GPU
#include "../kernel/GPUEvaluator.cxx"
#endif

#undef splitFirst
#endif
