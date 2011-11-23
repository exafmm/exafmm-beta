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
#include "kernel.h"
#include "dataset.h"

//! Interface between tree and kernel
template<Equation equation>
class Evaluator : public Kernel<equation>, public Dataset<equation> {
protected:
  C_iter      CI0;                                              //!< icells.begin()
  C_iter      CIB;                                              //!< icells begin per call
  C_iter      CIE;                                              //!< icells end per call
  C_iter      CJ0;                                              //!< jcells.begin()
  C_iter      CJB;                                              //!< jcells begin per call
  C_iter      CJE;                                              //!< jcells end per call
  Pairs       pairs;                                            //!< Stack of interacting cell pairs
  Lists       listM2L;                                          //!< M2L interaction list
  Lists       listM2P;                                          //!< M2P interaction list
  Lists       listP2P;                                          //!< P2P interaction list
  real        timeM2L;                                          //!< M2L execution time
  real        timeM2P;                                          //!< M2P execution time
  real        timeP2P;                                          //!< P2P execution time

  int         Iperiodic;                                        //!< Periodic image flag (using each bit for images)
  const int   Icenter;                                          //!< Periodic image flag at center
  Maps        flagM2L;                                          //!< Existance of periodic image for M2L
  Maps        flagM2P;                                          //!< Existance of periodic image for M2P
  Maps        flagP2P;                                          //!< Existance of periodic image for P2P

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::BI0;                                  //!< Target bodies begin iterator
  using Kernel<equation>::BIN;                                  //!< Target bodies end iterator
  using Kernel<equation>::BJ0;                                  //!< Source bodies begin iterator
  using Kernel<equation>::BJN;                                  //!< Source bodies end iterator
  using Kernel<equation>::CI;                                   //!< Target cell iterator
  using Kernel<equation>::CJ;                                   //!< Source cell iterator
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::Xperiodic;                            //!< Coordinate offset of periodic image
  using Kernel<equation>::keysHost;                             //!< Offsets for rangeHost
  using Kernel<equation>::rangeHost;                            //!< Offsets for sourceHost
  using Kernel<equation>::constHost;                            //!< Constants on host
  using Kernel<equation>::sourceHost;                           //!< Sources on host
  using Kernel<equation>::targetHost;                           //!< Targets on host
  using Kernel<equation>::sourceBegin;                          //!< Define map for offset of source cells
  using Kernel<equation>::sourceSize;                           //!< Define map for size of source cells
  using Kernel<equation>::targetBegin;                          //!< Define map for offset of target cells
  using Kernel<equation>::NP2P;                                 //!< Number of P2P kernel calls
  using Kernel<equation>::NM2P;                                 //!< Number of M2P kernel calls
  using Kernel<equation>::NM2L;                                 //!< Number of M2L kernel calls
  using Kernel<equation>::P2M;                                  //!< Evaluate P2M kernel
  using Kernel<equation>::M2M;                                  //!< Evaluate M2M kernel
  using Kernel<equation>::M2M_CPU;                              //!< Evaluate M2M kernel on CPU
  using Kernel<equation>::M2L;                                  //!< Evaluate M2L kernel
  using Kernel<equation>::M2P;                                  //!< Evaluate M2P kernel
  using Kernel<equation>::P2P;                                  //!< Evaluate P2P kernel
  using Kernel<equation>::P2P_CPU;                              //!< Evaluate P2P kernel on CPU
  using Kernel<equation>::L2L;                                  //!< Evaluate L2L kernel
  using Kernel<equation>::L2P;                                  //!< Evaluate L2P kernel
  using Dataset<equation>::initSource;                          //!< Initialize source values
  using Dataset<equation>::initTarget;                          //!< Initialize target values

private:
//! Tree walk for treecode
  void treecode(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0) {                   // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Ci->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        testMACM2P(CI0+Ci->CHILD[i],Cj);                        //   Test multipole acceptance criteria for M2P kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        testMACM2P(Ci,CJ0+Cj->CHILD[i]);                        //   Test multipole acceptance criteria for M2P kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

//! Tree walk for FMM
  void FMM(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0 ) {                  // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R) ) {// If source is twig or target is larger
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        testMACM2L(CI0+Ci->CHILD[i],Cj);                        //   Test multipole acceptance criteria for M2L kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        testMACM2L(Ci,CJ0+Cj->CHILD[i]);                        //   Test multipole acceptance criteria for M2L kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

//! Tree walk for treecode-FMM hybrid
  void hybrid(C_iter Ci, C_iter Cj) {
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0 ) {                  // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        testMACP2P(Ci,Cj);                                      //  Test multipole acceptance criteria for P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R) ) {// If source is twig or target is larger
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        int Ni = (CI0+Ci->CHILD[i])->NLEAF;                     //   Number of target leafs
        int Nj = Cj->NLEAF;                                     //   Number of source leafs
        if( timeP2P*Nj < timeM2P && timeP2P*Ni*Nj < timeM2L ) { //   If P2P is fastest
          testMACP2P(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for P2P kernel
        } else if ( timeM2P < timeP2P*Nj && timeM2P*Ni < timeM2L ) {// If M2P is fastest
          testMACM2P(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for M2P kernel
        } else {                                                //   If M2L is fastest
          testMACM2L(CI0+Ci->CHILD[i],Cj);                      //    Test multipole acceptance criteria for M2L kernel
        }                                                       //   End if for fastest kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        int Ni = Ci->NLEAF;                                     //   Number of target leafs
        int Nj = (CJ0+Cj->CHILD[i])->NLEAF;                     //   Number of source leafs
        if( timeP2P*Nj < timeM2P && timeP2P*Ni*Nj < timeM2L ) { //   If P2P is fastest
          testMACP2P(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for P2P kernel
        } else if ( timeM2P < timeP2P*Nj && timeM2P*Ni < timeM2L ) {// If M2P is fastest
          testMACM2P(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for M2P kernel
        } else {                                                //   If M2L is fastest
          testMACM2L(Ci,CJ0+Cj->CHILD[i]);                      //    Test multipole acceptance criteria for M2L kernel
        }                                                       //   End if for fastest kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

protected:
//! Get level from cell index
  int getLevel(bigint index) {
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void timeKernels();                                           //!< Time all kernels for auto-tuning

public:
//! Constructor
  Evaluator() : Icenter(1 << 13) {}
//! Destructor
  ~Evaluator() {}

//! Random distribution in [-1,1]^3 cube
  void random(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;   //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Random distribution on r = 1 sphere
  void sphere(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 - 1;             //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->X));                           //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Uniform distribution on [-1,1]^3 lattice (for debugging)
  void lattice(Bodies &bodies) {
    int level = int(log(bodies.size()*MPISIZE+1.)/M_LN2/3);     // Level of tree
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int d = 0, l = 0;                                         //  Initialize dimension and level
      int index = MPIRANK * bodies.size() + (B-bodies.begin()); //  Set index of body iterator
      vec<3,int> nx = 0;                                        //  Initialize 3-D cell index
      while( index != 0 ) {                                     //  Deinterleave bits while index is nonzero
        nx[d] += (index % 2) * (1 << l);                        //   Add deinterleaved bit to 3-D cell index
        index >>= 1;                                            //   Right shift the bits
        d = (d+1) % 3;                                          //   Increment dimension
        if( d == 0 ) l++;                                       //   If dimension is 0 again, increment level
      }                                                         //  End while loop for deinterleaving bits
      for( d=0; d!=3; ++d ) {                                   //  Loop over dimensions
        B->X[d] = -1 + (2 * nx[d] + 1.) / (1 << level);         //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
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

//! Get range of periodic images
  int getPeriodicRange() {
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

//! Create periodic images of bodies
  Bodies periodicBodies(Bodies &bodies) {
    Bodies jbodies;                                             // Vector for periodic images of bodies
    int prange = getPeriodicRange();                            // Get range of periodic images
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {//    Loop over bodies
            Body body = *B;                                     //     Copy current body
            body.X[0] += ix * 2 * R0;                           //     Shift x position
            body.X[1] += iy * 2 * R0;                           //     Shift y position
            body.X[2] += iz * 2 * R0;                           //     Shift z position
            jbodies.push_back(body);                            //     Push shifted body into jbodies
          }                                                     //    End loop over bodies
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    return jbodies;                                             // Return vector for periodic images of bodies
  }

//! Traverse tree to get interaction list
  void traverse(Cells &cells, Cells &jcells, int method) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * (IMAGES - 1);             //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    CI0 = cells.begin();                                        // Set begin iterator for target cells
    CJ0 = jcells.begin();                                       // Set begin iterator for source cells
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagM2L.resize(cells.size());                               // Resize M2L periodic image flag
    flagM2P.resize(cells.size());                               // Resize M2P periodic image flag
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      pairs.push(pair);                                         //  Push pair to stack
      while( !pairs.empty() ) {                                 //  While interaction stack is not empty
        pair = pairs.top();                                     //   Get interaction pair from top of stack
        pairs.pop();                                            //   Pop interaction stack
        switch (method) {                                       //   Swtich between methods
        case 0 : treecode(pair.first,pair.second); break;       //    0 : treecode
        case 1 : FMM(pair.first,pair.second);      break;       //    1 : FMM
        case 2 : hybrid(pair.first,pair.second);   break;       //    2 : hybrid
        }                                                       //   End switch between methods
      }                                                         //  End while loop for interaction stack
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
            pairs.push(pair);                                   //     Push pair to stack
            while( !pairs.empty() ) {                           //     While interaction stack is not empty
              pair = pairs.top();                               //      Get interaction pair from top of stack
              pairs.pop();                                      //      Pop interaction stack
              switch (method) {                                 //      Swtich between methods
              case 0 : treecode(pair.first,pair.second); break; //      0 : treecode
              case 1 : FMM(pair.first,pair.second);      break; //      1 : FMM
              case 2 : hybrid(pair.first,pair.second);   break; //      2 : hybrid
              }                                                 //      End switch between methods
            }                                                   //     End while loop for interaction stack
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-CI0].sort();                                 //  Sort interaction list
        listM2L[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
        listM2P[Ci-CI0].sort();                                 //  Sort interaction list
        listM2P[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
        listP2P[Ci-CI0].sort();                                 //  Sort interaction list
        listP2P[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
    }                                                           // Endif for periodic boundary condition
  }

//! Upward phase for periodic cells
  void upwardPeriodic(Cells &jcells) {
    Cells pccells, pjcells;                                     // Periodic jcells for M2L/M2P & M2M
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      C_iter C = pccells.end() - 1;                             //  Set previous periodic cell as source
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            Cell cell;                                          //     New periodic jcell for M2M
            cell.X[0] = C->X[0] + ix * 2 * C->R;                //     Set new x coordinate for periodic image
            cell.X[1] = C->X[1] + iy * 2 * C->R;                //     Set new y cooridnate for periodic image
            cell.X[2] = C->X[2] + iz * 2 * C->R;                //     Set new z coordinate for periodic image
            cell.M = C->M;                                      //     Copy multipoles to new periodic image
            pjcells.push_back(cell);                            //     Push cell into periodic jcell vector
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cell cell;                                                //  New periodic cell at next sublevel
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increase three times
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      CI = pccells.end() - 1;                                   //  Set current cell as target for M2M
      while( !pjcells.empty() ) {                               //  While there are periodic jcells remaining
        CJ = pjcells.end() - 1;                                 //   Set current jcell as source for M2M
        M2M_CPU();                                              //   Perform M2M_CPU kernel
        pjcells.pop_back();                                     //   Pop last element from periodic jcell vector
      }                                                         //  End while for remaining periodic jcells
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              cell.X[0]  = CI->X[0] + ix * 2 * CI->R;           //      Set new x coordinate for periodic image
              cell.X[1]  = CI->X[1] + iy * 2 * CI->R;           //      Set new y cooridnate for periodic image
              cell.X[2]  = CI->X[2] + iz * 2 * CI->R;           //      Set new z coordinate for periodic image
              cell.M     = CI->M;                               //      Copy multipoles to new periodic image
              cell.NLEAF = cell.NCHILD = 0;                     //      Initialize NLEAF & NCHILD
              jcells.push_back(cell);                           //      Push cell into periodic jcell vector
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }

  void setSourceBody();                                         //!< Set source buffer for bodies
  void setSourceCell(bool isM);                                 //!< Set source buffer for cells
  void setTargetBody(Lists lists, Maps flags);                  //!< Set target buffer for bodies
  void setTargetCell(Lists lists, Maps flags);                  //!< Set target buffer for cells
  void getTargetBody(Lists &lists);                             //!< Get body values from target buffer
  void getTargetCell(Lists &lists, bool isM);                   //!< Get cell values from target buffer
  void clearBuffers();                                          //!< Clear GPU buffers

  void testMACP2P(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for P2P kernel
  void testMACM2L(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for M2L kernel
  void testMACM2P(C_iter Ci, C_iter Cj);                        //!< Test multipole acceptance criteria for M2P kernel
  void traversePeriodic(Cells &cells, Cells &jcells, int method);//!< Traverse tree for periodic cells
  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate P2P kernel (all pairs)
  void evalP2M(Cells &twigs);                                   //!< Evaluate P2M kernel
  void evalM2M(Cells &cells);                                   //!< Evaluate M2M kernel
  void evalM2L(Cells &cells, bool kernel=false);                //!< Evaluate M2L kernel
  void evalM2P(Cells &cells, bool kernel=false);                //!< Evaluate M2P kernel
  void evalP2P(Cells &cells, bool kernel=false);                //!< Evaluate P2P kernel (near field)
  void evalL2L(Cells &cells);                                   //!< Evaluate L2L kernel
  void evalL2P(Cells &cells);                                   //!< Evaluate L2P kernel
};
#if cpu
#include "../kernel/cpuEvaluator.cxx"
#elif gpu
#include "../kernel/gpuEvaluator.cxx"
#endif

#endif
