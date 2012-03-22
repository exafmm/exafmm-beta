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
#include "../include/kernel.h"
#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT > Cj->RCRIT)

template<Equation equation>
class Evaluator : public Kernel<equation> {
private:
  real timeM2L;                                                 //!< M2L execution time
  real timeM2P;                                                 //!< M2P execution time
  real timeP2P;                                                 //!< P2P execution time

protected:
  bool TOPDOWN;                                                 //!< Flag to indicate top down tree construction
  int  MAXLEVEL;                                                //!< Maximum depth of tree
  real NM2L;                                                    //!< Number of M2L kernel calls
  real NM2P;                                                    //!< Number of M2P kernel calls
  real NP2P;                                                    //!< Number of P2P kernel calls

  int       Iperiodic;                                          //!< Periodic image flag (using each bit for images)
  const int Icenter;                                            //!< Periodic image flag at center

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::writeTrace;                           //!< Write traces of all events
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::Ci0;                                  //!< Begin iterator for target cells
  using Kernel<equation>::Cj0;                                  //!< Begin iterator for source cells
  using Kernel<equation>::P2M;                                  //!< Evaluate P2M kernel
  using Kernel<equation>::M2M;                                  //!< Evaluate M2M kernel
  using Kernel<equation>::M2L;                                  //!< Evaluate M2L kernel
  using Kernel<equation>::M2P;                                  //!< Evaluate M2P kernel
  using Kernel<equation>::P2P;                                  //!< Evaluate P2P kernel
  using Kernel<equation>::L2L;                                  //!< Evaluate L2L kernel
  using Kernel<equation>::L2P;                                  //!< Evaluate L2P kernel

private:
  inline void approximate(C_iter Ci, C_iter Cj) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {
      evalP2P(Ci,Cj);
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {
      evalM2P(Ci,Cj);
    } else {
      evalM2L(Ci,Cj);
    }
#elif TREECODE
    evalM2P(Ci,Cj);
#else
    evalM2L(Ci,Cj);
#endif
  }

//! Get iterator for root cell
  C_iter getRootCell(Cells &cells) {
    if( TOPDOWN ) {                                             // If tree was constructed top down
      return cells.begin();                                     //  Root cell is stored first
    } else {                                                    // If tree was constructed bottom up
      return cells.end() - 1;                                   //  Root cell is stored last
    }                                                           // End if for tree construction
  }

//! Set error controlled theta and corresponding critical cell radius
  void setRcrit(Cells &cells) {
    C_iter root = getRootCell(cells);
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(root->M[0]),1.0/3);
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      real a = c * pow(std::abs(C->M[0]),1.0/3);
      real x = 1.0 / THETA;
      for( int i=0; i<5; ++i ) {
        real f = x * x - 2 * x + 1 - a * pow(x,-P);
        real df = (P + 2) * x - 2 * (P + 1) + P / x;
        x -= f / df;
      }
      C->RCRIT *= x;
    }
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
  void upwardPass(Cells &cells) {
    startTimer("Upward pass");
    evalP2M(cells);
    evalM2M(cells,cells);
#ifndef SPHERICAL
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass",printNow);
  }

  void traverse(Cells &icells, Cells &jcells) {
    Ci0 = icells.begin();                                       // Begin iterator for target cells
    Cj0 = jcells.begin();                                       // Begin iterator for source cells
    Pair pair;                                                  // Form pair of cells
    if( TOPDOWN ) {                                             // If tree was constructed top down
      pair = make_pair(icells.begin(),jcells.begin());          //  The first cell is the root cell
    } else {                                                    // If tree was constructed bottom up
      pair = make_pair(icells.end()-1,jcells.end()-1);          //  The last cell is the root cell
    }                                                           // Endif for tree construction
    PairQueue pairQueue;                                        // Queue of interacting cell pairs
    pairQueue.push(pair);                                       // Push pair of root cells to queue
#if QUARK
    Quark *quark = QUARK_New(4);                                // Initialize QUARK object
#endif
    while( !pairQueue.empty() ) {                               // While dual traversal queue is not empty
      pair = pairQueue.front();                                 //  Get interaction pair from front of queue
      pairQueue.pop();                                          //  Pop dual traversal queue
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
      if( pairQueue.size() > 100 ) {                            //  When queue size reaches threshold
        while( !pairQueue.empty() ) {                           //   While dual traversal queue is not empty
          pair = pairQueue.front();                             //    Get interaction pair from front of queue
          pairQueue.pop();                                      //    Pop dual traversal queue
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

//! Traverse tree for periodic cells
  void traversePeriodic(Cells &icells, Cells &jcells) {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Set periodic coordinate offset
    Iperiodic = Icenter;                                        // Set periodic flag to center
    Cells ccells(1), pcells(1), ncells(27);                     // Create cells
    C_iter C = ccells.begin();                                  // ccells is periodic center cell
    C_iter Cj = pcells.begin();                                 // pcells is periodic source cell
    C_iter Ci;                                                  // Initialize iterator for periodic target cell
    if( TOPDOWN ) {                                             // If tree was constructed top down
      C = jcells.begin();                                       //  The first cell is the first periodic cell
    } else {                                                    // If tree was constructed bottom up
      C = jcells.end() - 1;                                     //  The last cell is the first periodic cell
    }                                                           // Endif for tree construction
    C->CHILD = 0;                                               // Set child cells for periodic center cell
    C->NCHILD = 27;                                             // Set number of child cells for periodic center cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              for( int cx=-1; cx<=1; ++cx ) {                   //      Loop over x periodic direction (child)
                for( int cy=-1; cy<=1; ++cy ) {                 //       Loop over y periodic direction (child)
                  for( int cz=-1; cz<=1; ++cz ) {               //        Loop over z periodic direction (child)
                    Cj->X[0]  = C->X[0] + (ix * 6 + cx * 2) * C->R;//     Set new x coordinate for periodic image
                    Cj->X[1]  = C->X[1] + (iy * 6 + cy * 2) * C->R;//     Set new y cooridnate for periodic image
                    Cj->X[2]  = C->X[2] + (iz * 6 + cz * 2) * C->R;//     Set new z coordinate for periodic image
                    Cj->M     = C->M;                           //         Copy multipoles to new periodic image
#if TREECODE
                    for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {//  Loop over target cells
                      if( Ci->NCHILD == 0 ) {                   //         If cell is twig
                        evalM2P(Ci,Cj);                         //          Perform M2P kernel
                      }                                         //         Endif for twig
                    }                                           //        End loop over target cells
#else
                   if( TOPDOWN ) {                              //        If tree was constructed top down
                     Ci = icells.begin();                       //         Set first cell as periodic target
                   } else {                                     //        If tree was constructed bottom up
                     Ci = icells.end() - 1;                     //         Set last cell as periodic target
                   }                                            //        Endif for tree constructed
                   evalM2L(Ci,Cj);                              //        Perform M2P kernel
#endif
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      C_iter Cp = ncells.begin();                               //  Iterator for periodic neighbor cells
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++Cp ) {                 //    Loop over z periodic direction
            Cp->X[0] = C->X[0] + ix * 2 * C->R;                 //     Set new x coordinate for periodic image
            Cp->X[1] = C->X[1] + iy * 2 * C->R;                 //     Set new y cooridnate for periodic image
            Cp->X[2] = C->X[2] + iz * 2 * C->R;                 //     Set new z coordinate for periodic image
            Cp->M    = C->M;                                    //     Copy multipoles to new periodic image
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      evalM2M(ccells,ncells);                                   //  Evaluate periodic M2M kernels for this sublevel
      C->R *= 3;                                                //  Increases center cell size three times
    }                                                           // End loop over sublevels of tree
    stopTimer("Traverse periodic",printNow);                    // Stop timer
  }

  void downwardPass(Cells &cells) {
    startTimer("Downward pass");
    evalL2L(cells);
    evalL2P(cells);
    stopTimer("Downward pass",printNow);
  }

  void timeKernels();

public:
  Evaluator() : timeM2L(0), timeM2P(0), timeP2P(0),
                TOPDOWN(false), MAXLEVEL(0),
                NM2L(0), NM2P(0), NP2P(0),
                Iperiodic(0), Icenter(1 << 13) {}
  ~Evaluator() {}

  void printTreeData(Cells &cells) {
    C_iter root = getRootCell(cells);
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Root center          : " << root->X              << std::endl;
    std::cout << "Root radius          : " << root->R              << std::endl;
    std::cout << "Bodies               : " << root->NDLEAF         << std::endl;
    std::cout << "Cells                : " << cells.size()         << std::endl;
    std::cout << "Tree depth           : " << MAXLEVEL             << std::endl;
    std::cout << "Total charge         : " << std::abs(root->M[0]) << std::endl;
    std::cout << "M2L calls            : " << NM2L                 << std::endl;
    std::cout << "M2P calls            : " << NM2P                 << std::endl;
    std::cout << "P2P calls            : " << NP2P                 << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
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

//! Use multipole acceptance criteria to determine whether to approximate, do P2P, or subdivide
  void interact(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    vect dist = Ci->X - Cj->X - Xperiodic;                      // Distance vector from source to target
    real R2 = norm(dist);                                       // R^2
    if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {      // If distance if far enough
      approximate(Ci,Cj);                                       //  Use approximate kernels, e.g. M2L, M2P
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {             // If both cells are leafs
      evalP2P(Ci,Cj);                                           //  Use P2P
    } else {                                                    // If cells are close but not leafs 
      Pair pair(Ci,Cj);                                         //  Form a pair of cell iterators
      pairQueue.push(pair);                                     //  Push pair to queue
    }                                                           // End if for multipole acceptance
  }

#if QUARK
  inline void interact(C_iter Ci, C_iter Cj, Quark *quark);     //!< interact() function using QUARK
#endif
  void direct(Bodies &ibodies, Bodies &jbodies);                //!< Evaluate direct summation
  inline void evalP2M(Cells &cells);                            //!< Evaluate all P2M kernels
  inline void evalM2M(Cells &cells, Cells &jcells);             //!< Evaluate all M2M kernels
  void evalM2L(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalP2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  inline void evalL2L(Cells &cells);                            //!< Evaluate all L2L kernels
  inline void evalL2P(Cells &cells);                            //!< Evaluate all L2P kernels

};

#include "CPUEvaluator.cxx"

#undef splitFirst
#endif
