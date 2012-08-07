#ifndef evaluator_h
#define evaluator_h
#include "cartesian.h"
#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT >= Cj->RCRIT)
#if COUNT
#define count(N) N++
#else
#define count(N)
#endif

class Evaluator : public Kernel {
private:
  real_t timeP2P;                                               //!< P2P execution time
  real_t timeM2P;                                               //!< M2P execution time
  real_t timeM2L;                                               //!< M2L execution time

protected:
  real_t NP2P;                                                  //!< Number of P2P kernel calls
  real_t NM2P;                                                  //!< Number of M2P kenrel calls
  real_t NM2L;                                                  //!< Number of M2L kernel calls

private:
//! Calculate Bmax
  real_t getBmax(vec3 const&X, C_iter C) const {
    real_t rad = C->R;                                          // Radius of cell
    real_t dx = rad+std::abs(X[0]-C->X[0]);                     // Add x distance from center of mass
    real_t dy = rad+std::abs(X[1]-C->X[1]);                     // Add y distance from center of mass
    real_t dz = rad+std::abs(X[2]-C->X[2]);                     // Add z distance from center of mass
    return std::sqrt( dx*dx + dy*dy + dz*dz );                  // Return scalar distance
  }

//! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj, bool mutual) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {// If P2P is fastest
      P2P(Ci,Cj,mutual);                                        //  P2P kernel
      count(NP2P);                                              //  Increment P2P counter
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {// If M2P is fastest
      M2P(Ci,Cj,mutual);                                        //  M2P kernel
      count(NM2P);                                              //  Increment M2P counter
    } else {                                                    // If M2L is fastest
      M2L(Ci,Cj,mutual);                                        //  M2L kernel
      count(NM2L);                                              //  Increment M2L counter
    }                                                           // End if for fastest kernel
#elif TREECODE
    M2P(Ci,Cj,mutual);                                          // M2P kernel
    count(NM2P);                                                // Increment M2P counter
#else
    M2L(Ci,Cj,mutual);                                          // M2L kernel
    count(NM2L);                                                // Increment M2L counter
#endif
  }

//! Split cell and push children to queue
  void splitCell(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    if(splitFirst(Ci,Cj)) {                                     // If first cell is larger or equal
      for( C_iter CC=Ci0+Ci->CHILD; CC!=Ci0+Ci->CHILD+Ci->NCHILD; ++CC ) {// Loop over first cell's children
        Pair pair(CC,Cj);                                       //   Form a pair of cell iterators
        pairQueue.push_back(pair);                              //   Push pair to queue
      }                                                         //  End loop over first cell's children
    } else {                                                    // If second cell is larger
      for( C_iter CC=Cj0+Cj->CHILD; CC!=Cj0+Cj->CHILD+Cj->NCHILD; ++CC ) {// Loop over second cell's children
        Pair pair(Ci,CC);                                       //   Form a pair of cell iterators
        pairQueue.push_back(pair);                              //   Push pair to queue
      }                                                         //  End loop over second cell's children
    }                                                           // End if for which cell to split
  }

//! Split cell and call function recursively for child
#if MTHREADS
  void splitCell(C_iter Ci, C_iter Cj, bool mutual) {
    if(splitFirst(Ci,Cj)) {                                     // If first cell is larger or equal
      task_group tg;                                            //  Create task group
      for( C_iter CC=Ci0+Ci->CHILD; CC!=Ci0+Ci->CHILD+Ci->NCHILD; ++CC ) {// Loop over first cell's children 
        if( CC->NDLEAF > 10000 ) tg.run([=]{applyMAC(CC,Cj,mutual);});// Create a new task and recurse
        else applyMAC(CC,Cj,mutual);                            //    Recurse without creating new task
      }                                                         //   End loop over first cell's children
      tg.wait();                                                //  Wait for task group to finish
    } else {                                                    // If second cell is larger
      for( C_iter CC=Cj0+Cj->CHILD; CC!=Cj0+Cj->CHILD+Cj->NCHILD; ++CC ) {// Loop over second cell's children
        applyMAC(Ci,CC,mutual);                                 //   Recurse without creating new task
      }                                                         //  End loop over second cell's children
    }                                                           // End if for which cell to split
  }
#elif OPENMP
  void splitCell(C_iter Ci, C_iter Cj, bool mutual) {
    if(splitFirst(Ci,Cj)) {                                     // If first cell is larger or equal
      for( C_iter CC=Ci0+Ci->CHILD; CC!=Ci0+Ci->CHILD+Ci->NCHILD; ++CC ) {// Loop over first cell's children 
        if( CC->NDLEAF > 10000 ) {
#pragma omp task
          applyMAC(CC,Cj,mutual);                               // Create a new task and recurse
        }
        else applyMAC(CC,Cj,mutual);                            //    Recurse without creating new task
      }                                                         //   End loop over first cell's children
    } else {                                                    // If second cell is larger
      for( C_iter CC=Cj0+Cj->CHILD; CC!=Cj0+Cj->CHILD+Cj->NCHILD; ++CC ) {// Loop over second cell's children
        applyMAC(Ci,CC,mutual);                                 //   Recurse without creating new task
      }                                                         //  End loop over second cell's children
    }                                                           // End if for which cell to split
  }
#endif

  void applyMAC(C_iter Ci, C_iter Cj, PairQueue &pairQueue, bool mutual) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;
    real_t R2 = norm(dX);
#if DUAL
    {
#else
    if(Ci->RCRIT != Cj->RCRIT) {
      splitCell(Ci,Cj,pairQueue);
    } else {
#endif
      if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
        approximate(Ci,Cj,mutual);
      } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
        if( Cj->NCLEAF == 0 ) {
          approximate(Ci,Cj,mutual);
        } else {
          P2P(Ci,Cj,mutual);
          count(NP2P);
        }
      } else {
        splitCell(Ci,Cj,pairQueue);
      }
    }
  }

#if MTHREADS
  void applyMAC(C_iter Ci, C_iter Cj, bool mutual) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;
    real_t R2 = norm(dX);
#if DUAL
    {
#else
    if(Ci->RCRIT != Cj->RCRIT) {
      splitCell(Ci,Cj,mutual);
    } else {
#endif
      if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
        approximate(Ci,Cj,mutual);
      } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
        if( Cj->NCLEAF == 0 ) {
          approximate(Ci,Cj,mutual);
        } else {
          P2P(Ci,Cj,mutual);
          count(NP2P);
        }
      } else {
        splitCell(Ci,Cj,mutual);
      }
    }
  }
#elif OPENMP
  void applyMAC(C_iter Ci, C_iter Cj, bool mutual) {
    vec3 dX = Ci->X - Cj->X - Xperiodic;
    real_t R2 = norm(dX);
#if DUAL
    {
#else
    if(Ci->RCRIT != Cj->RCRIT) {
      splitCell(Ci,Cj,mutual);
    } else {
#endif
      if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
        approximate(Ci,Cj,mutual);
      } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
        if( Cj->NCLEAF == 0 ) {
          approximate(Ci,Cj,mutual);
        } else {
          P2P(Ci,Cj,mutual);
          count(NP2P);
        }
      } else {
        splitCell(Ci,Cj,mutual);
      }
    }
  }
#endif

protected:
//! Set center of expansion to center of mass
  void setCenter(C_iter C) const {
    real_t m = 0;                                               // Initialize mass
    vec3 X = 0;                                                 // Initialize coordinates
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {        // Loop over leafs
      m += B->SRC;                                              //  Accumulate mass
      X += B->X * B->SRC;                                       //  Accumulate dipole
    }                                                           // End loop over leafs
    for( C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; ++c ) {// Loop over cell's children
      m += std::abs(c->M[0]);                                   //  Accumulate mass
      X += c->X * std::abs(c->M[0]);                            //  Accumulate dipole
    }                                                           // End loop over child cells
    X /= m;                                                     // Center of mass
#if USE_BMAX
    C->R = getBmax(X,C);                                        // Use Bmax as cell radius
#endif
#if COMcenter
    C->X = X;                                                   // Use center of mass as center of expansion
#endif
  }

//! Push iterm to cell queue
  void pushCell(C_iter C, CellQueue &cellQueue) {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {                      // If cell has no child or has few leafs
      P2P(C);                                                   //  P2P kernel
      count(NP2P);                                              //  Increment P2P counter
    } else {                                                    // Else
      cellQueue.push(C);                                        //  Push cell to queue
    }                                                           // End if for pushing cell
  }

//! Dual tree traversal of a single tree using a queue of cells
  void traverse(CellQueue &cellQueue) {
    PairQueue pairQueue;                                        // Queue of cell pairs
    while( !cellQueue.empty() ) {                               // While cell queue is not empty
      C_iter C = cellQueue.front();                             //  Get cell from front of queue
      cellQueue.pop();                                          //  Pop this front item from queue
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {// Loop over target cell's children
        pushCell(Ci,cellQueue);                                 //   Push target cell's child to queue
        for( C_iter Cj=Ci+1; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {//  Loop over upper half of source cell's children
          applyMAC(Ci,Cj,pairQueue,true);                       //    Apply multipole acceptance criterion to pair 
        }                                                       //   End loop over source cell's children
      }                                                         //  End loop over target cell's children
      traverse(pairQueue,true);                                 //  Traverse tree for all the pairs that were created
    }                                                           // End while loop over cell queue
  }

//! Dual tree traversal of a pair of trees using a queue of pairs
  void traverse(PairQueue &pairQueue, bool mutual=false) {
#if MTHREADS
    Pair pair = pairQueue.front();
    pairQueue.pop_back();
    applyMAC(pair.first,pair.second,mutual);
#elif OPENMP
    Pair pair = pairQueue.front();
    pairQueue.pop_back();
#pragma omp parallel
#pragma omp single
    applyMAC(pair.first,pair.second,mutual);
#else
#if QUARK
    Quark *quark = QUARK_New(12);
#endif
    while( !pairQueue.empty() ) {
      Pair pair = pairQueue.front();
      pairQueue.pop_front();
      applyMAC(pair.first,pair.second,pairQueue,mutual);
#if QUARK
      if( int(pairQueue.size()) > Ci0->NDLEAF / 50 ) {
        while( !pairQueue.empty() ) {
          pair = pairQueue.front();
          pairQueue.pop_front();
          traverseBranch(pair.first,pair.second,quark,mutual);
        }
      }
#endif // QUARK
    }
#if QUARK
    QUARK_Delete(quark);
    writeTrace();
#endif
#endif // MTHREADS
  }

//! Tree traversal of periodic cells
  void traversePeriodic(real_t R) {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Periodic coordinate offset
    Cells pcells(28);                                           // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *Cj0;                                                 // Copy values from source root
    Ci->CHILD = 0;                                              // Child cells for periodic center cell
    Ci->NCHILD = 27;                                            // Number of child cells for periodic center cell
    C_iter C0 = Cj0;                                            // Placeholder for Cj0
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              for( int cx=-1; cx<=1; ++cx ) {                   //      Loop over x periodic direction (child)
                for( int cy=-1; cy<=1; ++cy ) {                 //       Loop over y periodic direction (child)
                  for( int cz=-1; cz<=1; ++cz ) {               //        Loop over z periodic direction (child)
                    Xperiodic[0] = (ix * 3 + cx) * 2 * R;       //         Coordinate offset for x periodic direction
                    Xperiodic[1] = (iy * 3 + cy) * 2 * R;       //         Coordinate offset for y periodic direction
                    Xperiodic[2] = (iz * 3 + cz) * 2 * R;       //         Coordinate offset for z periodic direction
                    M2L(Ci0,Ci,false);                          //         Perform M2L kernel
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cj0 = pcells.begin();                                     //  Redefine Cj0 for M2M
      C_iter Cj = Cj0;                                          //  Iterator for periodic neighbor cells
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++Cj ) {                 //    Loop over z periodic direction
            Cj->X[0] = Ci->X[0] + ix * 2 * R;                   //     Set new x coordinate for periodic image
            Cj->X[1] = Ci->X[0] + iy * 2 * R;                   //     Set new y cooridnate for periodic image
            Cj->X[2] = Ci->X[0] + iz * 2 * R;                   //     Set new z coordinate for periodic image
            Cj->M    = Ci->M;                                   //     Copy multipoles to new periodic image
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Ci->M = 0;                                                //  Reset multipoles of periodic parent
      real_t Rmax = 0;                                          //  Dummy parameter for calling M2M
      setCenter(Ci);                                            //  Set center of mass for periodic parent
      M2M(Ci,Rmax);                                             //  Evaluate periodic M2M kernels for this sublevel
      R *= 3;                                                   //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
    stopTimer("Traverse periodic",printNow);                    // Stop timer
  }

public:
  Evaluator() : NP2P(0), NM2P(0), NM2L(0) {}
  ~Evaluator() {}

//! Traverse branch of tree
  void traverseBranch(C_iter Ci, C_iter Cj, bool mutual) {
    PairQueue privateQueue;                                     // Queue of interacting cell pairs
    Pair pair(Ci,Cj);                                           // Form a pair of cell iterators
    privateQueue.push_back(pair);                               // Push pair to queue
    while( !privateQueue.empty() ) {                            // While dual traversal queue is not empty
      pair = privateQueue.front();                              //  Get interaction pair from front of queue
      privateQueue.pop_front();                                 //  Pop dual traversal queue
      applyMAC(pair.first,pair.second,privateQueue,mutual);     //    Apply multipole acceptance criterion to pair
    }                                                           // End while loop for traversal queue
  }

  void timeKernels(bool mutual=true) {
    Bodies ibodies(1000), jbodies(1000);
    for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {
      Bi->X = 0;
      Bj->X = 1;
    }
    Cells cells;
    cells.resize(2);
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;
    Ci->X = 0;
    Ci->NDLEAF = 10;
    Ci->LEAF = ibodies.begin();
    Ci->M = 0;
    Ci->L = 0;
    Cj->X = 1;
    Cj->NDLEAF = 1000;
    Cj->LEAF = jbodies.begin();
    Cj->M = 0;
    startTimer("P2P kernel");
    P2P(Ci,Cj,mutual);
    timeP2P = stopTimer("P2P kernel") / 10000;
    startTimer("M2L kernel");
    for( int i=0; i!=1000; ++i ) M2L(Ci,Cj);
    timeM2L = stopTimer("M2L kernel") / 1000;
    startTimer("M2P kernel");
    for( int i=0; i!=100; ++i ) M2P(Ci,Cj,mutual);
    timeM2P = stopTimer("M2P kernel") / 1000;
  }

#if QUARK
  inline void traverseBranch(C_iter Ci, C_iter Cj, Quark *quark, bool mutual);
#endif

};

#if QUARK
inline void traverseQuark(Quark *quark) {
  Evaluator *E;
  C_iter Ci, Cj, Ci0, Cj0;
  bool mutual;
  quark_unpack_args_6(quark,E,Ci,Cj,Ci0,Cj0,mutual);
  ThreadTrace beginTrace;
  E->startTracer(beginTrace);
  E->traverseBranch(Ci,Cj,mutual);
  E->stopTracer(beginTrace,0x0000ff);
}

void Evaluator::traverseBranch(C_iter Ci, C_iter Cj, Quark *quark, bool mutual) {
  char string[256];
  sprintf(string,"%d %d",int(Ci-Ci0),int(Cj-Cj0));
  Quark_Task_Flags tflags = Quark_Task_Flags_Initializer;
  QUARK_Task_Flag_Set(&tflags,TASK_LABEL,intptr_t(string) );
  if( mutual ) {
    QUARK_Insert_Task(quark,traverseQuark,&tflags,
                      sizeof(Evaluator),this,NODEP,
                      sizeof(Cell),&*Ci,OUTPUT,
                      sizeof(Cell),&*Cj,OUTPUT,
                      sizeof(Cell),&*Ci0,NODEP,
                      sizeof(Cell),&*Cj0,NODEP,
                      sizeof(bool),&mutual,VALUE,
                      0);
  } else {
    QUARK_Insert_Task(quark,traverseQuark,&tflags,
                      sizeof(Evaluator),this,NODEP,
                      sizeof(Cell),&*Ci,OUTPUT,
                      sizeof(Cell),&*Cj,NODEP,
                      sizeof(Cell),&*Ci0,NODEP,
                      sizeof(Cell),&*Cj0,NODEP,
                      sizeof(bool),&mutual,VALUE,
                      0);
  }
}
#endif // QUARK
#undef splitFirst
#endif
