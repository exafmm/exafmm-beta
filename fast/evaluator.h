#ifndef evaluator_h
#define evaluator_h
#include "cartesian.h"
#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT >= Cj->RCRIT)

class Evaluator : public Kernel {
private:
  real timeP2P;
  real timeM2P;
  real timeM2L;

protected:
  bool    TOPDOWN;
  C_iter  ROOT, ROOT2;

public:
  real NP2P;
  real NM2P;
  real NM2L;

private:
  real getBmax(vect const&X, C_iter C) const {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  inline void approximate(C_iter Ci, C_iter Cj, bool mutual=true) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {
      P2P(Ci,Cj,mutual);
      NP2P++;
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {
      M2P(Ci,Cj,mutual);
      NM2P++;
    } else {
      M2L(Ci,Cj,mutual);
      NM2L++;
    }
#elif TREECODE
    M2P(Ci,Cj,mutual);
    NM2P++;
#else
    M2L(Ci,Cj,mutual);
    NM2L++;
#endif
  }

  void traverseBranch(C_iter Ci, C_iter Cj, bool mutual=true) {
    PairQueue privateQueue;
    Pair pair(Ci,Cj);
    privateQueue.push_back(pair);
    while( !privateQueue.empty() ) {
      pair = privateQueue.front();
      privateQueue.pop_front();
      if(splitFirst(pair.first,pair.second)) {
        C_iter C = pair.first;
        for( C_iter CC=Ci0+C->CHILD; CC!=Ci0+C->CHILD+C->NCHILD; ++CC ) {
          applyMAC(CC,pair.second,privateQueue,mutual);
        }
      } else {
        C_iter C = pair.second;
        for( C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; ++CC ) {
          applyMAC(pair.first,CC,privateQueue,mutual);
        }
      }
    }
  }

protected:
  void setRootCell(Cells &cells) {
    Ci0 = cells.begin();
    Cj0 = cells.begin();
    if( TOPDOWN ) {
      ROOT = Ci0;
    } else {
      ROOT = cells.end() - 1;
    }
  }

  void setRootCell(Cells &icells, Cells &jcells) {
    Ci0 = icells.begin();
    Cj0 = jcells.begin();
    if( TOPDOWN ) {
      ROOT  = Ci0;
      ROOT2 = Cj0;
    } else {
      ROOT  = icells.end() - 1;
      ROOT2 = jcells.end() - 1;
    }
  }

  void setCenter(C_iter C) const {
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += B->SRC;
      X += B->X * B->SRC;
    }
    for( C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; ++c ) {
      m += std::abs(c->M[0]);
      X += c->X * std::abs(c->M[0]);
    }
    X /= m;
#if USE_BMAX
    C->R = getBmax(X,C);
#endif
#if COMcenter
    C->X = X;
#endif
  }

  void setRcrit(Cells &cells) {
#if ERROR_OPT
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(ROOT->M[0]),1.0/3);
#endif
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      real x = 1.0 / THETA;
#if ERROR_OPT
      real a = c * pow(std::abs(C->M[0]),1.0/3);
      for( int i=0; i<5; ++i ) {
        real f = x * x - 2 * x + 1 - a * pow(x,-P);
        real df = (P + 2) * x - 2 * (P + 1) + P / x;
        x -= f / df;
      }
#endif
      C->RCRIT *= x;
    }
  }

  void pushCell(C_iter C, CellQueue &cellQueue) {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {
      P2P(C);
      NP2P++;
    } else {
      cellQueue.push(C);
    }
  }

  void traverse(CellQueue &cellQueue) {
    PairQueue pairQueue;
    while( !cellQueue.empty() ) {
      C_iter C = cellQueue.front();
      cellQueue.pop();
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        pushCell(Ci,cellQueue);
        for( C_iter Cj=Ci+1; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
          applyMAC(Ci,Cj,pairQueue);
        }
      }
      traverse(pairQueue,true);
    }
  }

  void traverse(PairQueue &pairQueue, bool mutual=false) {
#if MTHREADS
    Pair pair = pairQueue.front();
    pairQueue.pop_back();
    recursiveTraverse(pair.first,pair.second,mutual);
#else
#if QUARK
    Quark *quark = QUARK_New(12);
#endif
    while( !pairQueue.empty() ) {
      Pair pair = pairQueue.front();
      pairQueue.pop_front();
      if(splitFirst(pair.first,pair.second)) {
        C_iter C = pair.first;
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
          applyMAC(Ci,pair.second,pairQueue,mutual);
        }
      } else {
        C_iter C = pair.second;
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
          applyMAC(pair.first,Cj,pairQueue,mutual);
        }
      }
#if QUARK
      if( int(pairQueue.size()) > ROOT->NDLEAF / 100 ) {
        while( !pairQueue.empty() ) {
          pair = pairQueue.front();
          pairQueue.pop_front();
          traverseBranch(pair.first,pair.second,quark,mutual);
        }
      }
#else
      if( int(pairQueue.size()) > ROOT->NDLEAF / 100 ) {
//#pragma omp parallel for schedule(dynamic)
        for( int i=0; i<int(pairQueue.size()); i++ ) {
          traverseBranch(pairQueue[i].first,pairQueue[i].second,mutual);
        }
        pairQueue.clear();
      }
#endif // QUARK
    }
#if QUARK
    QUARK_Delete(quark);
    writeTrace();
#endif
#endif // MTHREADS
  }

#if MTHREADS
  void recursiveTraverse(C_iter Ci, C_iter Cj, bool mutual) {
    vect dX = Ci->X - Cj->X - Xperiodic;
    real R2 = norm(dX);
#if DUAL
    {
#else
    if(Ci->RCRIT != Cj->RCRIT) {
      if(splitFirst(Ci,Cj)) {
        C_iter C = Ci;
        task_group tg;
        for( C_iter CC=Ci0+C->CHILD; CC!=Ci0+C->CHILD+C->NCHILD; ++CC ) {
          if( CC->NDLEAF > 10000 ) tg.run([=]{recursiveTraverse(CC,Cj,mutual);});
          else recursiveTraverse(CC,Cj,mutual);
        }
        tg.wait();
      } else {
        C_iter C = Cj;
        for( C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; ++CC ) {
          recursiveTraverse(Ci,CC,mutual);
        }
      }
    } else {
#endif
      if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
        approximate(Ci,Cj,mutual);
      } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
        P2P(Ci,Cj,mutual);
        NP2P++;
      } else {
        if(splitFirst(Ci,Cj)) {
          C_iter C = Ci;
          task_group tg;
          for( C_iter CC=Ci0+C->CHILD; CC!=Ci0+C->CHILD+C->NCHILD; ++CC ) {
            if( CC->NDLEAF > 10000 ) tg.run([=]{recursiveTraverse(CC,Cj,mutual);});
            else recursiveTraverse(CC,Cj,mutual);
          }
          tg.wait();
        } else {
          C_iter C = Cj;
          for( C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; ++CC ) {
            recursiveTraverse(Ci,CC,mutual);
          }
        }
      }
    }
  }
#endif // MTHREADS

//! Traverse tree for periodic cells
  void traversePeriodic() {
    startTimer("Traverse periodic");                            // Start timer
    Xperiodic = 0;                                              // Periodic coordinate offset
    real R = R0;                                                // Radius at current level
    Cells pcells(28);                                           // Create cells
    C_iter Ci = pcells.end()-1;                                 // Last cell is periodic parent cell
    *Ci = *ROOT2;                                               // Copy values from source root
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
                    Xperiodic[1] = (iy * 3 + cy) * 2 * R;       //         Coordinate offset for y periodic directio
                    Xperiodic[2] = (iz * 3 + cz) * 2 * R;       //         Coordinate offset for y periodic directio
                    M2L(ROOT,Ci,false);                         //         Perform M2L kernel
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
      real Rmax = 0;                                            //  Dummy parameter for calling M2M
      setCenter(Ci);                                            //  Set center of mass for periodic parent
      M2M(Ci,Rmax);                                             //  Evaluate periodic M2M kernels for this sublevel
      R *= 3;                                                   //  Increase center cell size three times
      Cj0 = C0;                                                 //  Reset Cj0 back
    }                                                           // End loop over sublevels of tree
    stopTimer("Traverse periodic",printNow);                    // Stop timer
  }

public:
  Evaluator() : TOPDOWN(false), NP2P(0), NM2P(0), NM2L(0) {}
  ~Evaluator() {}

  void applyMAC(C_iter Ci, C_iter Cj, PairQueue &pairQueue, bool mutual=true) {
    vect dX = Ci->X - Cj->X - Xperiodic;
    real R2 = norm(dX);
#if DUAL
    {
#else
    if(Ci->RCRIT != Cj->RCRIT) {
      Pair pair(Ci,Cj);
      pairQueue.push_back(pair);
    } else {
#endif
      if(R2 > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
        approximate(Ci,Cj,mutual);
      } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
        P2P(Ci,Cj,mutual);
        NP2P++;
      } else {
        Pair pair(Ci,Cj);
        pairQueue.push_back(pair);
      }
    }
  }

#if QUARK
  inline void traverseBranch(C_iter Ci, C_iter Cj, Quark *quark, bool mutual=true);
#endif

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

};

#if QUARK
inline void traverseQuark(Quark *quark) {
  Evaluator *E;
  C_iter CI, CJ, Ci0, Cj0;
  bool mutual;
  quark_unpack_args_6(quark,E,CI,CJ,Ci0,Cj0,mutual);
  ThreadTrace beginTrace;
  E->startTracer(beginTrace);
  PairQueue privateQueue;
  Pair pair(CI,CJ);
  privateQueue.push_back(pair);
  while( !privateQueue.empty() ) {
    pair = privateQueue.front();
    privateQueue.pop_front();
    if(splitFirst(pair.first,pair.second)) {
      C_iter C = pair.first;
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        E->applyMAC(Ci,pair.second,privateQueue,mutual);
      }
    } else {
      C_iter C = pair.second;
      for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
        E->applyMAC(pair.first,Cj,privateQueue,mutual);
      }
    }
  }
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
