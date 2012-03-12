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

class Evaluator : public Kernel<Laplace> {
private:
  real timeP2P;
  real timeM2P;
  real timeM2L;

protected:
  bool    TOPDOWN;

private:
  real getBmax(vect const&X, C_iter C) const {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  inline void approximate(C_iter Ci, C_iter Cj) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {
      P2P(Ci,Cj);
      NP2P++;
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {
      M2P(Ci,Cj);
      NM2P++;
    } else {
      M2L(Ci,Cj);
      NM2L++;
    }
#elif TREECODE
    M2P(Ci,Cj);
    NM2P++;
#else
    M2L(Ci,Cj);
    NM2L++;
#endif
  }

public:
  void interact(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    vect dX = Ci->X - Cj->X;
    real Rq = norm(dX);
    if(Rq > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
      approximate(Ci,Cj);
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      P2P(Ci,Cj);
      NP2P++;
    } else {
      Pair pair(Ci,Cj);
      pairQueue.push(pair);
    }
  }

#if QUARK
  inline void interact(C_iter Ci, C_iter Cj, Quark *quark);
#endif

protected:
  C_iter setRootCell(Cells &cells) {
    if( TOPDOWN ) {
      return cells.begin();
    } else {
      return cells.end() - 1;
    }
  }

  void setCenter(C_iter C) const {
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += std::abs(B->SRC);
      X += B->X * std::abs(B->SRC);
    }
    for( C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; ++c ) {
      m += std::abs(c->M[0]);
      X += c->X * std::abs(c->M[0]);
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  void setRcrit(Cells &cells) {
    C_iter root = setRootCell(cells);
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

  void traverse(Cells &icells, Cells &jcells) {
    Ci0 = icells.begin();
    Cj0 = jcells.begin();
    Pair pair;
    if( TOPDOWN ) {
      pair = make_pair(icells.begin(),jcells.begin());
    } else {
      pair = make_pair(icells.end()-1,jcells.end()-1);
    }
    PairQueue pairQueue;
    pairQueue.push(pair);
#if QUARK
    Quark *quark = QUARK_New(4);
#endif
    while( !pairQueue.empty() ) {
      Pair Cij = pairQueue.front();
      pairQueue.pop();
      if(splitFirst(Cij.first,Cij.second)) {
        C_iter C = Cij.first;
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
          interact(Ci,Cij.second,pairQueue);
        }
      } else {
        C_iter C = Cij.second;
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
          interact(Cij.first,Cj,pairQueue);
        }
      }
#if QUARK
      if( pairQueue.size() > 100 ) {
        while( !pairQueue.empty() ) {
          Cij = pairQueue.front();
          pairQueue.pop();
          interact(Cij.first,Cij.second,quark);
        }
      }
#endif
    }
#if QUARK
    QUARK_Delete(quark);
    writeTrace();
#endif
  }

public:
  Evaluator() {}
  ~Evaluator() {}

  void timeKernels() {
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
    P2P(Ci,Cj);
    timeP2P = stopTimer("P2P kernel") / 10000;
    startTimer("M2L kernel");
    for( int i=0; i!=1000; ++i ) M2L(Ci,Cj);
    timeM2L = stopTimer("M2L kernel") / 1000;
    startTimer("M2P kernel");
    for( int i=0; i!=100; ++i ) M2P(Ci,Cj);
    timeM2P = stopTimer("M2P kernel") / 1000;
  }

};

#if QUARK
inline void interactQuark(Quark *quark) {
  Evaluator *E;
  C_iter CI, CJ, Ci0, Cj0;
  quark_unpack_args_5(quark,E,CI,CJ,Ci0,Cj0);
  ThreadTrace beginTrace;
  E->startTracer(beginTrace);
  PairQueue privateQueue;
  Pair pair(CI,CJ);
  privateQueue.push(pair);
  while( !privateQueue.empty() ) {
    Pair Cij = privateQueue.front();
    privateQueue.pop();
    if(splitFirst(Cij.first,Cij.second)) {
      C_iter C = Cij.first;
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        E->interact(Ci,Cij.second,privateQueue);
      }
    } else {
      C_iter C = Cij.second;
      for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
        E->interact(Cij.first,Cj,privateQueue);
      }
    }
  }
  E->stopTracer(beginTrace,0x0000ff);
}

void Evaluator::interact(C_iter Ci, C_iter Cj, Quark *quark) {
  char string[256];
  sprintf(string,"%d %d",int(Ci-Ci0),int(Cj-Cj0));
  Quark_Task_Flags tflags = Quark_Task_Flags_Initializer;
  QUARK_Task_Flag_Set(&tflags,TASK_LABEL,intptr_t(string) );
  QUARK_Insert_Task(quark,interactQuark,&tflags,
                    sizeof(Evaluator),this,NODEP,
                    sizeof(Cell),&*Ci,OUTPUT,
                    sizeof(Cell),&*Cj,NODEP,
                    sizeof(Cell),&*Ci0,NODEP,
                    sizeof(Cell),&*Cj0,NODEP,
                    0);
}
#endif
#undef splitFirst
#endif
