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
  bool TOPDOWN;
  real NM2L;                                                    //!< Number of M2L kernel calls
  real NM2P;                                                    //!< Number of M2P kernel calls
  real NP2P;                                                    //!< Number of P2P kernel calls

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::writeTrace;                           //!< Write traces of all events
  using Kernel<equation>::Ci0;                                  //!< icells.begin()
  using Kernel<equation>::Cj0;                                  //!< jcells.begin()
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
      NP2P++;
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {
      evalM2P(Ci,Cj);
      NM2P++;
    } else {
      evalM2L(Ci,Cj);
      NM2L++;
    }
#elif TREECODE
    evalM2P(Ci,Cj);
    NM2P++;
#else
    evalM2L(Ci,Cj);
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
      evalP2P(Ci,Cj);
      NP2P++;
    } else {
      Pair pair(Ci,Cj);
      pairQueue.push(pair);
    }
  }

#if QUARK
  inline void interact(C_iter Ci, C_iter Cj, Quark *quark);
#endif

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

protected:
  C_iter setRootCell(Cells &cells) {
    if( TOPDOWN ) {
      return cells.begin();
    } else {
      return cells.end() - 1;
    }
  }

  void upwardPass(Cells &cells) {
    startTimer("Upward pass");
    evalP2M(cells);                                             // Evaluate all P2M kernels
    evalM2M(cells,cells);                                       // Evaluate all M2M kernels
#if Cartesian
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass",printNow);
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

  void timeKernels();

public:
  Evaluator() : NM2L(0), NM2P(0), NP2P(0) {}
  ~Evaluator() {}

  inline void evalP2M(Cells &cells);                            //!< Evaluate all P2M kernels
  inline void evalM2M(Cells &cells, Cells &jcells);             //!< Evaluate all M2M kernels
  void evalM2L(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalP2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalL2L(Cells &cells);                                   //!< Evaluate all L2L kernels
  void evalL2P(Cells &cells);                                   //!< Evaluate all L2P kernels

};

#include "CPUEvaluator.cxx"

#undef splitFirst
#endif
