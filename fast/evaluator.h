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
#if Cartesian
#include "cartesian.h"
#elif Spherical
#include "spherical.h"
#endif

class Evaluator : public Kernel {
private:
  typedef std::pair<C_iter,C_iter> Pair;
  mutable std::stack<C_iter> selfStack;
  mutable std::stack<Pair> pairStack;
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

  void interact(C_iter C) {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {
      P2P(C);
      NP2P++;
    } else {
      selfStack.push(C);
    }
  }

  void interact(C_iter Ci, C_iter Cj, bool mutual=true) {
    vect dX = Ci->X - Cj->X;
    real Rq = norm(dX);
    if(Rq > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
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
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      P2P(Ci,Cj,mutual);
      NP2P++;
    } else {
      Pair pair(Ci,Cj);
      pairStack.push(pair);
    }
  }

  bool splitFirst(C_iter Ci, C_iter Cj) const {
    return Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT > Cj->RCRIT);
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
      m += B->SRC[0];
      X += B->X * B->SRC[0];
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
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(ROOT->M[0]),1.0/3);
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

  void traverse() {
    interact(ROOT);
    while(!selfStack.empty()) {
      C_iter C = selfStack.top();
      selfStack.pop();
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        interact(Ci);
        for( C_iter Cj=Ci+1; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
          interact(Ci,Cj);
        }
      }
      while(!pairStack.empty()) {
        Pair Cij = pairStack.top();
        pairStack.pop();
        if(splitFirst(Cij.first,Cij.second)) {
          C = Cij.first;
          for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
            interact(Ci,Cij.second);
          }
        } else {
          C = Cij.second;
          for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
            interact(Cij.first,Cj);
          }
        }
      }
    }
  }

  void traverse(bool mutual) {
    for( C_iter Cj=Cj0+ROOT->CHILD; Cj!=Cj0+ROOT->CHILD+ROOT->NCHILD; ++Cj ) {
      Pair pair(ROOT,Cj);
      pairStack.push(pair);
    }
    while(!pairStack.empty()) {
      Pair Cij = pairStack.top();
      pairStack.pop();
      if(splitFirst(Cij.first,Cij.second)) {
        C_iter C = Cij.first;
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
          interact(Ci,Cij.second,mutual);
        }
      } else {
        C_iter C = Cij.second;
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
          interact(Cij.first,Cj,mutual);
        }
      }
    }
  }

public:
  Evaluator() : NP2P(0), NM2L(0) {}
  ~Evaluator() {}

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
    startTimer("P2P kernel   ");
    P2P(Ci,Cj,mutual);
    timeP2P = stopTimer("P2P kernel   ") / 10000;
    startTimer("M2L kernel   ");
    for( int i=0; i!=1000; ++i ) M2L(Ci,Cj);
    timeM2L = stopTimer("M2L kernel   ") / 1000;
    startTimer("M2P kernel   ");
    for( int i=0; i!=100; ++i ) M2P(Ci,Cj,mutual);
    timeM2P = stopTimer("M2P kernel   ") / 1000;
  }

};

#endif
