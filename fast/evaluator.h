#ifndef evaluator_h
#define evaluator_h
#if CART
#include "kernel.h"
#elif SPHE
#include "kernel2.h"
#endif

class Evaluator : public Kernel {
private:
  typedef std::pair<C_iter,C_iter> Pair;
  mutable std::stack<C_iter> selfStack;
  mutable std::stack<Pair> pairStack;

protected:
  unsigned NCELL;
  Bodies BODIES;// Deleteing this slows down the code.
  C_iter CN;

private:
  void perform(C_iter C) const {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {
      P2P(C);
    } else {
      selfStack.push(C);
    }
  }

  void perform(C_iter Ci, C_iter Cj, bool mutual=true) const {
    vect dX = Ci->X - Cj->X;
    real Rq = norm(dX);
    if(Rq > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
      M2L(Ci,Cj,mutual);
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      P2P(Ci,Cj,mutual);
    } else {
      Pair pair(Ci,Cj);
      pairStack.push(pair);
    }
  }

  bool split_first(C_iter Ci, C_iter Cj) const
  {
    return Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT > Cj->RCRIT);
  }

protected:
  void set_rcrit() {
#if SPHE
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(CN->M[0].real(),1.0/3);
#else
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(CN->M[0],1.0/3);
#endif
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
#if SPHE
      real a = c * pow(C->M[0].real(),1.0/3);
#else
      real a = c * pow(C->M[0],1.0/3);
#endif
      real x = 1.0 / THETA;
      for( int i=0; i<5; ++i ) {
        real f = x * x - 2 * x + 1 - a * pow(x,-P);
        real df = (P + 2) * x - 2 * (P + 1) + P / x;
        x -= f / df;
      }
      C->RCRIT *= x;
    }
  }

  void upward() {
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=C0+NCELL-1; C!=C0-1; --C ) {
      setCenter(C);
      P2M(C);
      M2M(C);
    }
#if CART
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      for( int i=1; i<MCOEF; ++i ) C->M[i] /= C->M[0];
    }
#endif
    set_rcrit();
  }

  void downward(C_iter C) const {
    L2L(C);
    L2P(C);
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      downward(c);
    }
  }

  void traverse() const {
    perform(CN);
    while(!selfStack.empty()) {
      C_iter C = selfStack.top();
      selfStack.pop();
      for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
        perform(Ci);
        for( C_iter Cj=Ci+1; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          perform(Ci,Cj);
        }
      }
      while(!pairStack.empty()) {
        Pair Cij = pairStack.top();
        pairStack.pop();
        if(split_first(Cij.first,Cij.second)) {
          C = Cij.first;
          for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
            perform(Ci,Cij.second);
          }
        } else {
          C = Cij.second;
          for( C_iter Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
            perform(Cij.first,Cj);
          }
        }
      }
    }
  }

  void traverse(bool mutual) const {
    for( C_iter Cj=C0+CN->CHILD; Cj!=C0+CN->CHILD+CN->NCHILD; ++Cj ) {
      Pair pair(CN,Cj);
      pairStack.push(pair);
    }
    while(!pairStack.empty()) {
      Pair Cij = pairStack.top();
      pairStack.pop();
      if(split_first(Cij.first,Cij.second)) {
        C_iter C = Cij.first;
        for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
          perform(Ci,Cij.second,mutual);
        }
      } else {
        C_iter C = Cij.second;
        for( C_iter Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          perform(Cij.first,Cj,mutual);
        }
      }
    }
  }

public:
  Evaluator() {}
  ~Evaluator() {}
};

#endif
