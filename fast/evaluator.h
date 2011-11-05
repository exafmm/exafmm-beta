#ifndef evaluator_h
#define evaluator_h
#if CART
#include "cartesian.h"
#elif SPHE
#include "spherical.h"
#endif

class Evaluator : public Kernel {
private:
  typedef std::pair<C_iter,C_iter> Pair;
  mutable std::stack<C_iter> selfStack;
  mutable std::stack<Pair> pairStack;

protected:
  int      MAXLEVEL;
  unsigned NCELL;
  Bodies   BODIES;// Deleteing this slows down the code.
  C_iter   ROOT;

private:
  real getBmax(vect const&X, C_iter C) const {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  void interact(C_iter C) const {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {
      P2P(C);
    } else {
      selfStack.push(C);
    }
  }

  void interact(C_iter Ci, C_iter Cj, bool mutual=true) const {
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

  bool split_first(C_iter Ci, C_iter Cj) const {
    return Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT > Cj->RCRIT);
  }

protected:
  void setCenter(C_iter C) const {
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += B->SRC[0];
      X += B->X * B->SRC[0];
    }
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      m += std::abs(c->M[0]);
      X += c->X * std::abs(c->M[0]);
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  void setRcrit() {
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(ROOT->M[0]),1.0/3);
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
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

  void traverse() const {
    interact(ROOT);
    while(!selfStack.empty()) {
      C_iter C = selfStack.top();
      selfStack.pop();
      for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
        interact(Ci);
        for( C_iter Cj=Ci+1; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          interact(Ci,Cj);
        }
      }
      while(!pairStack.empty()) {
        Pair Cij = pairStack.top();
        pairStack.pop();
        if(split_first(Cij.first,Cij.second)) {
          C = Cij.first;
          for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
            interact(Ci,Cij.second);
          }
        } else {
          C = Cij.second;
          for( C_iter Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
            interact(Cij.first,Cj);
          }
        }
      }
    }
  }

  void traverse(bool mutual) const {
    for( C_iter Cj=C0+ROOT->CHILD; Cj!=C0+ROOT->CHILD+ROOT->NCHILD; ++Cj ) {
      Pair pair(ROOT,Cj);
      pairStack.push(pair);
    }
    while(!pairStack.empty()) {
      Pair Cij = pairStack.top();
      pairStack.pop();
      if(split_first(Cij.first,Cij.second)) {
        C_iter C = Cij.first;
        for( C_iter Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
          interact(Ci,Cij.second,mutual);
        }
      } else {
        C_iter C = Cij.second;
        for( C_iter Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          interact(Cij.first,Cj,mutual);
        }
      }
    }
  }

public:
  Evaluator() {}
  ~Evaluator() {}
};

#endif
