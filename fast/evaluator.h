#ifndef evaluator_h
#define evaluator_h
#include <kernel.h>
#include <stack>
#include <iomanip>

class Evaluator : public Kernel {
private:
  typedef Pair<Cell*,Cell*> pair;
  mutable std::stack<Cell*> selfStack;
  mutable Stack<Cell*,Cell*> pairStack;

protected:
  unsigned LEVEL;

private:
  void perform(Cell *C) const {
    if(C->NCHILD == 0 || C->NDLEAF < 64) {
      P2P(C);
    } else {
      selfStack.push(C);
    }
  }

  void perform(Cell *Ci, Cell *Cj, bool mutual=true) const {
    vect dX = Ci->X - Cj->X;
    real Rq = norm(dX);
    if(Rq > (Ci->RCRIT+Cj->RCRIT)*(Ci->RCRIT+Cj->RCRIT)) {
      M2L(Ci,Cj,mutual);
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {
      P2P(Ci,Cj,mutual);
    } else {
      pairStack.push(Ci,Cj);
    }
  }

  bool split_first(Cell *Ci, Cell *Cj) const
  {
    return Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->RCRIT > Cj->RCRIT);
  }

  void set_rcrit() {
    real c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(C0->M[0],1.0/3);
    for( Cell* C=C0; C!=C0+NCELLS; ++C ) {
      real a = c * pow(C->M[0],1.0/3);
      real x = 1.0 / THETA;
      for(int i=0; i<5; i++) {
        real f = x * x - 2 * x + 1 - a * pow(x,-P);
        real df = (P + 2) * x - 2 * (P + 1) + P / x;
        x -= f / df;
      }
      C->RCRIT *= x;
    }
  }

  void upward() {
    for( Cell* C=C0; C!=C0+NCELLS; ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( Cell *C=C0+NCELLS-1; C!=C0-1; --C ) {
      real bmax = 0;
      real dmax = 0;
      P2M(C,dmax,bmax);
      M2M(C,dmax,bmax);
    }
    for( Cell* C=C0; C!=C0+NCELLS; ++C ) {
      C->M[1] *= 0.5/C->M[0];
      C->M[2] *= 0.5/C->M[0];
      C->M[3] *= 0.5/C->M[0];
      C->M[4] *= 0.5/C->M[0];
      C->M[5] *= 0.5/C->M[0];
      C->M[6] *= 0.5/C->M[0];
    }
    set_rcrit();
  }

  void UpdateLeafs() {
    for( B_iter B=LEAFS.begin(); B!=LEAFS.end(); ++B ) {      // Loop over bodies
      B->SRC[0] = BODIES[B->IBODY].SRC[0];
      B->TRG = 0;
    }
  }

  void downward(Cell *C) const {
    L2L(C);
    L2P(C);
    for( Cell *c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      downward(c);
    }
  }

  void UpdateBodies() {
    for( B_iter B=LEAFS.begin(); B!=LEAFS.end(); ++B ) {      // Loop over bodies
      BODIES[B->IBODY].TRG[0] = B->TRG[0];
      BODIES[B->IBODY].TRG[1] = B->TRG[1];
      BODIES[B->IBODY].TRG[2] = B->TRG[2];
      BODIES[B->IBODY].TRG[3] = B->TRG[3];
    }
  }

  void write() const {
    std::cout<<" root center:           "<<C0->X            <<'\n';
    std::cout<<" root radius:           "<<RAD              <<'\n';
    std::cout<<" bodies loaded:         "<<C0->NDLEAF       <<'\n';
    std::cout<<" total scal:            "<<C0->M[0]         <<'\n';
    std::cout<<" cells used:            "<<NCELLS           <<'\n';
    std::cout<<" maximum level:         "<<LEVEL            <<'\n';
  }

protected:
  void traverse() const {
    perform(C0);
    while(!selfStack.empty()) {
      Cell *C = selfStack.top();
      selfStack.pop();
      for( Cell *Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
        perform(Ci);
        for( Cell *Cj=Ci+1; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          perform(Ci,Cj);
        }
      }
      while(!pairStack.empty()) {
        pair Cij = pairStack.pop();
        if(split_first(Cij.first,Cij.second)) {
          C = Cij.first;
          for( Cell *Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
            perform(Ci,Cij.second);
          }
        } else {
          C = Cij.second;
          for( Cell *Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
            perform(Cij.first,Cj);
          }
        }
      }
    }
  }

  void traverse(bool mutual) const {
    for( Cell *Cj=C0+C0->CHILD; Cj!=C0+C0->CHILD+C0->NCHILD; ++Cj ) {
      pairStack.push(C0,Cj);
    }
    while(!pairStack.empty()) {
      pair Cij = pairStack.pop();
      if(split_first(Cij.first,Cij.second)) {
        Cell *C = Cij.first;
        for( Cell *Ci=C0+C->CHILD; Ci!=C0+C->CHILD+C->NCHILD; ++Ci ) {
          perform(Ci,Cij.second,mutual);
        }
      } else {
        Cell *C = Cij.second;
        for( Cell *Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; ++Cj ) {
          perform(Cij.first,Cj,mutual);
        }
      }
    }
  }

public:
  Evaluator(Bodies &bodies, real rad, unsigned L, unsigned nleafs, unsigned ncells)
    : Kernel( bodies,rad,nleafs,ncells ), pairStack( 16*L-12 ), LEVEL ( L ) {}
  ~Evaluator() {}

  void exact() {
    UpdateLeafs();
    P2P(C0);
    for( B_iter B=LEAFS.begin(); B!=LEAFS.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
    UpdateBodies();
  }

  void approximate() {
    double tic,toc;
    tic = get_time();
    UpdateLeafs();
    upward();
    toc = get_time();
    std::cout << "upward : " << toc-tic << std::endl;
    tic = get_time();
    traverse();
    toc = get_time();
    std::cout << "intrct : " << toc-tic << std::endl;
    tic = get_time();
    for( Cell *C=C0+C0->CHILD; C!=C0+C0->CHILD+C0->NCHILD; ++C ) {
      downward(C);
    }
    toc = get_time();
    std::cout << "downwd : " << toc-tic << std::endl;
#ifndef MANY
    write();
#endif
    UpdateBodies();
  }
};

#endif
