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
    if(C->NCCELL == 0 || C->NDLEAF < 64) {
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
    } else if(Ci->NCCELL == 0 && Cj->NCCELL == 0) {
      P2P(Ci,Cj,mutual);
    } else {
      pairStack.push(Ci,Cj);
    }
  }

  bool split_first(Cell *Ci, Cell *Cj) const
  {
    return Cj->NCCELL == 0 || (Ci->NCCELL != 0 && Ci->RCRIT > Cj->RCRIT);
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
    for( Leaf* Li=L0; Li!=L0+NLEAFS; ++Li ) {
      Li->Q = BODIES[Li->I].SRC[0];
      Li->TRG = 0;
    }
  }

  void downward(Cell *C) const {
    L2L(C);
    L2P(C);
    for( Cell *c=C->FCCELL; c!=C->FCCELL+C->NCCELL; ++c ) {
      downward(c);
    }
  }

  void UpdateBodies() {
    for( Leaf* Li=L0; Li!=L0+NLEAFS; ++Li ) {
      BODIES[Li->I].TRG[0] = Li->TRG[0];
      BODIES[Li->I].TRG[1] = Li->TRG[1];
      BODIES[Li->I].TRG[2] = Li->TRG[2];
      BODIES[Li->I].TRG[3] = Li->TRG[3];
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
      for( Cell *Ci=C->FCCELL; Ci!=C->FCCELL+C->NCCELL; ++Ci ) {
        perform(Ci);
        for( Cell *Cj=Ci+1; Cj!=C->FCCELL+C->NCCELL; ++Cj ) {
          perform(Ci,Cj);
        }
      }
      while(!pairStack.empty()) {
        pair Cij = pairStack.pop();
        if(split_first(Cij.first,Cij.second)) {
          C = Cij.first;
          for( Cell *Ci=C->FCCELL; Ci!=C->FCCELL+C->NCCELL; ++Ci ) {
            perform(Ci,Cij.second);
          }
        } else {
          C = Cij.second;
          for( Cell *Cj=C->FCCELL; Cj!=C->FCCELL+C->NCCELL; ++Cj ) {
            perform(Cij.first,Cj);
          }
        }
      }
    }
  }

  void traverse(bool mutual) const {
    for( Cell *Cj=C0->FCCELL; Cj!=C0->FCCELL+C0->NCCELL; ++Cj ) {
      pairStack.push(C0,Cj);
    }
    while(!pairStack.empty()) {
      pair Cij = pairStack.pop();
      if(split_first(Cij.first,Cij.second)) {
        Cell *C = Cij.first;
        for( Cell *Ci=C->FCCELL; Ci!=C->FCCELL+C->NCCELL; ++Ci ) {
          perform(Ci,Cij.second,mutual);
        }
      } else {
        Cell *C = Cij.second;
        for( Cell *Cj=C->FCCELL; Cj!=C->FCCELL+C->NCCELL; ++Cj ) {
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
    for( Leaf* Li=L0; Li!=L0+NLEAFS; ++Li ) Li->TRG /= Li->Q;
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
    for( Cell *C=C0->FCCELL; C!=C0->FCCELL+C0->NCCELL; ++C ) {
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
