#include "kernel.h"

int main() {
  int const numBodies(100);
  Bodies bodiesI(numBodies);
  Bodies bodiesI2(numBodies);
  Bodies bodiesJ(numBodies);
  Cells  cells;
  cells.resize(4);
  Kernel K;
  K.initialize();

  for( int it=0; it!=10; ++it ) {
    real dist = (1 << it) / 2;
    for( B_iter B=bodiesI.begin(); B!=bodiesI.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->pos[d] = -rand() / (1. + RAND_MAX) - dist;
      }
      B->pot = 0;
    }
    for( B_iter B=bodiesJ.begin(); B!=bodiesJ.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->pos[d] = rand() / (1. + RAND_MAX);
      }
      B->scal = 1.0 / bodiesJ.size();
    }

    C_iter C = cells.begin();
    C->NLEAF = numBodies;
    C->LEAF = bodiesJ.begin();
    C->X = 0.5;
    C->M = 0;
    C++;
    C->X = 1;
    C->M = 0;
    C++;
    C->X = -1 - dist;
    C->L = 0;
    C++;
    C->NLEAF = numBodies;
    C->LEAF = bodiesI.begin();
    C->X = -0.5 - dist;
    C->L = 0;

    C = cells.begin();
    K.P2M(C);

    K.M2M(C+1,C);

    K.M2L(C+2,C+1);

    K.L2L(C+3,C+2);

    K.L2P(C+3);

//    K.M2P(C+3,C+1);

    bodiesI2 = bodiesI;
    for( B_iter B=bodiesI2.begin(); B!=bodiesI2.end(); ++B ) {
      B->pot = 0;
    }
    K.P2P(bodiesI2.begin(),bodiesI2.end(),bodiesJ.begin(),bodiesJ.end());

    B_iter B  = bodiesI.begin();
    B_iter B2 = bodiesI2.begin();
    real err(0), rel(0);
    for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
      err += (B->pot - B2->pot) * (B->pot - B2->pot);
      rel += B2->pot * B2->pot;
    }
    std::cout << dist << " " << std::sqrt(err/rel) << std::endl;
  }
  K.finalize();
}
