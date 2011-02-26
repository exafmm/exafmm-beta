#include "types.h"
#include "evaluator.h"

int main() {
  const int numBodies = 100;
  Bodies bodiesI(numBodies);
  Bodies bodiesI2(numBodies);
  Bodies bodiesJ(numBodies);
  Cells  icells;
  Cells  jcells;
  Evaluator E;
  E.initialize();

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

    Cell cell;
    cell.NLEAF  = numBodies;
    cell.LEAF   = bodiesJ.begin();
    cell.X      = 0.5;
    cell.M      = 0;
    cell.PARENT = 1;
    jcells.push_back(cell);
    E.evalP2M(jcells);
    cell.X      = 1;
    cell.M      = 0;
    jcells.push_back(cell);
    E.evalM2M(jcells);
    jcells.erase(jcells.begin());
    cell.X      = -1 - dist;
    cell.L      = 0;
    icells.push_back(cell);
    E.addM2L(jcells.begin());
    E.evalM2L(icells);
    cell.NLEAF  = numBodies;
    cell.LEAF   = bodiesI.begin();
    cell.X      = -0.5 - dist;
    cell.L      = 0;
    cell.NCHILD = 0;
    cell.PARENT = 1;
    icells.insert(icells.begin(),cell);
    E.evalL2L(icells);
    icells.pop_back();
//    E.evalL2P(icells);
    E.addM2P(jcells.begin());
    E.evalM2P(icells);
    icells.clear();
    jcells.clear();

    bodiesI2 = bodiesI;
    for( B_iter B=bodiesI2.begin(); B!=bodiesI2.end(); ++B ) {
      B->pot = 0;
    }
    E.evalP2P(bodiesI2,bodiesJ);

    B_iter B  = bodiesI.begin();
    B_iter B2 = bodiesI2.begin();
    real err = 0, rel = 0;
    for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
      err += (B->pot - B2->pot) * (B->pot - B2->pot);
      rel += B2->pot * B2->pot;
    }
    std::cout << dist << " " << std::sqrt(err/rel) << std::endl;
  }
  E.finalize();
}
