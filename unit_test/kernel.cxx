#include "dataset.h"
#include "evaluator.h"

int main() {
  const int numBodies = 100;
  std::string kernelName = "Laplace";
  Bodies ibodies(numBodies);
  Bodies ibodies2(numBodies);
  Bodies jbodies(numBodies);
  Cells  icells;
  Cells  jcells;
  Dataset D;
  Evaluator E;
  E.setKernel(kernelName);
  E.initialize();
  D.kernelName = kernelName;
  E.preCalculation();

  for( int it=0; it!=10; ++it ) {
    real dist = (1 << it) / 2;
    for( B_iter B=ibodies.begin(); B!=ibodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->X[d] = -rand() / (1. + RAND_MAX) - dist;
      }
    }
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->X[d] = rand() / (1. + RAND_MAX);
      }
    }
    D.initSource(jbodies);
    bool IeqJ = false;
    D.initTarget(ibodies,IeqJ);

    Cell cell;
    cell.NLEAF    = numBodies;
    cell.LEAF     = jbodies.begin();
    cell.X        = 0.5;
    cell.M        = 0;
    cell.ICELL    = 8;
    cell.NCHILD   = 0;
    cell.PARENT   = 1;
    jcells.push_back(cell);
    E.evalP2M(jcells);
    cell.X        = 1;
    cell.M        = 0;
    cell.ICELL    = 0;
    cell.NCHILD   = 1;
    cell.CHILD[0] = 0;
    jcells.push_back(cell);
    E.evalM2M(jcells);
    jcells.erase(jcells.begin());
    cell.X        = -1 - dist;
    cell.L        = 0;
    cell.ICELL    = 0;
    icells.push_back(cell);
    E.addM2L(jcells.begin());
    E.evalM2L(icells);
    cell.NLEAF    = numBodies;
    cell.LEAF     = ibodies.begin();
    cell.X        = -0.5 - dist;
    cell.L        = 0;
    cell.ICELL    = 1;
    cell.NCHILD   = 0;
    cell.PARENT   = 1;
    icells.insert(icells.begin(),cell);
    E.evalL2L(icells);
    icells.pop_back();
    E.evalL2P(icells);
    E.addM2P(jcells.begin());
//    E.evalM2P(icells);
    icells.clear();
    jcells.clear();

    ibodies2 = ibodies;
    D.initTarget(ibodies2,IeqJ);
    E.evalP2P(ibodies2,jbodies);

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    D.evalError(ibodies,ibodies2,diff1,norm1,diff2,norm2);
    std::cout << "Distance      : " << dist << std::endl;
    D.printError(diff1,norm1,diff2,norm2);
  }
  E.postCalculation();
  E.finalize();
}
