#include "serialfmm.h"

int main() {
  const int numBodies = 1000;
  const real xmax = 100.0;
  const real ksize = 11.0;
  const real alpha = 0.1;
  const real sigma = .25 / M_PI;
  IMAGES = 3;
  THETA = 1 / sqrt(4);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells, jcells;
  SerialFMM<Laplace> FMM;
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  srand48(2);
  real average = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int d=0; d!=3; ++d ) {
      B->X[d] = drand48() * xmax;
    }
    B->SRC = drand48();
    average += B->SRC;
    B->TRG = 0;
  }
  average /= numBodies;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->SRC -= average;
  }
  FMM.stopTimer("Set bodies   ",FMM.printNow);
  FMM.eraseTimer("Set bodies   ");

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies,xmax/2,xmax/2);
  FMM.stopTimer("Set domain   ",FMM.printNow);
  FMM.eraseTimer("Set domain   ");

  FMM.bottomup(bodies,cells);
  FMM.setEwald(ksize,alpha,sigma);
  jcells = cells;
  FMM.Ewald(bodies,cells,jcells);

  Bodies bodies2 = bodies;
  FMM.initTarget(bodies);
  FMM.startTimer("Downward     ");
  FMM.downward(cells,jcells);
  FMM.stopTimer("Downward     ",FMM.printNow);
  FMM.eraseTimer("Downward     ");

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2,true);
  FMM.printError(diff1,norm1,diff2,norm2);
  FMM.finalize();
}
