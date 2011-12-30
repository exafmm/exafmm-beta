#include "serialfmm.h"

int main() {
  const int numBodies = 1000;
  const real xmax = 100.0;
  const real ksize = 11.0;
  const real alpha = 0.1;
  const real sigma = .25 / M_PI;
  IMAGES = 2;
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells, jcells;
  SerialFMM<Laplace> FMM;
  FMM.printNow = true;

  srand48(2);
  // set positions and types
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

  FMM.setDomain(bodies,xmax/2,xmax/2);
  FMM.bottomup(bodies,cells);
  FMM.setEwald(ksize,alpha,sigma);
  jcells = cells;

  FMM.startTimer("Ewald        ");
  FMM.Ewald(bodies,cells,jcells);
  FMM.stopTimer("Ewald        ",FMM.printNow);
  FMM.eraseTimer("Ewald        ");

  FMM.startTimer("Direct       ");
  jbodies = FMM.periodicBodies(bodies);
  FMM.buffer = bodies;
  FMM.initTarget(FMM.buffer);
  FMM.evalP2P(FMM.buffer,jbodies);
  FMM.stopTimer("Direct       ",FMM.printNow);
  FMM.eraseTimer("Direct       ");

  // Calculate error
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  FMM.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2,true);
  FMM.printError(diff1,norm1,diff2,norm2);
}
