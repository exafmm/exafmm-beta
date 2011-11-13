#include "evaluator.h"

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Evaluator<Laplace> FMM;
  FMM.initialize();
  FMM.preCalculation();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  FMM.sphere(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies);
    FMM.stopTimer("Set periodic ",FMM.printNow);
  } else {
    jbodies = bodies;
  }

  FMM.startTimer("Direct GPU   ");
  FMM.evalP2P(bodies,jbodies);
  FMM.stopTimer("Direct GPU   ",FMM.printNow);

  FMM.startTimer("Direct CPU   ");
  bool onCPU = true;
  bodies.resize(numTarget);
  Bodies bodies2 = bodies;
  FMM.initTarget(bodies2);
  FMM.evalP2P(bodies2,jbodies,onCPU);
  FMM.stopTimer("Direct CPU   ",FMM.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  FMM.printError(diff1,norm1,diff2,norm2);
  FMM.postCalculation();
  FMM.finalize();
}
