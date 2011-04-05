#include "dataset.h"
#include "evaluator.h"

int main() {
  const int numBodies = 1000;
  std::string kernelName = "Laplace";
  Bodies bodies(numBodies);
  Bodies jbodies;
  Dataset D;
  Evaluator E;
  E.setKernel(kernelName);
  E.initialize();
  D.kernelName = kernelName;
  E.preCalculation();
  E.printNow = true;

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.stopTimer("Set bodies   ",E.printNow);

  E.startTimer("Set domain   ");
  E.setDomain(bodies);
  E.stopTimer("Set domain   ",E.printNow);

  if( IMAGES != 0 ) {
    E.startTimer("Set periodic ");
    jbodies = E.periodicBodies(bodies);
    E.stopTimer("Set periodic ",E.printNow);
  } else {
    jbodies = bodies;
  }

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,jbodies);
  E.stopTimer("Direct GPU   ",E.printNow);

  E.startTimer("Direct CPU   ");
  bool onCPU = true;
  Bodies bodies2 = bodies;
  D.initTarget(bodies2);
  E.evalP2P(bodies2,jbodies,onCPU);
  E.stopTimer("Direct CPU   ",E.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  D.printError(diff1,norm1,diff2,norm2);
  E.postCalculation();
  E.finalize();
}
