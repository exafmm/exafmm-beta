#include "dataset.h"
#include "serialfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Bodies bodies2;
  Bodies jbodies;
  Cells cells,jcells;
  Dataset dataset;
  dataset.kernelName = "Laplace";
  SerialFMM FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.random(bodies,1,1);
  bodies2 = bodies;
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies2);
    FMM.stopTimer("Set periodic ",FMM.printNow);
  } else {
    jbodies = bodies2;
  }

  FMM.startTimer("Direct sum   ");
  bodies2.resize(numTarget);
  FMM.evalP2P(bodies2,jbodies);
  FMM.stopTimer("Direct sum   ",FMM.printNow);
  FMM.eraseTimer("Direct sum   ");

  FMM.startTimer("Set domain   ");
  dataset.initTarget(bodies);
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);
#else
  FMM.bottomup(bodies,cells);
#endif

  jcells = cells;
  FMM.startTimer("Downward     ");
  FMM.downward(cells,jcells,1);
  FMM.stopTimer("Downward     ",FMM.printNow);

  FMM.startTimer("Unsort bodies");
  std::sort(bodies.begin(),bodies.end());
  FMM.stopTimer("Unsort bodies",FMM.printNow);
  FMM.writeTime();
  FMM.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  bodies.resize(numTarget);
  dataset.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  dataset.printError(diff1,norm1,diff2,norm2);
  FMM.finalize();
}
