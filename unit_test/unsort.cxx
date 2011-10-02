#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  std::string kernelName = "Laplace";
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Bodies bodies2;
  Bodies jbodies;
  Cells cells,jcells;
  Dataset D;
  TreeConstructor T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,1,1);
  bodies2 = bodies;
  T.stopTimer("Set bodies   ",T.printNow);

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies = T.periodicBodies(bodies2);
    T.stopTimer("Set periodic ",T.printNow);
  } else {
    jbodies = bodies2;
  }

  T.startTimer("Direct sum   ");
  bodies2.resize(numTarget);
  T.evalP2P(bodies2,jbodies);
  T.stopTimer("Direct sum   ",T.printNow);
  T.eraseTimer("Direct sum   ");

  T.startTimer("Set domain   ");
  D.initTarget(bodies);
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  jcells = cells;
  T.startTimer("Downward     ");
  T.downward(cells,jcells,1);
  T.stopTimer("Downward     ",T.printNow);

  T.startTimer("Unsort bodies");
  std::sort(bodies.begin(),bodies.end());
  T.stopTimer("Unsort bodies",T.printNow);
  T.writeTime();
  T.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  bodies.resize(numTarget);
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  D.printError(diff1,norm1,diff2,norm2);
  T.finalize();
}
