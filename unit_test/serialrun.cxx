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
  FMM.stopTimer("Set bodies   ",FMM.printNow);
  FMM.eraseTimer("Set bodies   ");

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);
  FMM.eraseTimer("Set domain   ");

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);
#else
  FMM.bottomup(bodies,cells);
#endif
  jcells = cells;
  FMM.startTimer("Downward     ");
  FMM.downward(cells,jcells,1);
  FMM.stopTimer("Downward     ",FMM.printNow);
  FMM.eraseTimer("Downward     ");

  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies);
    FMM.stopTimer("Set periodic ",FMM.printNow);
    FMM.eraseTimer("Set periodic ");
  } else {
    jbodies = bodies;
  }

  FMM.startTimer("Direct sum   ");
  bodies.resize(numTarget);
  FMM.buffer = bodies;
  dataset.initTarget(FMM.buffer);
  FMM.evalP2P(FMM.buffer,jbodies);
  FMM.stopTimer("Direct sum   ",FMM.printNow);
  FMM.eraseTimer("Direct sum   ");
  FMM.writeTime();
  FMM.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  dataset.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2);
  dataset.printError(diff1,norm1,diff2,norm2);
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(FMM.getR0(),FMM.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  FMM.finalize();
}
