#include "dataset.h"
#include "serialfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000000;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Cells cells;
  Dataset dataset;
  dataset.kernelName = "Laplace";
  SerialFMM FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.random(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);
#else
  FMM.bottomup(bodies,cells);
#endif

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(FMM.getR0(),FMM.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  FMM.finalize();
}
