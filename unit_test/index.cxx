#include "dataset.h"
#include "serialfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 1000000;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset dataset;
  dataset.kernelName = "Laplace";
  SerialFMM FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.sphere(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  FMM.startTimer("Set index    ");
  FMM.BottomUp::setIndex(bodies);
  FMM.stopTimer("Set index    ",FMM.printNow);

  FMM.buffer.resize(bodies.size());
  FMM.sortBodies(bodies,FMM.buffer);

  bigint oldIndex(bodies[0].ICELL);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    assert( oldIndex <= B->ICELL );
    oldIndex = B->ICELL;
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(FMM.getR0(),FMM.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  FMM.finalize();
}
