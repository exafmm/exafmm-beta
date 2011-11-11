#include "partition.h"
#include "dataset.h"
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
  Partition FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  if( MPIRANK == 0 ) FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.random(bodies,MPIRANK+1);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setGlobDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  FMM.BottomUp::setIndex(bodies);
  FMM.binBodies(bodies,0);

  FMM.buffer.resize(bodies.size());
  FMM.sortBodies(bodies,FMM.buffer);

  FMM.startTimer("Nth element  ");
  bigint nthGlobal = numBodies * MPISIZE / 3;
  bigint iSplit = FMM.nth_element(bodies,nthGlobal);
  int nthLocal = FMM.splitBodies(bodies,iSplit);
  FMM.stopTimer("Nth element  ",FMM.printNow);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->ICELL = B-bodies.begin() > nthLocal;
  }

#ifdef VTK
  if( MPIRANK == 0 ) {
    int Ncell = 0;
    vtkPlot vtk;
    vtk.setDomain(FMM.getR0(),FMM.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
  FMM.finalize();
}
