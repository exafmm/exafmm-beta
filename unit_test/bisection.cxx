#include "partition.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
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
  FMM.stopTimer("Set domain   ");

  FMM.startTimer("Set index    ");
  FMM.BottomUp::setIndex(bodies);
  FMM.stopTimer("Set index    ");

  FMM.buffer.resize(bodies.size());
  FMM.sortBodies(bodies,FMM.buffer);

  FMM.bisection(bodies);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->ICELL = 0;
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  if( MPIRANK == 0 ) {
    vtk.setDomain(FMM.getR0(),FMM.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    FMM.shiftBodies(bodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  FMM.finalize();
}
