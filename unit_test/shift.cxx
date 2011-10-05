#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset D;
  D.kernelName = "Laplace";
  Partition T;
  T.setKernel(D.kernelName);
  T.initialize();
  if( MPIRANK == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  if( MPIRANK % 2 == 0 ) {
    D.random(bodies,MPIRANK+1);
  } else {
    bodies.resize(50000);
    D.sphere(bodies,MPIRANK+1);
  }
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->ICELL = 0;

  int Ncell = 0;
  vtkPlot vtk;
  if( MPIRANK == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  T.startTimer("Shift bodies ");
  for( int i=1; i!=MPISIZE; ++i ) {
    T.shiftBodies(bodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  T.stopTimer("Shift bodies ",T.printNow);
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
