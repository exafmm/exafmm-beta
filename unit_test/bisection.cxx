#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  std::string kernelName = "Laplace";
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset D;
  Partition T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  if( MPIRANK == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,MPIRANK+1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ");

  T.startTimer("Set index    ");
  T.BottomUp::setIndex(bodies);
  T.stopTimer("Set index    ");

  T.startTimer("Sort bodies  ");
  T.buffer.resize(bodies.size());
  T.sortBodies(bodies,T.buffer);
  T.stopTimer("Sort bodies  ",T.printNow);

  T.bisection(bodies);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->ICELL = 0;
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  if( MPIRANK == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    T.shiftBodies(bodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
