#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 1000000;
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
  D.random(bodies,MPIRANK+1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

  T.startTimer("Set index    ");
  T.BottomUp::setIndex(bodies);
  T.binBodies(bodies,0);
  T.stopTimer("Set index    ",T.printNow);

  T.buffer.resize(bodies.size());
  T.sortBodies(bodies,T.buffer);

  T.startTimer("Nth element  ");
  bigint nthGlobal = numBodies * MPISIZE / 3;
  bigint iSplit = T.nth_element(bodies,nthGlobal);
  int nthLocal = T.splitBodies(bodies,iSplit);
  T.stopTimer("Nth element  ",T.printNow);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->ICELL = B-bodies.begin() > nthLocal;
  }

#ifdef VTK
  if( MPIRANK == 0 ) {
    int Ncell = 0;
    vtkPlot vtk;
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
