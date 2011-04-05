#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  std::string kernelName = "Laplace";
  Bodies bodies(numBodies);
  Dataset D;
  Partition T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  if( T.commRank() == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,T.commRank()+1);
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
  if( T.commRank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=T.commSize(); ++i ) {
    T.shiftBodies(bodies);
    if( T.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( T.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
