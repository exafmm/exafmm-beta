#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  Bodies bodies(numBodies);
  Dataset D;
  Partition T;
  bool print = false;
  if( T.commRank() == 0 ) print = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,T.commRank()+1);
  T.stopTimer("Set bodies   ",print);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",print);

  T.startTimer("Set index    ");
  T.BottomUp::setIndex(bodies);
  T.stopTimer("Set index    ",print);

  T.startTimer("Sort index   ");
  T.buffer.resize(bodies.size());
  T.sort(bodies,T.buffer);
  T.stopTimer("Sort index   ",print);

  T.startTimer("Partition    ");
  T.bisection(bodies);
  T.stopTimer("Partition    ",print);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->I = 0;
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
}
