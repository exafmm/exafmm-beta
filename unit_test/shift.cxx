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
  if( T.commRank() % 2 == 0 ) {
    D.random(bodies,T.commRank()+1);
  } else {
    bodies.resize(50000);
    D.sphere(bodies,T.commRank()+1);
  }
  T.stopTimer("Set bodies   ",print);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",print);

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;

  int Ncell = 0;
  vtkPlot vtk;
  if( T.commRank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  T.startTimer("Shift bodies ");
  for( int i=1; i!=T.commSize(); ++i ) {
    T.shiftBodies(bodies);
    if( T.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  T.stopTimer("Shift bodies ",print);
  if( T.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
