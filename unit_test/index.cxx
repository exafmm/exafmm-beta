#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 1000000;
  Bodies bodies(numBodies);
  Dataset D;
  TreeConstructor T;

  T.startTimer("Set bodies   ");
  D.sphere(bodies);
  T.stopTimer("Set bodies   ",true);

  T.startTimer("Set domain   ");
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",true);

  T.startTimer("Set index    ");
  T.BottomUp::setIndex(bodies);
  T.stopTimer("Set index    ",true);

  T.startTimer("Sort index   ");
  T.buffer.resize(bodies.size());
  T.sort(bodies,T.buffer);
  T.stopTimer("Sort index   ",true);

  bigint oldIndex(bodies[0].I);
  int b = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {
    assert( oldIndex <= B->I );
    oldIndex = B->I;
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
