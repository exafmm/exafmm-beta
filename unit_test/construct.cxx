#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000000;
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  TreeConstructor T;

  T.startTimer("Set bodies   ");
  D.sphere(bodies);
  T.stopTimer("Set bodies   ",true);

  T.startTimer("Set domain   ");
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",true);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
