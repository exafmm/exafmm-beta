#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells,jcells;
  Dataset D;
  TreeConstructor T;
  T.printNow = true;

  T.startTimer("Set bodies   ");
  D.sphere(bodies,1,1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  jcells = cells;
  T.startTimer("Downward     ");
  T.downward(cells,jcells,1);
  T.stopTimer("Downward     ",T.printNow);

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies = T.periodicBodies(bodies);
    T.stopTimer("Set periodic ",T.printNow);
  } else {
    jbodies = bodies;
  }

  T.startTimer("Direct sum   ");
  Evaluator E;
  T.buffer = bodies;
  D.initTarget(T.buffer);
  T.evalP2P(T.buffer,jbodies);
  T.stopTimer("Direct sum   ",T.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  D.evalError(bodies,T.buffer,diff1,norm1,diff2,norm2);
  D.printError(diff1,norm1,diff2,norm2);
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
