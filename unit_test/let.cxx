#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  if( T.commRank() == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,T.commRank()+1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

  T.bisection(bodies);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  T.commBodies(cells);

  T.commCells(bodies,cells);

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
    Body body;
    body.I = 1;
    body.X  = C->X;
    body.scal = 0;
    bodies.push_back(body);
  }

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
