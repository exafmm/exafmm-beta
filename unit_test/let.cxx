#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  const int numBodies = 100000;
  tic = get_time();
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  bool print = true;
  if( T.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.random(bodies,T.commRank()+1);
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setGlobDomain(bodies);
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  T.bisection(bodies);
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;

#ifdef TOPDOWN
  T.topdown(bodies,cells,print);
#else
  T.bottomup(bodies,cells,print);
#endif

  tic = get_time();
  T.commBodies(cells);
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  T.commCells(bodies,cells);
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
    Body body;
    body.I = 1;
    body.pos  = C->X;
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
