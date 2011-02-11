#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(100000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies buffer(numBodies);
  Dataset D(bodies);
  LocalEssentialTree P(bodies);
  bool print(true);
  if( P.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(P.commRank()+1);
  D.random();
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  P.setGlobDomain();
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  P.bisection(buffer);
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;

#ifdef TOPDOWN
  P.topdown(buffer,print);
#else
  P.bottomup(buffer,print);
#endif

  tic = get_time();
  P.commBodies();
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  P.commCells(buffer);
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( B_iter B=buffer.begin(); B!=buffer.end(); ++B ) B->I = 1;
  bodies.insert(bodies.end(),buffer.begin(),buffer.end());

  int Ncell(0);
  vtkPlot vtk;
  if( P.commRank() == 0 ) {
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=P.commSize(); ++i ) {
    P.shiftBodies(buffer);
    if( P.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  toc = get_time();
  if(print) std::cout << "Shift bodies  : " << toc-tic << std::endl;
  if( P.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
