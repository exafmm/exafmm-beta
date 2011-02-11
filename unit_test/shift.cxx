#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(100000);
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D(bodies);
  Partition P(bodies);
  bool print(true);
  if( P.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(P.commRank()+1);
  if( P.commRank() % 2 == 0 ) {
    D.random();
  } else {
    bodies.resize(50000);
    D.sphere();
  }
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  P.setGlobDomain();
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;

  int Ncell(0);
  vtkPlot vtk;
  if( P.commRank() == 0 ) {
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=P.commSize(); ++i ) {
    P.shiftBodies();
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
