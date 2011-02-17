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
  Dataset D;
  Partition P;
  bool print(true);
  if( P.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.random(bodies,P.commRank()+1);
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  P.setGlobDomain(bodies);
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  P.BottomUp::setIndex(bodies);
  toc = get_time();
  if(print) std::cout << "Set index     : " << toc-tic << std::endl;

  tic = get_time();
  P.buffer.resize(bodies.size());
  P.sort(bodies,P.buffer);
  toc = get_time();
  if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

  tic = get_time();
  P.bisection(bodies);
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->I = 0;
  }

#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  if( P.commRank() == 0 ) {
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=P.commSize(); ++i ) {
    P.shiftBodies(bodies);
    if( P.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( P.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
