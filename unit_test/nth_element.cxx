#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(1000000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies buffer(numBodies);
  Dataset D(bodies);
  Partition P(bodies);
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
  P.BottomUp::setIndex();
  P.binBodies(0);
  toc = get_time();
  if(print) std::cout << "Set index     : " << toc-tic << std::endl;

  tic = get_time();
  P.sort(bodies,buffer);
  toc = get_time();
  if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

  tic = get_time();
  bigint nthGlobal = numBodies * P.commSize() / 3;
  bigint iSplit = P.nth_element(bodies,nthGlobal);
  int nthLocal = P.splitBodies(iSplit);
  toc = get_time();
  if(print) std::cout << "Nth element   : " << toc-tic << std::endl;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->I = B-bodies.begin() > nthLocal;
  }

#ifdef VTK
  if( P.commRank() == 0 ) {
    int Ncell(0);
    vtkPlot vtk;
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
}
