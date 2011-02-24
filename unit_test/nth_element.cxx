#include "partition.h"
#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  const int numBodies = 1000000;
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D;
  Partition T;
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
  T.BottomUp::setIndex(bodies);
  T.binBodies(bodies,0);
  toc = get_time();
  if(print) std::cout << "Set index     : " << toc-tic << std::endl;

  tic = get_time();
  T.buffer.resize(bodies.size());
  T.sort(bodies,T.buffer);
  toc = get_time();
  if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

  tic = get_time();
  bigint nthGlobal = numBodies * T.commSize() / 3;
  bigint iSplit = T.nth_element(bodies,nthGlobal);
  int nthLocal = T.splitBodies(bodies,iSplit);
  toc = get_time();
  if(print) std::cout << "Nth element   : " << toc-tic << std::endl;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->I = B-bodies.begin() > nthLocal;
  }

#ifdef VTK
  if( T.commRank() == 0 ) {
    int Ncell = 0;
    vtkPlot vtk;
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
}
