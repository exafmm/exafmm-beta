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
  TreeConstructor T(bodies);
  Dataset D(bodies);
  Partition mpi(bodies);
  bool print(true);
  if( mpi.rank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(mpi.rank()+1);
  D.random();
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  mpi.setGlobDomain(T.setR0(),T.setX0());
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  T.BottomUp::setIndex();
  mpi.binBodies(T.Ibody,0);
  toc = get_time();
  if(print) std::cout << "Set index     : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(T.Ibody,bodies,buffer);
  toc = get_time();
  if(print) std::cout << "Sort index    : " << toc-tic << std::endl;

  tic = get_time();
  bigint nthGlobal = numBodies * mpi.size() / 3;
  bigint iSplit = mpi.nth_element(&T.Ibody[0],numBodies,nthGlobal);
  int nthLocal = mpi.splitBodies(T.Ibody,iSplit);
  toc = get_time();
  if(print) std::cout << "Nth element   : " << toc-tic << std::endl;
  mpi.print(iSplit);
  mpi.print(&T.Ibody[0],numBodies-10,numBodies);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B )
    T.Ibody[B-bodies.begin()] = B-bodies.begin() > nthLocal;

#ifdef VTK
  if( mpi.rank() == 0 ) {
    int Ncell(0);
    vtkPlot vtk;
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
}
