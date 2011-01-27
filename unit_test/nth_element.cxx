#include "partition.h"
#include "dataset.h"
#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(1000000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies bodies2(numBodies);
  BottomUpTreeConstructor T(bodies);
  Dataset D(bodies);
  Partition mpi(bodies);
  toc = get_time();
  mpi.print("Allocate      : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  srand(mpi.rank()+1);
  D.random();
  toc = get_time();
  mpi.print("Set bodies    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setDomain();
  mpi.setDomain();
  toc = get_time();
  mpi.print("Set domain    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setIndex();
  mpi.binBodies(T.Ibody,0);
  toc = get_time();
  mpi.print("Set index     : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2);
  toc = get_time();
  mpi.print("Sort index    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  bigint nthGlobal = numBodies * mpi.size() / 3;
  bigint iSplit = mpi.nth_element(&T.Ibody[0],numBodies,nthGlobal);
  int nthLocal = mpi.splitBodies(T.Ibody,iSplit);
  toc = get_time();
  mpi.print("Nth element   : ",0);
  mpi.print(toc-tic);
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
