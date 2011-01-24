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
  Partition mpi;
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
  toc = get_time();
  mpi.print("Set domain    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setMorton();
  toc = get_time();
  mpi.print("Set Morton    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2);
  toc = get_time();
  mpi.print("Sort Morton   : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  bigint median = numBodies * mpi.size() / 2;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B )
    T.Ibody[B-bodies.begin()] = mpi.rank()*numBodies/2 + (B-bodies.begin());
  median = mpi.nth_element(T.Ibody,numBodies,median);
  toc = get_time();
  mpi.print("Nth element   : ",0);
  mpi.print(toc-tic);
  mpi.print(median);
  mpi.print(T.Ibody,numBodies-10,numBodies);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B )
    T.Ibody[B-bodies.begin()] = T.Ibody[B-bodies.begin()] > median;

#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
