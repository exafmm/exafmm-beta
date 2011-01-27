#include "partition.h"
#include "dataset.h"
#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(100000);
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
  if( mpi.rank() % 2 == 0 ) {
    D.random();
  } else {
    bodies.resize(50000);
    D.sphere();
  }
  toc = get_time();
  mpi.print("Set bodies    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setDomain();
  mpi.setDomain();
  toc = get_time();
  mpi.print("Set domain    : ",0);
  mpi.print(toc-tic);

#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  if( mpi.rank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=mpi.size(); ++i ) {
    mpi.shiftBodies(bodies2);
    if( mpi.rank() == 0 )
      vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  }
  toc = get_time();
  mpi.print("Shift bodies  : ",0);
  mpi.print(toc-tic);
  if( mpi.rank() == 0 )
    vtk.plot(Ncell);
#endif
}
