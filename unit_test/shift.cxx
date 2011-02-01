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
  if( mpi.rank() % 2 == 0 ) {
    D.random();
  } else {
    bodies.resize(50000);
    D.sphere();
  }
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  mpi.setGlobDomain(T.setR0(),T.setX0());
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  if( mpi.rank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=mpi.size(); ++i ) {
    mpi.shiftBodies(T.Ibody,buffer);
    if( mpi.rank() == 0 )
      vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  }
  toc = get_time();
  if(print) std::cout << "Shift bodies  : " << toc-tic << std::endl;
  if( mpi.rank() == 0 )
    vtk.plot(Ncell);
#endif
}
