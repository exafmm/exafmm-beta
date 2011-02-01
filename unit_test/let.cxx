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
  Bodies bodies2(numBodies);
  Dataset D(bodies);
  LocalEssentialTree mpi(bodies);
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
  mpi.setGlobDomain(mpi.setR0(),mpi.setX0());
  toc = get_time();
  mpi.print("Set domain    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.multisection(mpi.Ibody,bodies2);
  toc = get_time();
  mpi.print("Multisection  : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.BottomUpTreeConstructor::setIndex();
  toc = get_time();
  mpi.print("Set index     : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.sort(mpi.Ibody,bodies,bodies2);
  toc = get_time();
  mpi.print("Sort index    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.prune();
  toc = get_time();
  mpi.print("Prune tree    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.BottomUpTreeConstructor::grow(bodies2);
  toc = get_time();
  mpi.print("Grow tree     : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.sort(mpi.Ibody,bodies,bodies2,false);
  toc = get_time();
  mpi.print("Sort descend  : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.link();
  toc = get_time();
  mpi.print("Link cells    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  mpi.commBodies();
  toc = get_time();
  mpi.print("Comm bodies   : ",0);
  mpi.print(toc-tic);

  std::fill(mpi.Ibody.begin(),mpi.Ibody.end(),0);
#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  if( mpi.rank() == 0 ) {
    vtk.setDomain(mpi.getR0(),mpi.getX0());
    vtk.setGroupOfPoints(mpi.Ibody,bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=mpi.size(); ++i ) {
    mpi.shiftBodies(mpi.Ibody,bodies2);
    if( mpi.rank() == 0 )
      vtk.setGroupOfPoints(mpi.Ibody,bodies,Ncell);
  }
  toc = get_time();
  mpi.print("Shift bodies  : ",0);
  mpi.print(toc-tic);
  if( mpi.rank() == 0 )
    vtk.plot(Ncell);
#endif
}
