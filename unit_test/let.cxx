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
  Bodies buffer(numBodies);
  Dataset D(bodies);
  LocalEssentialTree mpi(bodies);
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
  mpi.setGlobDomain(mpi.setR0(),mpi.setX0());
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  mpi.multisection(mpi.Ibody,buffer);
  toc = get_time();
  if(print) std::cout << "Multisection  : " << toc-tic << std::endl;

#ifdef TOPDOWN
  mpi.topdown(buffer,print);
#else
  mpi.bottomup(buffer,print);
#endif

  tic = get_time();
  mpi.commBodies();
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  mpi.commCells();
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

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
    mpi.shiftBodies(mpi.Ibody,buffer);
    if( mpi.rank() == 0 )
      vtk.setGroupOfPoints(mpi.Ibody,bodies,Ncell);
  }
  toc = get_time();
  if(print) std::cout << "Shift bodies  : " << toc-tic << std::endl;
  if( mpi.rank() == 0 )
    vtk.plot(Ncell);
#endif
}
