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
  TreeConstructor T;
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere(bodies);
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain(bodies);
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  T.BottomUp::setIndex(bodies);
  toc = get_time();
  std::cout << "Set index     : " << toc-tic << std::endl;

  tic = get_time();
  T.buffer.resize(bodies.size());
  T.sort(bodies,T.buffer);
  toc = get_time();
  std::cout << "Sort index    : " << toc-tic << std::endl;

  bigint oldIndex(bodies[0].I);
  int b = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {
    assert( oldIndex <= B->I );
    oldIndex = B->I;
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
