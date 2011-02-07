#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(10000000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies buffer(numBodies);
  TreeConstructor T(bodies);
  Dataset D(bodies);
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere();
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain();
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

#ifdef TOPDOWN
  T.topdown(buffer);
#else
  T.bottomup(buffer);
#endif

#ifdef VTK
  T.sort(T.Ibody,bodies,buffer);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
