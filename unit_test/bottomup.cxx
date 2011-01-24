#include "dataset.h"
#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(10000000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies bodies2(numBodies);
  BottomUpTreeConstructor T(bodies);
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

  tic = get_time();
  T.setMorton();
  toc = get_time();
  std::cout << "Set Morton    : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2);
  toc = get_time();
  std::cout << "Sort Morton   : " << toc-tic << std::endl;

  tic = get_time();
  T.prune();
  toc = get_time();
  std::cout << "Prune tree    : " << toc-tic << std::endl;

  tic = get_time();
  T.grow(bodies2);
  toc = get_time();
  std::cout << "Grow tree     : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2,false);
  toc = get_time();
  std::cout << "Sort descend  : " << toc-tic << std::endl;

  tic = get_time();
  T.link();
  toc = get_time();
  std::cout << "Link cells    : " << toc-tic << std::endl;

#ifdef VTK
  T.sort(T.Ibody,bodies,bodies2);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
