#include "dataset.h"
#include "bottomup.h"
#include "kernel.h"
#ifdef VTK
#include "vtk.h"
#endif

int main()
{
  double tic,toc;
  int const Nbody=10000;
  tic = get_time();
  Bodies bodies(Nbody);
  Bodies bodies2(Nbody);
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

  tic = get_time();
  T.P2M();
  toc = get_time();
  std::cout << "P2M           : " << toc-tic << std::endl;

  tic = get_time();
  T.M2M();
  toc = get_time();
  std::cout << "M2M           : " << toc-tic << std::endl;

  tic = get_time();
  T.evaluate(1);
  toc = get_time();
  std::cout << "Evaluate      : " << toc-tic << std::endl;

  tic = get_time();
  Kernel K;
  bodies2 = bodies;
  K.P2P(bodies2.begin(),bodies2.end());
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = bodies2.begin();
  real err(0), rel(0);
  for( int i=0; i!=Nbody; ++i,++B,++B2 ) {
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
#ifdef VTK
  T.sort(T.Ibody,bodies,bodies2);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
