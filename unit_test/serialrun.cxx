#include "dataset.h"
#include "construct.h"
#include "kernel.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(10000);
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
  buffer = bodies;
  K.P2P(buffer.begin(),buffer.end());
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = buffer.begin();
  real err(0),rel(0);
  for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
    B->pot -= B->scal / std::sqrt(EPS2);                        //  Initialize body values
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
#ifdef VTK
  T.sort(bodies,buffer);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
