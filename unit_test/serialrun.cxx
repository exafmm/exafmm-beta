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
  Cells cells;
  Dataset D;
  TreeConstructor T;
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere(bodies,1,1);
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain(bodies);
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  tic = get_time();
  T.setCI0(cells.begin());
  T.setCJ0(cells.begin());
  T.upward(cells);
  toc = get_time();
  std::cout << "Upward        : " << toc-tic << std::endl;

  tic = get_time();
  T.downward(cells,1);
  toc = get_time();
  std::cout << "Downward      : " << toc-tic << std::endl;

  tic = get_time();
  Kernel K;
  T.buffer = bodies;
  K.P2P(T.buffer.begin(),T.buffer.end());
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = T.buffer.begin();
  real err(0),rel(0);
  for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
    B->pot -= B->scal / std::sqrt(EPS2);
#ifdef DEBUG
    std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;
#endif
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
#ifdef VTK
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
