#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  TreeConstructor T;
  T.printNow = true;

  T.startTimer("Set bodies   ");
  D.sphere(bodies,1,1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  T.startTimer("Downward     ");
  T.downward(cells,cells,1);
  T.stopTimer("Downward     ",T.printNow);

  T.startTimer("Direct sum   ");
  Evaluator E;
  T.buffer = bodies;
  for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
    B->pot = -B->scal / std::sqrt(EPS2);
  }
  T.evalP2P(T.buffer,T.buffer);
  T.stopTimer("Direct sum   ",T.printNow);

  B_iter B  = bodies.begin();
  B_iter B2 = T.buffer.begin();
  real err = 0, rel = 0;
  for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
    B->pot  -= B->scal / std::sqrt(EPS2);
#ifdef DEBUG
    std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;
#endif
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
