#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells,jcells;
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

  jcells = cells;
  T.downward(cells,jcells,1);

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies = T.periodicBodies(bodies);
    T.stopTimer("Set periodic ",T.printNow);
  } else {
    jbodies = bodies;
  }

  T.startTimer("Direct sum   ");
  Evaluator E;
  T.buffer = bodies;
  for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
    B->pot = -B->Q / std::sqrt(EPS2);
    B->acc = 0;
  }
  T.evalP2P(T.buffer,jbodies);
  T.stopTimer("Direct sum   ",T.printNow);

  B_iter B  = bodies.begin();
  B_iter B2 = T.buffer.begin();
  real potDiff = 0, potNorm = 0, accDiff = 0, accNorm = 0;
  for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
    B->pot  -= B->Q / std::sqrt(EPS2);
#ifdef DEBUG
    std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;
#endif
    potDiff += (B->pot - B2->pot) * (B->pot - B2->pot);
    potNorm += B2->pot * B2->pot;
    accDiff += norm(B->acc - B2->acc);
    accNorm += norm(B2->acc);
  }
  std::cout << "Error (pot)   : " << std::sqrt(potDiff/potNorm) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(accDiff/accNorm) << std::endl;
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
