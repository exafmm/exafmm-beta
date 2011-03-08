#include "dataset.h"
#include "construct.h"
#include "kernel.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 1000;
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  Evaluator E;
  TreeConstructor T;

  for( int it=0; it!=17; ++it ) {
    numBodies = int(pow(10,(it+24)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.sphere(bodies,1,1);
    T.setDomain(bodies);
#ifdef TOPDOWN
    T.topdown(bodies,cells,false);
#else
    T.bottomup(bodies,cells,false);
#endif

    T.startTimer("Downward     ");
    T.downward(cells,cells,1);
    T.stopTimer("Downward     ");

    T.startTimer("Direct sum   ");
    T.buffer = bodies;
    for( B_iter B=T.buffer.begin(); B!=T.buffer.end(); ++B ) {
      B->pot = -B->scal / std::sqrt(EPS2);
    }
    E.evalP2P(T.buffer,T.buffer);
    T.stopTimer("Direct sum   ");
    T.printAllTime();
    T.eraseTimer("Direct sum   ");
    T.writeTime();
    T.resetTimer();

    B_iter B  = bodies.begin();
    B_iter B2 = T.buffer.begin();
    real err = 0, rel = 0;
    for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
      B->pot -= B->scal / std::sqrt(EPS2);
#ifdef DEBUG
      std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;
#endif
      err += (B->pot - B2->pot) * (B->pot - B2->pot);
      rel += B2->pot * B2->pot;
    }
    std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
    cells.clear();
  }
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
