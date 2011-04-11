#include "dataset.h"
#include "construct.h"
#include "kernel.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 1000;
  std::string kernelName = "Laplace";
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells;
  Dataset D;
  TreeConstructor T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+24)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.random(bodies,1,1);
    T.startTimer("FMM          ");
    T.setDomain(bodies);
    cells.clear();
#ifdef TOPDOWN
    T.topdown(bodies,cells);
#else
    T.bottomup(bodies,cells);
#endif
    T.downward(cells,cells,1);
    T.stopTimer("FMM          ",true);
    T.eraseTimer("FMM          ");

    T.startTimer("Direct sum   ");
    T.buffer = bodies;
#if 1
    D.initTarget(T.buffer);
    if( IMAGES != 0 ) {
      jbodies = T.periodicBodies(T.buffer);
    } else {
      jbodies = T.buffer;
    }
    T.evalP2P(T.buffer,jbodies);
    D.writeTarget(T.buffer);
#else
    D.readTarget(T.buffer);
#endif
    T.stopTimer("Direct sum   ",true);
    T.eraseTimer("Direct sum   ");
    T.writeTime();
    T.resetTimer();

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    D.evalError(bodies,T.buffer,diff1,norm1,diff2,norm2);
    D.printError(diff1,norm1,diff2,norm2);
  }
#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  T.finalize();
}
