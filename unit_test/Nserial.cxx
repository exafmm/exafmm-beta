#include "dataset.h"
#include "construct.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 10000;
  int numTarget = 100;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies, jbodies;
  Cells cells;
  Dataset D;
  D.kernelName = "Laplace";
  TreeConstructor T;
  T.setKernel(D.kernelName);
  T.initialize();

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+32)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.sphere(bodies,1,1);
    T.startTimer("FMM          ");
    T.setDomain(bodies);
    cells.clear();
#ifdef TOPDOWN
    T.topdown(bodies,cells);
#else
    T.bottomup(bodies,cells);
#endif
    T.downward(cells,cells,2);
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
    T.buffer.resize(numTarget);
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
    bodies.resize(numTarget);
    D.evalError(bodies,T.buffer,diff1,norm1,diff2,norm2);
    D.printError(diff1,norm1,diff2,norm2);
  }
  T.finalize();
}
