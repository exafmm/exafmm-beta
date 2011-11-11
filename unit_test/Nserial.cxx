#include "dataset.h"
#include "serialfmm.h"
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
  Dataset dataset;
  dataset.kernelName = "Laplace";
  SerialFMM FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+32)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    dataset.sphere(bodies,1,1);
    FMM.startTimer("FMM          ");
    FMM.setDomain(bodies);
    cells.clear();
#ifdef TOPDOWN
    FMM.topdown(bodies,cells);
#else
    FMM.bottomup(bodies,cells);
#endif
    FMM.downward(cells,cells,2);
    FMM.stopTimer("FMM          ",true);
    FMM.eraseTimer("FMM          ");

    FMM.startTimer("Direct sum   ");
    FMM.buffer = bodies;
#if 1
    dataset.initTarget(FMM.buffer);
    if( IMAGES != 0 ) {
      jbodies = FMM.periodicBodies(FMM.buffer);
    } else {
      jbodies = FMM.buffer;
    }
    FMM.buffer.resize(numTarget);
    FMM.evalP2P(FMM.buffer,jbodies);
    dataset.writeTarget(FMM.buffer);
#else
    dataset.readTarget(FMM.buffer);
#endif
    FMM.stopTimer("Direct sum   ",true);
    FMM.eraseTimer("Direct sum   ");
    FMM.writeTime();
    FMM.resetTimer();

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    bodies.resize(numTarget);
    dataset.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2);
    dataset.printError(diff1,norm1,diff2,norm2);
  }
  FMM.finalize();
}
