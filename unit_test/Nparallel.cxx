#include "let.h"
#include "dataset.h"
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
  LocalEssentialTree T;
  T.setKernel(D.kernelName);
  T.initialize();
  bool printNow = MPIRANK == 0;

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+32)/8.0));
    if(printNow) std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.random(bodies,MPIRANK+1);
    T.startTimer("FMM          ");
    T.setGlobDomain(bodies);
    T.octsection(bodies);
    cells.clear();
#ifdef TOPDOWN
    T.topdown(bodies,cells);
#else
    T.bottomup(bodies,cells);
#endif
    T.commBodies(cells);
    jbodies = bodies;
    Cells jcells = cells;
    T.commCells(jbodies,jcells);

    T.downward(cells,jcells,1);
    T.stopTimer("FMM          ",printNow);
    T.eraseTimer("FMM          ");

    T.startTimer("Direct sum   ");
    Bodies bodies2 = bodies;
#if 1
    D.initTarget(bodies2);
    if( IMAGES != 0 ) {
      jbodies = T.periodicBodies(bodies2);
    } else {
      jbodies = bodies2;
    }
    bodies2.resize(numTarget);
    for( int i=0; i!=MPISIZE; ++i ) {
      T.shiftBodies(jbodies);
      T.evalP2P(bodies2,jbodies);
      if(T.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
    }
    D.writeTarget(bodies2);
#else
    D.readTarget(bodies2);
#endif
    T.stopTimer("Direct sum   ",printNow);
    T.eraseTimer("Direct sum   ");
    if(printNow) T.writeTime();
    if(printNow) T.resetTimer();

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
    bodies.resize(numTarget);
    D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
    MPI_Datatype MPI_TYPE = T.getType(diff1);
    MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
    if(printNow) D.printError(diff3,norm3,diff4,norm4);
  }
  T.finalize();
}
