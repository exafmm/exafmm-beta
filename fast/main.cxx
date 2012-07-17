#include "dataset.h"
#include "parallelfmm.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies;
  Cells cells;
  Dataset DATA;
  ParallelFMM FMM;
#if HYBRID
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
  numBodies = int(pow(10,(it+24)/8.0));
#else
  {
  FMM.printNow = MPIRANK == 0;
#if BUILD
  numBodies = 10000000;
#else
  numBodies = 1000000;
#endif
#endif // MANY
  if(FMM.printNow) std::cout << "N                    : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.cube(bodies);
  FMM.startTimer("FMM");

  FMM.octsection(bodies);
#if BOTTOMUP
  FMM.bottomup(bodies,cells);
#else
  FMM.topdown(bodies,cells);
#endif
#if BUILD
#else
  FMM.startPAPI();
#if IneJ
  Bodies jbodies = bodies;
  for( int irank=0; irank!=MPISIZE; irank++ ) {
    FMM.shiftBodies(jbodies);
    Cells jcells;
    FMM.bottomup(jbodies,jcells);
    FMM.evaluate(cells,jcells);
  }
#else
  FMM.evaluate(cells);
#endif
  FMM.stopPAPI();
  FMM.stopTimer("FMM",FMM.printNow);
  FMM.eraseTimer("FMM");
  FMM.writeTime();
  FMM.resetTimer();

  if (bodies.size() > 100) bodies.resize(100);
  Bodies bodies2 = bodies;
  DATA.initTarget(bodies2);
  FMM.startTimer("Direct sum");
  for( int i=0; i!=MPISIZE; ++i ) {
    FMM.shiftBodies(jbodies);
    FMM.direct(bodies2,jbodies);
    if(FMM.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  FMM.stopTimer("Direct sum",FMM.printNow);
  FMM.eraseTimer("Direct sum");
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = FMM.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(FMM.printNow) DATA.printError(diff3,norm3,diff4,norm4);
#endif // BUILD
  }
}
