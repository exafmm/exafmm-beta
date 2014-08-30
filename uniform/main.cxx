#include "ewald.h"
#if Serial
#include "serialfmm.h"
#else
#include "parallelfmm.h"
#endif

int main() {
#if Serial
  SerialFMM FMM;
#else
  ParallelFMM FMM;
#endif
  const int numBodies = 100000;
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int gatherLevel = 1;
  const int numImages = 0;
  FMM.allocate(numBodies, maxLevel, numImages);
  logger::verbose = FMM.MPIRANK == 0;
  if( FMM.printNow ) {
    printf("N       : %d\n",FMM.numBodies);
    printf("Levels  : %d\n",FMM.maxLevel);
    printf("Images  : %d\n",FMM.numImages);
    printf("------------------\n");
  }

  logger::startTimer("Partition");
  FMM.partitioner(gatherLevel);
  logger::stopTimer("Partition");

  for( int it=0; it<1; it++ ) {
    int ix[3] = {0, 0, 0};
    FMM.R0 = 0.5;
    for_3d FMM.RGlob[d] = FMM.R0 * FMM.numPartition[FMM.maxGlobLevel][d];
    FMM.getGlobIndex(ix,FMM.MPIRANK,FMM.maxGlobLevel);
    for_3d FMM.X0[d] = 2 * FMM.R0 * (ix[d] + .5);
    srand48(FMM.MPIRANK);
    real average = 0;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][0] = 2 * FMM.R0 * (drand48() + ix[0]);
      FMM.Jbodies[i][1] = 2 * FMM.R0 * (drand48() + ix[1]);
      FMM.Jbodies[i][2] = 2 * FMM.R0 * (drand48() + ix[2]);
      FMM.Jbodies[i][3] = (drand48() - .5) / FMM.numBodies;
      average += FMM.Jbodies[i][3];
    }
    average /= FMM.numBodies;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][3] -= average;
    }
  
    logger::startTimer("Grow tree");
    FMM.sortBodies();
    FMM.buildTree();
    logger::stopTimer("Grow tree");
  
    logger::startTimer("Upward pass");
    FMM.upwardPass();
    logger::stopTimer("Upward pass");
  
#if Serial
#else
    logger::startTimer("Comm LET bodies");
    FMM.P2PSend();
    FMM.P2PRecv();
    logger::stopTimer("Comm LET bodies");

    logger::startTimer("Comm LET cells");
    for( int lev=FMM.maxLevel; lev>0; lev-- ) {
      MPI_Barrier(MPI_COMM_WORLD);
      FMM.M2LSend(lev);
      FMM.M2LRecv(lev);
    }
    FMM.rootGather();
    logger::stopTimer("Comm LET cells");
    FMM.globM2M();
    FMM.globM2L();
#endif
  
    FMM.periodicM2L();

#if Serial
#else
    logger::startTimer("Downward pass");
    FMM.globL2L();
    logger::stopTimer("Downward pass");
#endif
  
    FMM.downwardPass();
    if( FMM.printNow ) printf("------------------\n");

    Bodies bodies(FMM.numBodies);
    B_iter B = bodies.begin();
    for (int b=0; b<FMM.numBodies; b++, B++) {
      for_3d B->X[d] = FMM.Jbodies[b][d];
      B->SRC = FMM.Jbodies[b][3];
    }
    logger::startTimer("Total Direct");
#if Serial
    FMM.direct();
#else
    FMM.globDirect();
#endif
    if( FMM.printNow ) printf("------------------\n");
    logger::stopTimer("Total Direct");
  }
  FMM.deallocate();

  logger::startTimer("Attach root");
  logger::stopTimer("Attach root", 0);
  logger::startTimer("Comm partition");
  logger::stopTimer("Comm partition", 0);
  logger::startTimer("Get bounds");
  logger::stopTimer("Get bounds", 0);
  logger::startTimer("Link LET");
  logger::stopTimer("Link LET", 0);
  logger::startTimer("Link tree");
  logger::stopTimer("Link tree", 0);
  logger::startTimer("Root to body");
  logger::stopTimer("Root to body", 0);
  logger::startTimer("Set LET");
  logger::stopTimer("Set LET", 0);
  logger::startTimer("Set LET size");
  logger::stopTimer("Set LET size", 0);
  logger::writeTime(FMM.MPIRANK);

}
