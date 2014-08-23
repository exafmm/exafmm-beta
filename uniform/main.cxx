#if Serial
#include "serialfmm.h"
#else
#include "parallelfmm.h"
#endif

int main() {
  double tic, toc;
#if Serial
  SerialFMM FMM;
#else
  ParallelFMM FMM;
#endif
  char fname[256];
  sprintf(fname,"oldtime%5.5d.dat",FMM.MPIRANK);
  std::ofstream fid(fname);
  const int numBodies = 100000;
  const int maxLevel = 4;
  const int gatherLevel = 1;
  const int numImages = 2;
  FMM.allocate(numBodies, maxLevel, numImages);
  //logger::verbose = true;
  //FMM.printNow = false;
  if( FMM.printNow ) {
    printf("N       : %d\n",FMM.numBodies);
    printf("Levels  : %d\n",FMM.maxLevel);
    printf("Images  : %d\n",FMM.numImages);
    printf("------------------\n");
  }

  logger::startTimer("Partition");
  tic = FMM.getTime();
  FMM.partitioner(gatherLevel);
  toc = FMM.getTime();
  if( FMM.printNow ) printf("Part    : %lf\n",toc-tic);
  logger::stopTimer("Partition");

  for( int it=0; it<1; it++ ) {
    int ix[3] = {0, 0, 0};
    FMM.R0 = .5;
    for_3d FMM.RGlob[d] = FMM.R0 * FMM.numPartition[FMM.maxGlobLevel][d];
    FMM.getGlobIndex(ix,FMM.MPIRANK,FMM.maxGlobLevel);
    for_3d FMM.X0[d] = 2 * FMM.R0 * (ix[d] + .5);
    srand48(FMM.MPIRANK);
    real average = 0;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][0] = drand48() + 2 * FMM.R0 * ix[0];
      FMM.Jbodies[i][1] = drand48() + 2 * FMM.R0 * ix[1];
      FMM.Jbodies[i][2] = drand48() + 2 * FMM.R0 * ix[2];
      FMM.Jbodies[i][3] = (drand48() - .5) / FMM.numBodies;
      average += FMM.Jbodies[i][3];
    }
    average /= FMM.numBodies;
    for( int i=0; i<FMM.numBodies; i++ ) {
      FMM.Jbodies[i][3] -= average;
    }
  
    logger::startTimer("Grow tree");
    tic = FMM.getTime();
    FMM.sortBodies();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "Sort " << toc-tic << std::endl;
    if( FMM.printNow ) printf("Sort    : %lf\n",toc-tic);
  
    tic = FMM.getTime();
    FMM.buildTree();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "Tree " << toc-tic << std::endl;
    if( FMM.printNow ) printf("Tree    : %lf\n",toc-tic);
    logger::stopTimer("Grow tree");
  
    logger::startTimer("Upward pass");
    tic = FMM.getTime();
    FMM.upwardPass();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "Upward " << toc-tic << std::endl;
    logger::stopTimer("Upward pass");
  
#if Serial
#else
    logger::startTimer("Comm LET bodies");
    tic = FMM.getTime();
    FMM.P2PSend();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "P2P Send " << toc-tic << std::endl;
    if( FMM.printNow ) printf("P2P Send: %lf @ lev: %d\n",toc-tic,FMM.maxLevel);

    tic = FMM.getTime();
    FMM.P2PRecv();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "P2P Recv " << toc-tic << std::endl;
    if( FMM.printNow ) printf("P2P Recv: %lf @ lev: %d\n",toc-tic,FMM.maxLevel);
    logger::stopTimer("Comm LET bodies");

    logger::startTimer("Comm LET cells");
    for( int lev=FMM.maxLevel; lev>0; lev-- ) {
      MPI_Barrier(MPI_COMM_WORLD);
      tic = FMM.getTime();
      FMM.M2LSend(lev);
      //toc = FMM.getTime();
      //if( FMM.printNow ) fid << std::setw(20) << std::left << "M2L Send " << toc-tic << std::endl;
      //if( FMM.printNow ) printf("M2L Send: %lf @ lev: %d\n",toc-tic,lev);

      //tic = FMM.getTime();
      FMM.M2LRecv(lev);
      toc = FMM.getTime();
      fid << toc-tic << std::endl;
      //if( FMM.printNow ) fid << std::setw(20) << std::left << "M2L Recv " << toc-tic << std::endl;
      if( FMM.printNow ) printf("M2L Recv: %lf @ lev: %d\n",toc-tic,lev);
    }

    tic = FMM.getTime();
    FMM.rootGather();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "Gather " << toc-tic << std::endl;
    if( FMM.printNow ) printf("Gather  : %lf\n",toc-tic);
    logger::stopTimer("Comm LET cells");
  
    FMM.globM2M();
  
    FMM.globM2L(fid);
#endif
  
    tic = FMM.getTime();
    FMM.periodicM2L();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "M2L Peri " << toc-tic << std::endl;
    if( FMM.printNow ) printf("M2L Peri: %lf\n",toc-tic);
  
#if Serial
#else
    logger::startTimer("Downward pass");
    tic = FMM.getTime();
    FMM.globL2L();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "L2L Glob " << toc-tic << std::endl;
    if( FMM.printNow ) printf("L2L Glob: %lf\n",toc-tic);
    logger::stopTimer("Downward pass");
#endif
  
    tic = FMM.getTime();
    FMM.downwardPass();
    toc = FMM.getTime();
    //if( FMM.printNow ) fid << std::setw(20) << std::left << "Downward " << toc-tic << std::endl;
    if( FMM.printNow ) {
      printf("Downward: %lf\n",toc-tic);
      printf("------------------\n");
    }

#if 1
    tic = FMM.getTime();
#if Serial
    FMM.direct();
#else
    FMM.globDirect();
#endif
    toc = FMM.getTime();
    if( FMM.printNow ) {
      printf("------------------\n");
      printf("Direct  : %lf\n",toc-tic);
    }
#endif
  }
  FMM.deallocate();

  logger::startTimer("Attach root");
  logger::stopTimer("Attach root");
  logger::startTimer("Comm partition");
  logger::stopTimer("Comm partition");
  logger::startTimer("Get bounds");
  logger::stopTimer("Get bounds");
  logger::startTimer("Link LET");
  logger::stopTimer("Link LET");
  logger::startTimer("Link tree");
  logger::stopTimer("Link tree");
  logger::startTimer("Root to body");
  logger::stopTimer("Root to body");
  logger::startTimer("Set LET");
  logger::stopTimer("Set LET");
  logger::startTimer("Set LET size");
  logger::stopTimer("Set LET size");
  logger::writeTime(FMM.MPIRANK);

}
