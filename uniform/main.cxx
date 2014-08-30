#include "base_mpi.h"
#include "args.h"
#include "dataset.h"
#include "ewald.h"
#include "traversal.h"
#include "tree_mpi.h"
#if Serial
#include "serialfmm.h"
#else
#include "parallelfmm.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Traversal traversal(args.nspawn, args.images);
#if Serial
  SerialFMM FMM;
#else
  ParallelFMM FMM;
#endif
  TreeMPI treeMPI(FMM.MPIRANK, FMM.MPISIZE, args.images);

  const int numBodies = args.numBodies;
  const int ncrit = args.ncrit;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int gatherLevel = 1;
  const int numImages = args.images;
  const real cycle = 1;
  FMM.allocate(numBodies, maxLevel, numImages);
  args.verbose &= FMM.MPIRANK == 0;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, PP);

  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startTimer("Partition");
  FMM.partitioner(gatherLevel);
  logger::stopTimer("Partition");

  for( int it=0; it<1; it++ ) {
    int ix[3] = {0, 0, 0};
    FMM.R0 = 0.5 * cycle;
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
    logger::stopTimer("Comm LET cells", 0);
    FMM.globM2M();
    FMM.globM2L();
#endif
  
    FMM.periodicM2L();

#if Serial
#else
    logger::startTimer("Downward pass");
    FMM.globL2L();
    logger::stopTimer("Downward pass", 0);
#endif
  
    FMM.downwardPass();
    logger::stopTimer("Total FMM", 0);

    Bodies bodies(FMM.numBodies);
    B_iter B = bodies.begin();
    for (int b=0; b<FMM.numBodies; b++, B++) {
      for_3d B->X[d] = FMM.Jbodies[b][d];
      B->SRC = FMM.Jbodies[b][3];
      for_4d B->TRG[d] = FMM.Ibodies[b][d];
    }
    Bodies jbodies = bodies;
    const int numTargets = 100;
    data.sampleBodies(bodies, numTargets);
    Bodies bodies2 = bodies;
    data.initTarget(bodies);
    
    logger::startTimer("Total Direct");
    logger::printTitle("MPI direct sum");
    for (int i=0; i<FMM.MPISIZE; i++) {
      if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << FMM.MPISIZE << std::endl;
      treeMPI.shiftBodies(jbodies);
      traversal.direct(bodies, jbodies, cycle);
    }
    logger::printTitle("FMM vs. direct");
#if Serial
    FMM.direct();
#else
    FMM.globDirect();
#endif
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
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
