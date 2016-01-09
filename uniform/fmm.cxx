#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree_tbb.h"
#include "dataset.h"
#include "ewald.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
#if EXAFMM_SERIAL
#include "../uniform/serialfmm.h"
#else
#include "../uniform/parallelfmm.h"
#endif
using namespace exafmm;

int main(int argc, char ** argv) {
  const int ksize = 11;
  const real_t cycle = 20 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = 20;

  Args args(argc, argv);
  args.ncrit = 32;
  args.images = 1;
  BaseMPI baseMPI;
  BoundBox boundBox(args.nspawn);
  BuildTree buildTree(args.ncrit, args.nspawn);
  Dataset data;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);
  Traversal traversal(args.nspawn, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
#if EXAFMM_SERIAL
  SerialFMM FMM;
#else
  ParallelFMM FMM;
#endif
  TreeMPI treeMPI(FMM.MPIRANK, FMM.MPISIZE, args.images);

  args.numBodies /= FMM.MPISIZE;
  const int numBodies = args.numBodies;
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int gatherLevel = 1;
  const int numImages = args.images;
  if (numImages > 0 && int(log2(FMM.MPISIZE)) % 3 != 0) {
    if (FMM.MPIRANK==0) printf("Warning: MPISIZE must be a power of 8 for periodic domain to be square\n");
  }

  FMM.allocate(numBodies, maxLevel, numImages);
  args.verbose &= FMM.MPIRANK == 0;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, EXAFMM_PP);

  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startTimer("Partition");
  FMM.partitioner(gatherLevel);
  logger::stopTimer("Partition");

  for( int it=0; it<1; it++ ) {
    int ix[3] = {0, 0, 0};
    FMM.R0 = 0.5 * cycle / FMM.numPartition[FMM.maxGlobLevel][0];
    for_3d FMM.RGlob[d] = FMM.R0 * FMM.numPartition[FMM.maxGlobLevel][d];
    FMM.getGlobIndex(ix,FMM.MPIRANK,FMM.maxGlobLevel);
    for_3d FMM.X0[d] = 2 * FMM.R0 * (ix[d] + .5);
    srand48(FMM.MPIRANK);
    real_t average = 0;
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

#if EXAFMM_SERIAL
#else
    logger::startTimer("Comm LET bodies");
    FMM.P2PSend();
    FMM.P2PRecv();
    logger::stopTimer("Comm LET bodies");
#endif
  
    logger::startTimer("Upward pass");
    FMM.upwardPass();
    logger::stopTimer("Upward pass");
  
#if EXAFMM_SERIAL
#else
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

#if EXAFMM_SERIAL
#else
    logger::startTimer("Downward pass");
    FMM.globL2L();
    logger::stopTimer("Downward pass", 0);
#endif
  
    FMM.downwardPass();
    /*
    logger::stopTimer("Total FMM", 0);

    Bodies bodies(FMM.numBodies);
    B_iter B = bodies.begin();
    for (int b=0; b<FMM.numBodies; b++, B++) {
      for_3d B->X[d] = FMM.Jbodies[b][d];
      B->SRC = FMM.Jbodies[b][3];
      for_4d B->TRG[d] = FMM.Ibodies[b][d];
    }
    Bodies jbodies = bodies;
    vec3 localDipole = upDownPass.getDipole(bodies, FMM.RGlob[0]);
    vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
    int numBodies = baseMPI.allreduceInt(bodies.size());
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
#ifndef EXAFMM_IJHPCA
#if 1
    logger::startTimer("Total Ewald");
    Bounds bounds = boundBox.getBounds(bodies);
    Bodies buffer = bodies;
    Cells cells = buildTree.buildTree(bodies, buffer, bounds);
    Bodies bodies2 = bodies;
    data.initTarget(bodies);
    for (int i=0; i<FMM.MPISIZE; i++) {
      if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << FMM.MPISIZE << std::endl;
      if (FMM.MPISIZE > 1) treeMPI.shiftBodies(jbodies);
      bounds = boundBox.getBounds(jbodies);
      buffer = jbodies;
      Cells jcells = buildTree.buildTree(jbodies, buffer, bounds);
      ewald.wavePart(bodies, jbodies);
      ewald.realPart(cells, jcells);
    }
    ewald.selfTerm(bodies);
#else
    logger::startTimer("Total Direct");
    const int numTargets = 100;
    data.sampleBodies(bodies, numTargets);
    Bodies bodies2 = bodies;
    data.initTarget(bodies);
    for (int i=0; i<FMM.MPISIZE; i++) {
      if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << FMM.MPISIZE << std::endl;
      if (FMM.MPISIZE > 1) treeMPI.shiftBodies(jbodies);
      traversal.direct(bodies, jbodies, cycle);
    }
    traversal.normalize(bodies);
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    logger::stopTimer("Total Direct");
    logger::resetTimer("Total Direct");
#endif
    Verify verify;
    double potSum = verify.getSumScalar(bodies);
    double potSum2 = verify.getSumScalar(bodies2);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    logger::printTitle("FMM vs. direct");
#if EXAFMM_SERIAL
    double potDif = (potSum - potSum2) * (potSum - potSum2);
    double potNrm = potSum * potSum;
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
#else
    double potSumGlob, potSumGlob2, accDifGlob, accNrmGlob;
    MPI_Reduce(&potSum,  &potSumGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potSum2, &potSumGlob2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accDif,  &accDifGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accNrm,  &accNrmGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
    double potNrmGlob = potSumGlob * potSumGlob;
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
#endif
#endif
    */
  }
  FMM.deallocate();

#ifndef EXAFMM_IJHPCA
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
#endif
}
