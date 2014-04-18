#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "ewald.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
#if VTK
#include "vtk.h"
#endif
#if MASS
#error Turn off MASS for this test
#endif
#if EXPANSION < 10
#error Use P >=10 for this test
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  BaseMPI baseMPI;
  Bodies bodies, bodies2, bodies3, jbodies, gbodies;
  BoundBox boundBox(args.nspawn);
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit, args.nspawn);
  BuildTree globalTree(1, args.nspawn);
  Cells cells, jcells;
  Dataset data;
  Partition partition(baseMPI.mpirank, baseMPI.mpisize);
  Traversal traversal(args.nspawn, args.images);
  TreeMPI treeMPI(baseMPI.mpirank, baseMPI.mpisize, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;
  num_threads(args.threads);

  const int ksize = 11;
  const real_t cycle = 2 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = cycle * alpha / 3;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);

  args.verbose &= baseMPI.mpirank == 0;
  logger::verbose = args.verbose;
  logger::printTitle("Ewald Parameters");
  args.print(logger::stringLength, P);
  ewald.print(logger::stringLength);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  //data.writeSources(bodies, baseMPI.mpirank);
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    localBounds = boundBox.getBounds(bodies);
    globalBounds = baseMPI.allreduceBounds(localBounds);
    localBounds = partition.octsection(bodies, globalBounds);
    bodies = treeMPI.commBodies(bodies);

    cells = localTree.buildTree(bodies, localBounds);
    upDownPass.upwardPass(cells);
    treeMPI.allgatherBounds(localBounds);
    treeMPI.setLET(cells, cycle);
    treeMPI.commBodies();
    treeMPI.commCells();

    traversal.initWeight(cells);
    traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
    if (args.graft) {
      treeMPI.linkLET();
      gbodies = treeMPI.root2body();
      jcells = globalTree.buildTree(gbodies, globalBounds);
      treeMPI.attachRoot(jcells);
      traversal.dualTreeTraversal(cells, jcells, cycle, false);
    } else {
      for (int irank=0; irank<baseMPI.mpisize; irank++) {
	treeMPI.getLET(jcells, (baseMPI.mpirank+irank)%baseMPI.mpisize);
	traversal.dualTreeTraversal(cells, jcells, cycle, false);
      }
    }
    upDownPass.downwardPass(cells);
    vec3 localDipole = upDownPass.getDipole(bodies,0);
    vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
    int numBodies = baseMPI.allreduceInt(bodies.size());
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
#if 1
    bodies2 = bodies;
    data.initTarget(bodies);
    logger::printTitle("Ewald Profiling");
    logger::startTimer("Total Ewald");
#if 1
    jbodies = bodies;
    for (int i=0; i<baseMPI.mpisize; i++) {
      if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << baseMPI.mpisize << std::endl;
      treeMPI.shiftBodies(jbodies);
      localBounds = boundBox.getBounds(jbodies);
      jcells = localTree.buildTree(jbodies, localBounds);
      ewald.wavePart(bodies, jbodies);
      ewald.realPart(cells, jcells);
    }
#else
    jbodies = treeMPI.allgatherBodies(bodies);
    jcells = localTree.buildTree(jbodies, globalBounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
#endif

    ewald.selfTerm(bodies);
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    logger::stopTimer("Total Ewald");
#else
    jbodies = bodies;
    const int numTargets = 100;
    bodies3 = bodies;
    data.sampleBodies(bodies, numTargets);
    bodies2 = bodies;
    data.initTarget(bodies);
    logger::startTimer("Total Direct");
    for (int i=0; i<baseMPI.mpisize; i++) {
      treeMPI.shiftBodies(jbodies);
      traversal.direct(bodies, jbodies, cycle);
      if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
    }
    traversal.normalize(bodies);
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    logger::stopTimer("Total Direct");
    bodies = bodies3;
#endif
    double potSum = verify.getSumScalar(bodies);
    double potSum2 = verify.getSumScalar(bodies2);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    double potSumGlob, potSumGlob2, accDifGlob, accNrmGlob;
    MPI_Reduce(&potSum,  &potSumGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potSum2, &potSumGlob2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accDif,  &accDifGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accNrm,  &accNrmGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    logger::printTitle("FMM vs. Ewald");
    double potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
    double potNrmGlob = potSumGlob * potSumGlob;
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
    localTree.printTreeData(cells);
    traversal.printTraversalData();
    logger::printPAPI();
    data.initTarget(bodies);
  }
#if VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->IBODY = 0;
  for (int irank=0; irank<baseMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(baseMPI.mpirank+irank)%baseMPI.mpisize);
    for (C_iter C=jcells.begin(); C!=jcells.end(); C++) {
      Body body;
      body.IBODY = 1;
      body.X     = C->X;
      body.SRC   = 0;
      jbodies.push_back(body);
    }
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(jbodies);
  for (int i=1; i<baseMPI.mpisize; i++) {
    treeMPI.shiftBodies(jbodies);
    vtk.setGroupOfPoints(jbodies);
  }
  if (baseMPI.mpirank == 0) {
    vtk.plot();
  }
#endif
  return 0;
}
