#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  BaseMPI baseMPI;
  Bodies bodies, bodies2, bodies3, jbodies, gbodies;
  BoundBox boundBox(args.nspawn);
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit, args.nspawn);
  BuildTree globalTree(1, args.nspawn);
  Cells cells, jcells, gcells;
  Dataset data;
  Partition partition(baseMPI.mpirank, baseMPI.mpisize);
  Traversal traversal(args.nspawn, args.images);
  TreeMPI treeMPI(baseMPI.mpirank, baseMPI.mpisize, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;
  num_threads(args.threads);

  const real_t cycle = 2 * M_PI;
  //args.numBodies /= baseMPI.mpisize;
  args.verbose &= baseMPI.mpirank == 0;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
#if IneJ
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] += M_PI;
    B->X[0] *= 0.5;
  }
  jbodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank+baseMPI.mpisize, baseMPI.mpisize);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    B->X[0] -= M_PI;
    B->X[0] *= 0.5;
  }
#endif
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    localBounds = boundBox.getBounds(bodies);
#if IneJ
    localBounds = boundBox.getBounds(jbodies, localBounds);
#endif
    globalBounds = baseMPI.allreduceBounds(localBounds);
    partition.bisection(bodies, globalBounds);
    bodies = treeMPI.commBodies(bodies);
#if IneJ
    partition.bisection(jbodies, globalBounds);
    jbodies = treeMPI.commBodies(jbodies);
#endif
    localBounds = boundBox.getBounds(bodies);
    cells = localTree.buildTree(bodies, localBounds);
    localBounds = boundBox.getBounds(cells, localBounds);
    upDownPass.upwardPass(cells);
#if IneJ
    localBounds = boundBox.getBounds(jbodies);
    jcells = localTree.buildTree(jbodies, localBounds);
    localBounds = boundBox.getBounds(jcells, localBounds);
    upDownPass.upwardPass(jcells);
#endif

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree
    treeMPI.allgatherBounds(localBounds);
#if IneJ
    treeMPI.setLET(jcells, cycle);
#else
    treeMPI.setLET(cells, cycle);
#endif
#pragma omp parallel sections
    {
#pragma omp section
      {
	treeMPI.commBodies();
	treeMPI.commCells();
      }
#pragma omp section
      {
	traversal.initWeight(cells);
#if IneJ
	traversal.dualTreeTraversal(cells, jcells, cycle, false);
#else
	traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
	jbodies = bodies;
#endif
      }
    }
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
#else
    jbodies = bodies;
    for (int irank=0; irank<baseMPI.mpisize; irank++) {
      treeMPI.shiftBodies(jbodies);
      jcells.clear();
      localBounds = boundBox.getBounds(jbodies);
      jcells = localTree.buildTree(jbodies, localBounds);
      upDownPass.upwardPass(jcells);
      traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
    }
#endif
    upDownPass.downwardPass(cells);

    logger::stopPAPI();
    logger::stopTimer("Total FMM", 0);
#if DIRECT
    logger::printTitle("MPI direct sum");
    const int numTargets = 100;
    bodies3 = bodies;
    data.sampleBodies(bodies, numTargets);
    bodies2 = bodies;
    data.initTarget(bodies);
    logger::startTimer("Total Direct");
    for (int i=0; i<baseMPI.mpisize; i++) {
      if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
      treeMPI.shiftBodies(jbodies);
      traversal.direct(bodies, jbodies, cycle);
    }
    traversal.normalize(bodies);
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    logger::stopTimer("Total Direct");
    logger::resetTimer("Total Direct");
    double potDif = verify.getDifScalar(bodies, bodies2);
    double potNrm = verify.getNrmScalar(bodies);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    double potDifGlob, potNrmGlob, accDifGlob, accNrmGlob;
    MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accDif, &accDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&accNrm, &accNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    logger::printTitle("FMM vs. direct");
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
    localTree.printTreeData(cells);
    traversal.printTraversalData();
    logger::printPAPI();
    bodies = bodies3;
    data.initTarget(bodies);
#endif
    logger::resetTimer("Total FMM");
#if WRITE_TIME
    logger::writeTime(baseMPI.mpirank);
#endif
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
