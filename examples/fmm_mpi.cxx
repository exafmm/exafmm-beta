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
#include "kernel_select.h"
using namespace exafmm;
vec3 Kernel::Xperiodic = 0;
double Kernel::eps2 = 0.0;
#if EXAFMM_HELMHOLTZ
complex_t Kernel::wavek = complex_t(10.,1.) / real_t(2 * M_PI);
#endif

int main(int argc, char ** argv) {
  const vec3 cycle = 2 * M_PI;
  Args args(argc, argv);
  BaseMPI baseMPI;
  Bodies bodies, bodies2, jbodies, gbodies, buffer;
  BoundBox boundBox(args.nspawn);
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit, args.nspawn);
  BuildTree globalTree(1, args.nspawn);
  Cells cells, jcells, gcells;
  Dataset data;
  Partition partition(baseMPI.mpirank, baseMPI.mpisize);
  TreeMPI<kernel> treeMPI(baseMPI.mpirank, baseMPI.mpisize, args.images);
  Traversal<kernel> traversal(args.nspawn, args.images);  
  UpDownPass<kernel> upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;
  num_threads(args.threads);

  kernel::setup();
  args.numBodies /= baseMPI.mpisize;
  args.verbose &= baseMPI.mpirank == 0;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  buffer.reserve(bodies.size());
  if (args.IneJ) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    jbodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank+baseMPI.mpisize, baseMPI.mpisize);
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
  }
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    localBounds = boundBox.getBounds(bodies);
    if (args.IneJ) {
      localBounds = boundBox.getBounds(jbodies, localBounds);
    }
    globalBounds = baseMPI.allreduceBounds(localBounds);
    partition.partitionHilbert(bodies, globalBounds);    
    bodies = treeMPI.commBodies(bodies);
    if (args.IneJ) {
      partition.partitionHilbert(bodies, globalBounds);      
      jbodies = treeMPI.commBodies(jbodies);
    }
    localBounds = boundBox.getBounds(bodies);
    cells = localTree.buildTree(bodies, buffer, localBounds);
    localBounds = boundBox.getBounds(cells, localBounds);
    upDownPass.upwardPass(cells);
    if (args.IneJ) {
      localBounds = boundBox.getBounds(jbodies);
      jcells = localTree.buildTree(jbodies, buffer, localBounds);
      localBounds = boundBox.getBounds(jcells, localBounds);
      upDownPass.upwardPass(jcells);
    }
#if 1 
  treeMPI.allgatherBounds(localBounds);
  if (args.IneJ) {  
    treeMPI.setSendLET(jcells, bodies);
  } else {
    treeMPI.setSendLET(cells, bodies);
  }
  traversal.initListCount(cells);
  traversal.initWeight(cells);
  if (args.IneJ) {
    traversal.traverse(cells, jcells, cycle, args.dual, false);
  } else {
    traversal.traverse(cells, cells, cycle, args.dual, args.mutual);
    jbodies = bodies;
  }
  treeMPI.dualTreeTraversalRemote(cells,bodies,baseMPI.mpirank,baseMPI.mpisize,args.nspawn, args.granularity);
  #pragma omp parallel sections
    {
#pragma omp section
      {
    treeMPI.flushAllRequests();
      }
#pragma omp section
      {
    upDownPass.downwardPass(cells);
    logger::stopPAPI();
    logger::stopTimer("Total FMM", 0);
      }
    }        
#else
#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree
    treeMPI.allgatherBounds(localBounds);
    if (args.IneJ) {
      treeMPI.setLET(jcells, cycle);
    } else {
      treeMPI.setLET(cells, cycle);
    }
#pragma omp parallel sections
    {
#pragma omp section
      {
	treeMPI.commBodies();
	treeMPI.commCells();
      }
#pragma omp section
      {
	traversal.initListCount(cells);
	traversal.initWeight(cells);
	if (args.IneJ) {
	  traversal.traverse(cells, jcells, cycle, args.dual, false);
	} else {
	  traversal.traverse(cells, cells, cycle, args.dual, args.mutual);
	  jbodies = bodies;
	}
      }
    }
    if (baseMPI.mpisize > 1) {
      if (args.graft) {
	treeMPI.linkLET();
	gbodies = treeMPI.root2body();
	jcells = globalTree.buildTree(gbodies, buffer, globalBounds);
	treeMPI.attachRoot(jcells);
	traversal.traverse(cells, jcells, cycle, args.dual, false);
      } else {
	for (int irank=0; irank<baseMPI.mpisize; irank++) {
	  treeMPI.getLET(jcells, (baseMPI.mpirank+irank)%baseMPI.mpisize);
	  traversal.traverse(cells, jcells, cycle, args.dual, false);
	}
      }
    }
#else
    jbodies = bodies;
    for (int irank=0; irank<baseMPI.mpisize; irank++) {
      treeMPI.shiftBodies(jbodies);
      jcells.clear();
      localBounds = boundBox.getBounds(jbodies);
      jcells = localTree.buildTree(jbodies, buffer, localBounds);
      upDownPass.upwardPass(jcells);
      traversal.traverse(cells, jcells, cycle, args.dual, args.mutual);
    }
#endif
    upDownPass.downwardPass(cells);
    logger::stopPAPI();
    logger::stopTimer("Total FMM", 0);
#endif
    logger::printTitle("MPI direct sum");
    const int numTargets = 100;
    buffer = bodies;
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
    bodies = buffer;
    data.initTarget(bodies);
    logger::resetTimer("Total FMM");
    if (args.write) {
      logger::writeTime(baseMPI.mpirank);
    }
    traversal.writeList(cells, baseMPI.mpirank);
  }
  return 0;
}
