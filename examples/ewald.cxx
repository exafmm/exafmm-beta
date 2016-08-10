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
#if EXAFMM_MASS
#error Turn off MASS for this test
#endif
#if EXAFMM_EXPANSION < 10
#error Use P >=10 for this test
#endif
using namespace exafmm;

int main(int argc, char ** argv) {
  const int ksize = 11;
  const vec3 cycle = 2 * M_PI;
  const real_t alpha = 10 / max(cycle);
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = max(cycle) / 2;
  Args args(argc, argv);
  args.numBodies = 1000;
  args.images = 3;
  BaseMPI baseMPI;
  Bodies bodies, bodies2, jbodies, gbodies, buffer;
  BoundBox boundBox(args.nspawn);
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit, args.nspawn);
  BuildTree globalTree(1, args.nspawn);
  Cells cells, jcells;
  Dataset data;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);
  Partition partition(baseMPI.mpirank, baseMPI.mpisize);
  Traversal traversal(args.nspawn, args.images, args.path);
  TreeMPI treeMPI(baseMPI.mpirank, baseMPI.mpisize, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify(args.path);
  num_threads(args.threads);

  kernel::eps2 = 0.0;
  args.verbose &= baseMPI.mpirank == 0;
  verify.verbose = args.verbose;
  logger::verbose = args.verbose;
  logger::path = args.path;
  logger::printTitle("Ewald Parameters");
  args.print(logger::stringLength, P);
  ewald.print(logger::stringLength);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X *= cycle / (2 * M_PI);
  }
  buffer.reserve(bodies.size());
  //data.writeSources(bodies, baseMPI.mpirank);
  bool pass = true;
  bool isTime = false;
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    localBounds = boundBox.getBounds(bodies);
    globalBounds = baseMPI.allreduceBounds(localBounds);
    localBounds = partition.octsection(bodies, globalBounds);
    bodies = treeMPI.commBodies(bodies);

    cells = localTree.buildTree(bodies, buffer, localBounds);
    upDownPass.upwardPass(cells);
    treeMPI.allgatherBounds(localBounds);
    treeMPI.setLET(cells, cycle);
    treeMPI.commBodies();
    treeMPI.commCells();

    traversal.initListCount(cells);
    traversal.initWeight(cells);
    traversal.traverse(cells, cells, cycle, args.dual, args.mutual);
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
    upDownPass.downwardPass(cells);
    vec3 localDipole = upDownPass.getDipole(bodies,0);
    vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
    int numBodies = baseMPI.allreduceInt(bodies.size());
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
    logger::stopPAPI();
    double totalFMM = logger::stopTimer("Total FMM");
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
      jcells = localTree.buildTree(jbodies, buffer, localBounds);
      ewald.wavePart(bodies, jbodies);
      ewald.realPart(cells, jcells);
    }
#else
    jbodies = treeMPI.allgatherBodies(bodies);
    jcells = localTree.buildTree(jbodies, buffer, globalBounds);
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
    buffer = bodies;
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
    bodies = buffer;
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
    double potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
    double potNrmGlob = potSumGlob * potSumGlob;
    double potRel = std::sqrt(potDifGlob/potNrmGlob);
    double accRel = std::sqrt(accDifGlob/accNrmGlob);
    logger::printTitle("FMM vs. Ewald");
    verify.print("Rel. L2 Error (pot)",potRel);
    verify.print("Rel. L2 Error (acc)",accRel);
    localTree.printTreeData(cells);
    traversal.printTraversalData();
    logger::printPAPI();
    data.initTarget(bodies);
    double totalFMMGlob;
    MPI_Reduce(&totalFMM, &totalFMMGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    totalFMMGlob /= baseMPI.mpisize;
    if (!baseMPI.mpirank) {
      if (!isTime)
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, potRel, accRel);
      else
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, totalFMMGlob);
    }
    MPI_Bcast(&pass, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (pass) {
      if (!isTime) {
        if (args.verbose) std::cout << "passed accuracy regression at t: " << t << std::endl; 
        t = -1;
        isTime = true;        
      } else {
        if (args.verbose) std::cout << "passed time regression at t: " << t << std::endl;
        break;
      }
    }
  }
  if (!pass) {
    if (args.verbose) {
      if(!isTime) std::cout << "failed accuracy regression" << std::endl;
      else std::cout << "failed time regression" << std::endl;
    }
    abort();
  }
  return 0;
}
