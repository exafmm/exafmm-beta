#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "ewald.h"
#include "kernel.h"
#include "logger.h"
#include "namespace.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
using namespace EXAFMM_NAMESPACE;

int main(int argc, char ** argv) {
  const int numBodies = 1000;
  const int ksize = 11;
  const vec3 cycle = 2 * M_PI;
  const real_t alpha = ksize / max(cycle);
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = max(cycle) / 2;
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);
  Args args(argc, argv);
  args.numBodies = numBodies;
  //args.images = 3;
  BaseMPI baseMPI;
  Bodies bodies, bodies2, jbodies, gbodies, buffer;
  BoundBox boundBox;
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit);
  BuildTree globalTree(1);
  Cells cells, jcells;
  Dataset data;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);
  Kernel kernel(args.P, eps2, wavek);
  Partition partition(baseMPI);
  Traversal traversal(kernel, args.theta, args.nspawn, args.images, args.path);
  TreeMPI treeMPI(kernel, baseMPI, args.theta, args.images);
  UpDownPass upDownPass(kernel);
  Verify verify(args.path);
  num_threads(args.threads);

  args.verbose &= baseMPI.mpirank == 0;
  verify.verbose = args.verbose;
  logger::verbose = args.verbose;
  logger::path = args.path;
  logger::printTitle("Ewald Parameters");
  args.print(logger::stringLength);
  ewald.print(logger::stringLength);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X *= cycle / (2 * M_PI);
  }
  buffer.reserve(bodies.size());
  bool pass = true;
  bool isTime = false;
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    int numIteration = 1;
    if (isTime) numIteration = 10;
    for (int it=0; it<numIteration; it++) {
      std::stringstream title;
      title << "Time average loop " << it;
      logger::printTitle(title.str());
      data.initTarget(bodies);
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
      traversal.traverse(cells, cells, cycle, args.dual);
      if (baseMPI.mpisize > 1) {
        if (args.graft) {
          treeMPI.linkLET();
          gbodies = treeMPI.root2body();
          jcells = globalTree.buildTree(gbodies, buffer, globalBounds);
          treeMPI.attachRoot(jcells);
          traversal.traverse(cells, jcells, cycle, args.dual);
        } else {
          for (int irank=0; irank<baseMPI.mpisize; irank++) {
            treeMPI.getLET(jcells, (baseMPI.mpirank+irank)%baseMPI.mpisize);
            traversal.traverse(cells, jcells, cycle, args.dual);
          }
        }
      }
      upDownPass.downwardPass(cells);
    }
    logger::stopPAPI();
    double totalFMM = logger::stopTimer("Total FMM");
    totalFMM /= numIteration;

    if (!isTime) {
      buffer = bodies;
#if 0 // 0: direct 1: Ewald
      vec3 localDipole = upDownPass.getDipole(bodies,0);
      vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
      int numBodies = baseMPI.allreduceInt(bodies.size());
      upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
      bodies2 = bodies;
      data.initTarget(bodies);
      logger::printTitle("Ewald Profiling");
      logger::startTimer("Total Ewald");
#if 1 // 0: gather bodies 1: shift bodies
      jbodies = bodies;
      for (int i=0; i<baseMPI.mpisize; i++) {
        if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << baseMPI.mpisize << std::endl;
        treeMPI.shiftBodies(jbodies);
        localBounds = boundBox.getBounds(jbodies);
        jcells = localTree.buildTree(jbodies, buffer, localBounds);
        ewald.wavePart(bodies, jbodies);
        ewald.realPart(cells, jcells);
      }
#else // gather bodies
      jbodies = treeMPI.allgatherBodies(bodies);
      jcells = localTree.buildTree(jbodies, buffer, globalBounds);
      ewald.wavePart(bodies, jbodies);
      ewald.realPart(cells, jcells);
#endif
      ewald.selfTerm(bodies);
      logger::printTitle("Total runtime");
      logger::printTime("Total FMM");
      logger::stopTimer("Total Ewald");
#else // direct
      jbodies = bodies;
      const int numTargets = 100;
      data.sampleBodies(bodies, numTargets);
      bodies2 = bodies;
      data.initTarget(bodies);
      logger::startTimer("Total Direct");
      for (int i=0; i<baseMPI.mpisize; i++) {
        treeMPI.shiftBodies(jbodies);
        traversal.direct(bodies, jbodies, cycle);
        if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
      }
      logger::printTitle("Total runtime");
      logger::printTime("Total FMM");
      logger::stopTimer("Total Direct");
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
      if (!baseMPI.mpirank) {
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, potRel, accRel);
      }
      MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
      if (pass) {
        if (verify.verbose) std::cout << "passed accuracy regression at t: " << t << std::endl;
        if (args.accuracy) break;
        t = -1;
        isTime = true;
      }
      bodies = buffer;
    } else {
      double totalFMMGlob;
      MPI_Reduce(&totalFMM, &totalFMMGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      totalFMMGlob /= baseMPI.mpisize;
      if (!baseMPI.mpirank) {
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, totalFMMGlob);
      }
      MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
      if (pass) {
        if (verify.verbose) std::cout << "passed time regression at t: " << t << std::endl;
        break;
      }
    }
    data.initTarget(bodies);
  }
  if (!pass) {
    if (verify.verbose) {
      if(!isTime) std::cout << "failed accuracy regression" << std::endl;
      else std::cout << "failed time regression" << std::endl;
    }
    abort();
  }
  return 0;
}
