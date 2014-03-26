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
  BoundBox boundBox(args.nspawn);
  BuildTree buildTree(args.ncrit, args.nspawn);
  Dataset data;
  Logger logger;
  Partition partition;
  Traversal traversal(args.nspawn, args.images);
  TreeMPI treeMPI(args.images);
  UpDownPass upDownPass(args.theta);
  Verify verify;

  const int ksize = 11;
  const real_t cycle = 2 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = cycle * alpha / 3;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);

  args.verbose &= baseMPI.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundBox.verbose = true;
    buildTree.verbose = true;
    upDownPass.verbose = true;
    traversal.verbose = true;
    ewald.verbose = true;
    treeMPI.verbose = true;
    verify.verbose = true;
  }
  logger.printTitle("Ewald Parameters");
  args.print(logger.stringLength, P);
  ewald.print(logger.stringLength);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  //data.writeSources(bodies, baseMPI.mpirank);
  Bounds localBounds = boundBox.getBounds(bodies);
  Bounds globalBounds = baseMPI.allreduceBounds(localBounds);
  localBounds = partition.octsection(bodies, globalBounds);
  bodies = treeMPI.commBodies(bodies);

  Cells cells = buildTree.buildTree(bodies, localBounds);
  upDownPass.upwardPass(cells);
  treeMPI.setLET(cells,localBounds,cycle);
  treeMPI.commBodies();
  treeMPI.commCells();

  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Cells jcells;
  for (int irank=1; irank<baseMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(baseMPI.mpirank+irank)%baseMPI.mpisize);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  upDownPass.downwardPass(cells);
  vec3 localDipole = upDownPass.getDipole(bodies,0);
  vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
  int numBodies = baseMPI.allreduceInt(bodies.size());
  upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
#if 1
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.printTitle("Ewald Profiling");
  logger.startTimer("Total Ewald");
#if 1
  Bodies jbodies = bodies;
  for (int i=0; i<baseMPI.mpisize; i++) {
    if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << baseMPI.mpisize << std::endl;
    treeMPI.shiftBodies(jbodies);
    localBounds = boundBox.getBounds(jbodies);
    Cells jcells = buildTree.buildTree(jbodies, localBounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
  }
#else
  Bodies jbodies = treeMPI.allgatherBodies(bodies);
  jcells = buildTree.buildTree(jbodies, globalBounds);
  ewald.wavePart(bodies, jbodies);
  ewald.realPart(cells, jcells);
#endif
  ewald.selfTerm(bodies);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Ewald");
#else
  Bodies jbodies = bodies;
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  for (int i=0; i<baseMPI.mpisize; i++) {
    treeMPI.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
  }
  traversal.normalize(bodies);
  upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
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
  logger.printTitle("FMM vs. Ewald");
  double potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
  double potNrmGlob = potSumGlob * potSumGlob;
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
  buildTree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();
#if VTK
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) B->IBODY = 0;
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    Body body;
    body.IBODY = 1;
    body.X     = C->X;
    body.SRC   = 0;
    bodies.push_back(body);
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(bodies);
  vtk.plot();
#endif
  return 0;
}
