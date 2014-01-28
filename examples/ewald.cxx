#include "tree_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "ewald.h"
#include "logger.h"
#include "traversal.h"
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
  Dataset data;
  Logger logger;
  Verify verify;

  const int ksize = 11;
  const real_t cycle = 2 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = cycle * alpha / 3;
  BoundBox boundbox(args.nspawn);
  BuildTree build(args.ncrit, args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn, args.images);
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);
  TreeMPI treeMPI(args.images);
  args.verbose &= treeMPI.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    build.verbose = true;
    pass.verbose = true;
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
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, treeMPI.mpirank, treeMPI.mpisize);
  //data.writeSources(bodies, treeMPI.mpirank);
  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = treeMPI.allreduceBounds(localBounds);
  localBounds = treeMPI.partition(bodies, globalBounds);
  bodies = treeMPI.commBodies(bodies);

  Cells cells = build.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
  treeMPI.setLET(cells,localBounds,cycle);
  treeMPI.commBodies();
  treeMPI.commCells();

  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Cells jcells;
  for (int irank=1; irank<treeMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(treeMPI.mpirank+irank)%treeMPI.mpisize);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  pass.downwardPass(cells);
  vec3 localDipole = pass.getDipole(bodies,0);
  vec3 globalDipole = treeMPI.allreduceVec3(localDipole);
  int numBodies = treeMPI.allreduceInt(bodies.size());
  pass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
#if 1
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.printTitle("Ewald Profiling");
  logger.startTimer("Total Ewald");
#if 1
  Bodies jbodies = bodies;
  for (int i=0; i<treeMPI.mpisize; i++) {
    if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << treeMPI.mpisize << std::endl;
    treeMPI.shiftBodies(jbodies);
    localBounds = boundbox.getBounds(jbodies);
    Cells jcells = build.buildTree(jbodies, localBounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
  }
#else
  Bodies jbodies = treeMPI.allgatherBodies(bodies);
  jcells = build.buildTree(jbodies, globalBounds);
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
  for (int i=0; i<treeMPI.mpisize; i++) {
    treeMPI.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << treeMPI.mpisize << std::endl;
  }
  traversal.normalize(bodies);
  pass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
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
  build.printTreeData(cells);
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
