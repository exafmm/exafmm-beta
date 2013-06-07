#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "ewald.h"
#include "logger.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#if VTK
#include "vtk.h"
#endif
#if MASS
#error Turn off MASS for this test
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;
  Sort sort;

  const int ksize = 11.;
  const real_t cycle = 2 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = cycle * alpha / 3;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
  Ewald ewald(ksize,alpha,sigma,cutoff,cycle);
  LocalEssentialTree LET(args.images);
  args.verbose &= LET.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    ewald.verbose = true;
    LET.verbose = true;
  }
  logger.printTitle("Ewald Parameters");
  args.print(logger.stringLength,P);
  ewald.print(logger.stringLength);
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution);
  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies, globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);

  Cells cells = tree.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
  LET.setLET(cells,localBounds,cycle);
  LET.commBodies();
  LET.commCells();

  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Cells jcells;
  for (int irank=1; irank<LET.mpisize; irank++) {
    LET.getLET(jcells,(LET.mpirank+irank)%LET.mpisize);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  pass.downwardPass(cells);
  vec3 localDipole = pass.getDipole(bodies,0);
  vec3 globalDipole = LET.allreduceVec3(localDipole);
  int numBodies = LET.allreduceInt(bodies.size());
  pass.dipoleCorrection(bodies,globalDipole,numBodies,cycle);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
#if 1
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.printTitle("Ewald Profiling");
  logger.startTimer("Total Ewald");
#if 1
  Bodies jbodies = bodies;
  for (int i=0; i<LET.mpisize; i++) {
    LET.shiftBodies(jbodies);
    if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << LET.mpisize << std::endl;
    localBounds = boundbox.getBounds(jbodies);
    Cells jcells = tree.buildTree(jbodies, localBounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
  }
#else
  Bodies jbodies = LET.allgatherBodies(bodies);
  jcells = tree.buildTree(jbodies, globalBounds);
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
  for (int i=0; i<LET.mpisize; i++) {
    LET.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << LET.mpisize << std::endl;
  }
  traversal.normalize(bodies);
  pass.dipoleCorrection(bodies,globalDipole,numBodies,cycle);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
#endif
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  ewald.evalError(bodies2, bodies, diff1, norm1, diff2, norm2);
  MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  logger.printTitle("FMM vs. Ewald");
  logger.printError(diff3, norm3, diff4, norm4);
  tree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();
#if VTK
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) B->ICELL = 0;
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    bodies.push_back(body);
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(bodies);
  vtk.plot();
#endif
  MPI_Finalize();
  return 0;
}
