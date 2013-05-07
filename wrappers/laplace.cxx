#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"

extern "C" void FMM(int n, double* x, double* q, double *p, double* f, int periodicflag) {
  Args args;
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Logger logger;
  Sort sort;

  args.numBodies = n;
  args.THETA = 0.6;
  args.NCRIT = 16;
  args.NSPAWN = 1000;
  args.IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  args.mutual = 1;
  args.verbose = 1;
  args.distribution = "external";

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.NSPAWN);
  BuildTree tree(args.NCRIT,args.NSPAWN);
  UpDownPass pass(args.THETA);
  Traversal traversal(args.NSPAWN,args.IMAGES);
  LocalEssentialTree LET(args.IMAGES);
  logger.verbose = LET.MPIRANK == 0;
  args.verbose &= logger.verbose;
  if (args.verbose) {
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    LET.verbose = true;
    logger.printTitle("Parameters");
  }
  if(LET.MPIRANK == 0) args.print(logger.stringLength,P);
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  if (args.verbose) logger.printTitle("Profiling");
  bodies.resize(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    B->SRC  = q[i];
    B->TRG[0] =  p[i];
    B->TRG[1] = -f[3*i+0];
    B->TRG[2] = -f[3*i+1];
    B->TRG[3] = -f[3*i+2];
    B->IBODY  = i;
  }
  logger.startTimer("Total FMM");

  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  Box box = boundbox.bounds2box(localBounds);
  logger.startPAPI();
  tree.buildTree(bodies, cells, box);
  pass.upwardPass(cells);
  LET.setLET(cells,localBounds,cycle);
  LET.commBodies();
  LET.commCells();
  traversal.dualTreeTraversal(cells, cells, cycle);
  for (int irank=1; irank<LET.MPISIZE; irank++) {
    LET.getLET(jcells,(LET.MPIRANK+irank)%LET.MPISIZE);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  pass.downwardPass(cells);

  LET.unpartition(bodies);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->ICELL = B->IBODY;
  }
  bodies = sort.sortBodies(bodies);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  if (logger.verbose) {
    logger.printTitle("Total runtime");
    logger.printTime("Total FMM");
  }

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i]     = B->SRC;
    p[i]     =  B->TRG[0];
    f[3*i+0] = -B->TRG[1];
    f[3*i+1] = -B->TRG[2];
    f[3*i+2] = -B->TRG[3];
  }
}
