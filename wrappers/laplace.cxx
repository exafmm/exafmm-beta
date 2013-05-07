#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"

extern "C" void FMM(int ni, double * xi, double * pi, double * fi, int nj, double * xj, double *qj, int periodicflag) {
  Args args;
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Logger logger;
  Sort sort;

  args.numBodies = ni;
  args.THETA = 0.6;
  args.NCRIT = 16;
  args.NSPAWN = 1000;
  args.IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  args.mutual = 0;
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
  logger.startTimer("Total FMM");
  bodies.resize(ni);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    B->SRC  = 1;
    B->TRG[0] = pi[i];
    B->TRG[1] = fi[3*i+0];
    B->TRG[2] = fi[3*i+1];
    B->TRG[3] = fi[3*i+2];
    B->IBODY  = i;
  }
  Bounds localBounds = boundbox.getBounds(bodies);
  jbodies.resize(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0] = xj[3*i+0];
    B->X[1] = xj[3*i+1];
    B->X[2] = xj[3*i+2];
    B->SRC  = qj[i];
  }
  localBounds = boundbox.getBounds(jbodies,localBounds);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  LET.partition(jbodies,globalBounds);
  jbodies = sort.sortBodies(jbodies);
  jbodies = LET.commBodies(jbodies);
  Box box = boundbox.bounds2box(localBounds);
  logger.startPAPI();
  tree.buildTree(bodies, cells, box);
  pass.upwardPass(cells);
  tree.buildTree(jbodies, jcells, box);
  pass.upwardPass(jcells);
  LET.setLET(jcells,localBounds,cycle);
  LET.commBodies();
  LET.commCells();
  traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
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
    pi[i]     = B->TRG[0];
    fi[3*i+0] = B->TRG[1];
    fi[3*i+1] = B->TRG[2];
    fi[3*i+2] = B->TRG[3];
  }
}
