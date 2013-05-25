#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"

extern "C" void fmm(int n, double * x, double * q, double * p, double * f, double cycle, int images) {
  Args args;
  Logger logger;
  Sort sort;

  args.numBodies = n;
  args.theta = 0.35;
  args.ncrit = 16;
  args.nspawn = 1000;
  args.images = images;
  args.mutual = 0;
  args.verbose = 1;
  args.distribution = "external";

  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
  LocalEssentialTree LET(args.images);
  logger.verbose = LET.mpirank == 0;
  args.verbose &= logger.verbose;
  if (args.verbose) {
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    LET.verbose = true;
    logger.printTitle("FMM Parameters");
  }
  if(LET.mpirank == 0) args.print(logger.stringLength,P);
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  if (args.verbose) logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    if( B->X[0] < -cycle/2 ) B->X[0] += cycle;
    if( B->X[1] < -cycle/2 ) B->X[1] += cycle;
    if( B->X[2] < -cycle/2 ) B->X[2] += cycle;
    if( B->X[0] >  cycle/2 ) B->X[0] -= cycle;
    if( B->X[1] >  cycle/2 ) B->X[1] -= cycle;
    if( B->X[2] >  cycle/2 ) B->X[2] -= cycle;
    B->SRC = q[i];
    B->TRG[0] = p[i];
    B->TRG[1] = f[3*i+0];
    B->TRG[2] = f[3*i+1];
    B->TRG[3] = f[3*i+2];
    B->IBODY = i;
  }
  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
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
  vec3 globalDipole = LET.allreduce(localDipole);
  int numBodies = LET.allreduce(bodies.size());
  pass.dipoleCorrection(bodies,globalDipole,numBodies,cycle);

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
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
}
