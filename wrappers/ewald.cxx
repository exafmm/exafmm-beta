#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "ewald.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"

extern "C" void ewald(int n, double * x, double * q, double * p, double * f,
                      int ksize, double alpha, double cycle) {
  Args args;
  Logger logger;
  Sort sort;

  args.numBodies = n;
  args.theta = 0.35;
  args.ncrit = 32;
  args.nspawn = 1000;
  args.images = 0;
  args.mutual = 0;
  args.verbose = 1;
  args.distribution = "external";

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
  logger.printTitle("Ewald Profiling");
  logger.startTimer("Total Ewald");
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

  Cells cells = tree.buildTree(bodies, globalBounds);
  Bodies jbodies = bodies;

  for (int i=0; i<LET.mpisize; i++) {
    LET.shiftBodies(jbodies);
    Cells jcells = tree.buildTree(jbodies, globalBounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
    if (args.verbose) std::cout << "Ewald loop           : " << i+1 << "/" << LET.mpisize << std::endl;
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->ICELL = B->IBODY;
  }
  bodies = sort.sortBodies(bodies);
  logger.stopPAPI();
  logger.stopTimer("Total Ewald");
  if (logger.verbose) {
    logger.printTitle("Total runtime");
    logger.printTime("Total Ewald");
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
}

extern "C" void ewald_(int * n, double * x, double * q, double * p, double * f,
                      int * ksize, double * alpha, double * cycle) {
  ewald(*n, x, q, p, f, *ksize, *alpha, *cycle);
}
