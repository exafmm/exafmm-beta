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

  const real_t cycle = 2 * M_PI;
  const real_t ksize = 11.;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t theta = .4;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
  Ewald ewald(ksize,alpha,sigma,theta,cycle);
  logger.verbose = true;
  if (args.verbose) {
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    ewald.verbose = true;
    logger.printTitle("Parameters");
  }
  args.print(logger.stringLength,P);
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  if(args.verbose) logger.printTitle("Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution);
  Bounds bounds = boundbox.getBounds(bodies);
  Cells cells = tree.buildTree(bodies, bounds);
  pass.upwardPass(cells);
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  pass.downwardPass(cells);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
#if 1
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Ewald");
  ewald.wavePart(bodies);
  ewald.realPart(cells,cells);
  ewald.dipoleCorrection(bodies,0,cycle);
  if (args.verbose) logger.printTitle("Total runtime");
  if (logger.verbose) logger.printTime("Total FMM");
  logger.stopTimer("Total Ewald", args.verbose);
#else
  Bodies jbodies = bodies;
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  traversal.direct(bodies, jbodies, cycle);
  traversal.normalize(bodies);
  if (args.verbose) logger.printTitle("Total runtime");
  if (logger.verbose) logger.printTime("Total FMM");
  logger.stopTimer("Total Direct", args.verbose);
#endif
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  ewald.evalError(bodies2, bodies, diff1, norm1, diff2, norm2);
  if (args.verbose) {
    logger.printTitle("FMM vs. Ewald");
    logger.printError(diff1, norm1, diff2, norm2);
    tree.printTreeData(cells);
    traversal.printTraversalData();
    logger.printPAPI();
  }
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
  return 0;
}
