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

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Cells cells;
  Dataset data;
  Logger logger;
  Sort sort;

  const real_t ksize = 11.;
  const real_t alpha = .1;
  const real_t sigma = .25 / M_PI;
  const real_t theta = .5;
  const real_t cycle = 100.;
  BoundBox boundbox(args.NSPAWN);
  BuildTree tree(args.NCRIT,args.NSPAWN);
  UpDownPass pass(args.THETA);
  Traversal traversal(args.NSPAWN,args.IMAGES);
  Ewald ewald(ksize,alpha,sigma,theta,cycle);
  logger.verbose = true;
  if (args.verbose) {
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
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
  Bodies bodies = data.initBodies(bodies, args.distribution);
  Bounds bounds = ewald.rescale(bodies);
  Box box = boundbox.bounds2box(bounds);
  tree.buildTree(bodies, cells, box);
  pass.upwardPass(cells);
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  pass.downwardPass(cells);
  if (args.verbose) logger.printTitle("Total runtime");
  logger.stopPAPI();
  logger.stopTimer("Total FMM", logger.verbose);
  boundbox.writeTime();
  tree.writeTime();
  pass.writeTime();
  traversal.writeTime();
  boundbox.resetTimer();
  tree.resetTimer();
  pass.resetTimer();
  traversal.resetTimer();
  logger.resetTimer();
#if 1
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Ewald");
  ewald.wavePart(bodies);
  ewald.realPart(cells,cells);
  logger.stopTimer("Total Ewald", args.verbose);
#else
  Bodies jbodies = bodies;
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  traversal.direct(bodies, jbodies, cycle);
  traversal.normalize(bodies);
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
