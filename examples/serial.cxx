#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
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
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, args.chargeSign);
  Bounds bounds = boundbox.getBounds(bodies);
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, args.chargeSign);
  bounds = boundbox.getBounds(jbodies,bounds);
#endif
  Box box = boundbox.bounds2box(bounds);
  Cells cells = tree.buildTree(bodies, box);                    // TODO : make it work without this
#if IneJ
  Cells jcells = tree.buildTree(jbodies, box);                  // TODO : make it work without this
#endif
  tree.resetTimer();
  logger.startTimer("Total FMM");
  logger.startPAPI();
  cells = tree.buildTree(bodies, box);
  pass.upwardPass(cells);
#if IneJ
  jcells = tree.buildTree(jbodies, box);
  pass.upwardPass(jcells);
  traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
#else
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Bodies jbodies = bodies;
  Cells jcells = cells;
#endif
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
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies2);
  logger.startTimer("Total Direct");
  traversal.direct(bodies2, jbodies, cycle);
  traversal.normalize(bodies2);
  logger.stopTimer("Total Direct", args.verbose);
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  data.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
  if (args.verbose) {
    logger.printTitle("FMM vs. direct");
    logger.printError(diff1, norm1, diff2, norm2);
    tree.printTreeData(cells);
    traversal.printTraversalData();
    logger.printPAPI();
  }
#if VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->ICELL = 0;
  for (C_iter C=jcells.begin(); C!=jcells.end(); C++) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(jbodies);
  vtk.plot();
#endif
  return 0;
}
