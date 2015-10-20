#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;

  const real_t eps2 = 0.0;
  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta,eps2);
  Traversal traversal(args.nspawn,args.images,eps2);
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength,P);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, 0);
  Bounds bounds = boundbox.getBounds(bodies);
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, 1);
  bounds = boundbox.getBounds(jbodies,bounds);
#endif
  Cells cells = tree.buildTree(bodies, bounds);
  pass.upwardPass(cells);
#if IneJ
  Cells jcells = tree.buildTree(jbodies, bounds);
  pass.upwardPass(jcells);
  traversal.dualTreeTraversal(cells, jcells, cycle);
#else
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Bodies jbodies = bodies;
#endif
  pass.downwardPass(cells);
  logger.printTitle("Total runtime");
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
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
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  traversal.direct(bodies, jbodies, cycle);
  traversal.normalize(bodies);
  logger.stopTimer("Total Direct");
  double diff1 = 0, norm1 = 0;
  data.evalError(bodies2, bodies, diff1, norm1);
  logger.printTitle("FMM vs. direct");
  logger.printError(diff1, norm1);
  tree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();
  return 0;
}
