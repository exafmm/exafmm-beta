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
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution);
  Bounds bounds = boundbox.getBounds(bodies);
  Cells cells = tree.buildTree(bodies, bounds);
  pass.upwardPass(cells);
#if SPLIT
  Bodies pbodies = data.getPositive(bodies);
  Bodies nbodies = data.getNegative(bodies);
  Bounds pbounds = boundbox.getBounds(pbodies);
  Bounds nbounds = boundbox.getBounds(nbodies);
  Cells pcells = tree.buildTree(pbodies, pbounds);
  Cells ncells = tree.buildTree(nbodies, nbounds);
  pass.upwardPass(pcells);
  pass.upwardPass(ncells);
#endif
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution);
  bounds = boundbox.getBounds(jbodies);
  Cells jcells = tree.buildTree(jbodies, bounds);
  pass.upwardPass(jcells);
  traversal.dualTreeTraversal(cells, jcells, cycle);
#else
#if SPLIT
  traversal.dualTreeTraversal(pcells, pcells, cycle, args.mutual);
  traversal.dualTreeTraversal(pcells, ncells, cycle);
  traversal.dualTreeTraversal(ncells, pcells, cycle);
  traversal.dualTreeTraversal(ncells, ncells, cycle, args.mutual);
#else
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
#endif
#endif
#if SPLIT
  pass.downwardPass(pcells);
  pass.downwardPass(ncells);
  bodies = pbodies;
  bodies.insert(bodies.end(),nbodies.begin(),nbodies.end());
#else
  pass.downwardPass(cells);
#endif  
  Bodies jbodies = bodies;
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
