#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset data;
  Logger logger;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.NSPAWN);
  BuildTree tree(args.NCRIT,args.NSPAWN);
  UpDownPass pass(args.IMAGES,args.THETA);
  Traversal traversal(args.NSPAWN,args.IMAGES);
  logger.printNow = true;
  boundbox.printNow = true;
  tree.printNow = true;
  pass.printNow = true;
  traversal.printNow = true;
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
#ifdef MANY
  for (int it=0; it<25; it++) {
    int numBodies = int(pow(10,(it+24)/8.0));
#else
  {
    int numBodies = args.numBodies;
#endif // MANY
    if (logger.printNow) std::cout << std::endl
      << "Num bodies           : " << numBodies << std::endl;
    std::cout << "--- Profiling --------------------" << std::endl;
    bodies.resize(numBodies);
    data.initBodies(bodies, args.distribution);
    Bounds bounds = boundbox.getBounds(bodies);
    Box box = boundbox.bounds2box(bounds);
    tree.buildTree(bodies, cells, box);                         //TODO : make it work without this
    tree.resetTimer();
    logger.startTimer("Total FMM");
    tree.buildTree(bodies, cells, box);
    pass.upwardPass(cells);
    logger.startPAPI();
    traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
    logger.stopPAPI();
    pass.downwardPass(cells);
    std::cout << "--- Total runtime ----------------" << std::endl;
    logger.stopTimer("Total FMM", logger.printNow);
    boundbox.writeTime();
    tree.writeTime();
    pass.writeTime();
    traversal.writeTime();
    boundbox.resetTimer();
    tree.resetTimer();
    pass.resetTimer();
    traversal.resetTimer();
    jbodies = bodies;
    if (int(bodies.size()) > args.numTarget) data.sampleBodies(bodies, args.numTarget);
    Bodies bodies2 = bodies;
    data.initTarget(bodies2);
    logger.startTimer("Total Direct");
    pass.direct(bodies2, jbodies, cycle);
    pass.normalize(bodies2);
    logger.stopTimer("Total Direct", logger.printNow);
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    data.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    if (logger.printNow) {
      data.printError(diff1, norm1, diff2, norm2);
      tree.printTreeData(cells);
      traversal.printTraversalData();
    }
  }
#ifdef VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->ICELL = 0;
  for (C_iter C=jcells.begin(); C!=jcells.end(); C++) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(M_PI,0);
  vtk.setGroupOfPoints(jbodies, Ncell);
  vtk.plot(Ncell);
#endif
  return 0;
}
