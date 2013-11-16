#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "verify.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;
  Verify verify;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree build(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    build.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    verify.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength, P, 0);
  for (int t = 0; t < args.repeat; t++) {
    logger.printTitle("FMM Profiling");
    logger.startTimer("Total FMM");
    logger.startPAPI();
    Bodies bodies = data.initBodies(args.numBodies, args.distribution, 0);
    Bounds bounds = boundbox.getBounds(bodies);
#if IneJ
    Bodies jbodies = data.initBodies(args.numBodies, args.distribution, 1);
    bounds = boundbox.getBounds(jbodies,bounds);
#endif
    Cells cells = build.buildTree(bodies, bounds);
    pass.upwardPass(cells);
#if IneJ
    Cells jcells = build.buildTree(jbodies, bounds);
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
    build.writeTime();
    pass.writeTime();
    traversal.writeTime();
    data.sampleBodies(bodies, args.numTargets);
    Bodies bodies2 = bodies;
    data.initTarget(bodies);
    logger.startTimer("Total Direct");
    traversal.direct(bodies, jbodies, cycle);
    traversal.normalize(bodies);
    logger.stopTimer("Total Direct");
    double potDif = verify.getDifScalar(bodies, bodies2);
    double potNrm = verify.getNrmScalar(bodies);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    logger.printTitle("FMM vs. direct");
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
    build.printTreeData(cells);
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
