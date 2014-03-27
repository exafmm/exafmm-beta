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
  BoundBox boundBox(args.nspawn);
  BuildTree buildTree(args.ncrit,args.nspawn);
  Dataset data;
  Logger logger;
  Traversal traversal(args.nspawn,args.images);
  UpDownPass upDownPass(args.theta);
  Verify verify;

  const real_t cycle = 2 * M_PI;
  if (args.verbose) {
    logger.verbose = true;
    boundBox.verbose = true;
    buildTree.verbose = true;
    traversal.verbose = true;
    upDownPass.verbose = true;
    verify.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength, P, 0);
  for (int t = 0; t < args.repeat; t++) {
#if DAG_RECORDER == 2
    dr_start(0);
#endif
    logger.printTitle("FMM Profiling");
    logger.startTimer("Total FMM");
    logger.startPAPI();
    Bodies bodies = data.initBodies(args.numBodies, args.distribution, 0);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    Bounds bounds = boundBox.getBounds(bodies);
#if IneJ
    Bodies jbodies = data.initBodies(args.numBodies, args.distribution, 1);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
    bounds = boundBox.getBounds(jbodies,bounds);
#endif
    Cells cells = buildTree.buildTree(bodies, bounds);
    upDownPass.upwardPass(cells);
#if IneJ
    Cells jcells = buildTree.buildTree(jbodies, bounds);
    upDownPass.upwardPass(jcells);
    traversal.dualTreeTraversal(cells, jcells, cycle);
#else
    traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
    Bodies jbodies = bodies;
#endif
    upDownPass.downwardPass(cells);
    logger.printTitle("Total runtime");
    logger.stopPAPI();
    logger.stopTimer("Total FMM");
#if WRITE_TIME
    boundBox.writeTime();
    buildTree.writeTime();
    upDownPass.writeTime();
    traversal.writeTime();
#endif
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
    buildTree.printTreeData(cells);
    traversal.printTraversalData();
    logger.printPAPI();
#if DAG_RECORDER == 2
    dr_stop();
#endif
  }
#if DAG_RECORDER == 2
  dr_dump();
#endif
#if VTK
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) B->IBODY = 0;
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    Body body;
    body.IBODY = 1;
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
