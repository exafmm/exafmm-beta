#include "tree_mpi.h"
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
  BuildTree build(args.ncrit, args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn, args.images);
  TreeMPI treeMPI(args.images);
  args.numBodies /= treeMPI.mpisize;
  args.verbose &= treeMPI.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    build.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    treeMPI.verbose = true;
    verify.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength, P, treeMPI.mpirank);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, treeMPI.mpirank, treeMPI.mpisize);
  Bounds localBounds = boundbox.getBounds(bodies);
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, treeMPI.mpirank+treeMPI.mpisize, treeMPI.mpisize);
  localBounds = boundbox.getBounds(jbodies,localBounds);
  Bounds globalBounds = treeMPI.allreduceBounds(localBounds);
  localBounds = treeMPI.partition(bodies,globalBounds);
  bodies = treeMPI.commBodies(bodies);
  treeMPI.partition(jbodies,globalBounds);
  jbodies = treeMPI.commBodies(jbodies);
  Cells cells = build.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
  Cells jcells = build.buildTree(jbodies, localBounds);
  pass.upwardPass(jcells);

  treeMPI.setLET(jcells,localBounds,cycle);
  treeMPI.commBodies();
  treeMPI.commCells();
  traversal.dualTreeTraversal(cells, jcells, cycle);
  for (int irank=1; irank<treeMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(treeMPI.mpirank+irank)%treeMPI.mpisize);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  pass.downwardPass(cells);

  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  logger.printTitle("MPI direct sum");
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  for (int i=0; i<treeMPI.mpisize; i++) {
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << treeMPI.mpisize << std::endl;
    treeMPI.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
  }
  traversal.normalize(bodies);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
  boundbox.writeTime(treeMPI.mpirank);
  build.writeTime(treeMPI.mpirank);
  pass.writeTime(treeMPI.mpirank);
  traversal.writeTime(treeMPI.mpirank);
  treeMPI.writeTime(treeMPI.mpirank);
  double potDif = verify.getDifScalar(bodies, bodies2);
  double potNrm = verify.getNrmScalar(bodies);
  double accDif = verify.getDifVector(bodies, bodies2);
  double accNrm = verify.getNrmVector(bodies);
  double potDifGlob, potNrmGlob, accDifGlob, accNrmGlob;
  MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accDif, &accDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accNrm, &accNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  logger.printTitle("FMM vs. direct");
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
  build.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();

#if VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->IBODY = 0;
  for (int irank=0; irank<treeMPI.mpisize; irank++) {
    treeMPI.gettreeMPI(jcells,(treeMPI.mpirank+irank)%treeMPI.mpisize);
    for (C_iter C=jcells.begin(); C!=jcells.end(); C++) {
      Body body;
      body.IBODY = 1;
      body.X     = C->X;
      body.SRC   = 0;
      jbodies.push_back(body);
    }
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(jbodies);
  for (int i=1; i<treeMPI.mpisize; i++) {
    treeMPI.shiftBodies(jbodies);
    vtk.setGroupOfPoints(jbodies);
  }
  if (treeMPI.mpirank == 0) {
    vtk.plot();
  }
#endif
  return 0;
}
