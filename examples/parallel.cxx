#include "tree_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "sort.h"
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
  Sort sort;
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
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, treeMPI.mpirank+treeMPI.mpisize, treeMPI.mpisize);
  localBounds = boundbox.getBounds(jbodies,localBounds);
#endif
  Bounds globalBounds = treeMPI.allreduceBounds(localBounds);
  localBounds = treeMPI.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = treeMPI.commBodies(bodies);
#if IneJ
  treeMPI.partition(jbodies,globalBounds);
  jbodies = sort.sortBodies(jbodies);
  jbodies = treeMPI.commBodies(jbodies);
#endif
  Cells cells = build.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
#if IneJ
  Cells jcells = build.buildTree(jbodies, localBounds);
  pass.upwardPass(jcells);
#endif

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree : Step 1
#if IneJ
  treeMPI.setLET(jcells,localBounds,cycle);
#else
  treeMPI.setLET(cells,localBounds,cycle);
#endif
  treeMPI.commBodies();
  treeMPI.commCells();
#if IneJ
  traversal.dualTreeTraversal(cells, jcells, cycle);
#else
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Bodies jbodies = bodies;
  Cells jcells;
#endif
  for (int irank=1; irank<treeMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(treeMPI.mpirank+irank)%treeMPI.mpisize);

#if 0 // Set to 1 for debugging full treeMPI communication : Step 2 (treeMPI must be set to full tree)
    treeMPI.shiftBodies(jbodies); // This will overwrite recvBodies. (define recvBodies2 in body_mpi.h to avoid this)
    Cells icells;
    build.buildTree(jbodies, icells);
    pass.upwardPass(icells);
    assert( icells.size() == jcells.size() );
    CellQueue Qi, Qj;
    Qi.push(icells.begin());
    Qj.push(jcells.begin());
    int ic=0;
    while (!Qi.empty()) {
      C_iter Ci=Qi.front(); Qi.pop();
      C_iter Cj=Qj.front(); Qj.pop();
      if (Ci->ICELL != Cj->ICELL) {
	std::cout << treeMPI.mpirank << " ICELL  : " << Ci->ICELL << " " << Cj->ICELL << std::endl;
	break;
      }
      if (Ci->NCHILD != Cj->NCHILD) {
	std::cout << treeMPI.mpirank << " NCHILD : " << Ci->NCHILD << " " << Cj->NCHILD << std::endl;
	break;
      }
      if (Ci->NCHILD == 0 && Cj->NCHILD == 0 && Ci->NBODY != Cj->NBODY) {
	std::cout << treeMPI.mpirank << " NBODY  : " << Ci->NBODY << " " << Cj->NBODY << std::endl;
	break;
      }
      real_t sumi = 0, sumj = 0;
      if (Ci->NCHILD == 0) {
	for (int i=0; i<Ci->NBODY; i++) {
	  B_iter Bi = Ci->BODY+i;
	  B_iter Bj = Cj->BODY+i;
	  sumi += Bi->X[0];
	  sumj += Bj->X[0];
	}
      }
      if (fabs(sumi-sumj)/fabs(sumi) > 1e-6) {
	std::cout << treeMPI.mpirank << " " << Ci->ICELL << " " << sumi << " " << sumj << std::endl;
      }
      assert( fabs(sumi-sumj)/fabs(sumi) < 1e-6 );
      for (int i=0; i<Ci->NCHILD; i++) Qi.push(icells.begin()+Ci->ICHILD+i);
      for (int i=0; i<Cj->NCHILD; i++) Qj.push(jcells.begin()+Cj->ICHILD+i);
      ic++;
    }
    assert( ic == int(icells.size()) );
#endif
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
#else
  for (int irank=0; irank<treeMPI.mpisize; irank++) {
    treeMPI.shiftBodies(jbodies);
    jcells.clear();
    localBounds = boundbox.getBounds(jbodies);
    jcells = build.buildTree(jbodies, localBounds);
    pass.upwardPass(jcells);
    traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
  }
#endif
  pass.downwardPass(cells);

#if 0
  treeMPI.unpartition(bodies);
  bodies = sort.sortBodies(bodies);
  bodies = treeMPI.commBodies(bodies);
  bodies = sort.unsort(bodies);
#endif
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
  boundbox.resetTimer();
  build.resetTimer();
  pass.resetTimer();
  traversal.resetTimer();
  treeMPI.resetTimer();
  logger.resetTimer();
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
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->ICELL = 0;
  for (int irank=0; irank<treeMPI.mpisize; irank++) {
    treeMPI.gettreeMPI(jcells,(treeMPI.mpirank+irank)%treeMPI.mpisize);
    for (C_iter C=jcells.begin(); C!=jcells.end(); C++) {
      Body body;
      body.ICELL = 1;
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
