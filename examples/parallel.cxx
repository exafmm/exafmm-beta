#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  BaseMPI baseMPI;
  BoundBox boundBox(args.nspawn);
  BuildTree buildTree(args.ncrit, args.nspawn);
  Dataset data;
  Logger logger;
  Partition partition;
  Traversal traversal(args.nspawn, args.images);
  TreeMPI treeMPI(args.images);
  UpDownPass upDownPass(args.theta);
  Verify verify;

  const real_t cycle = 2 * M_PI;

  args.numBodies /= baseMPI.mpisize;
  args.verbose &= baseMPI.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundBox.verbose = true;
    buildTree.verbose = true;
    upDownPass.verbose = true;
    traversal.verbose = true;
    treeMPI.verbose = true;
    verify.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength, P, baseMPI.mpirank);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] += M_PI;
    B->X[0] *= 0.5;
  }
  Bounds localBounds = boundBox.getBounds(bodies);
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank+baseMPI.mpisize, baseMPI.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] -= M_PI;
    B->X[0] *= 0.5;
  }
  localBounds = boundBox.getBounds(jbodies,localBounds);
#endif
  Bounds globalBounds = baseMPI.allreduceBounds(localBounds);
  localBounds = partition.octsection(bodies,globalBounds);
  bodies = treeMPI.commBodies(bodies);
#if IneJ
  partition.octsection(jbodies,globalBounds);
  jbodies = treeMPI.commBodies(jbodies);
#endif
  localBounds = boundBox.getBounds(bodies);
  Cells cells = buildTree.buildTree(bodies, localBounds);
  upDownPass.upwardPass(cells);
#if IneJ
  localBounds = boundBox.getBounds(jbodies);
  Cells jcells = buildTree.buildTree(jbodies, localBounds);
  upDownPass.upwardPass(jcells);
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
  for (int irank=1; irank<baseMPI.mpisize; irank++) {
    treeMPI.getLET(jcells,(baseMPI.mpirank+irank)%baseMPI.mpisize);

#if 0 // Set to 1 for debugging full treeMPI communication : Step 2 (treeMPI must be set to full tree)
    treeMPI.shiftBodies(jbodies); // This will overwrite recvBodies. (define recvBodies2 in body_mpi.h to avoid this)
    Cells icells;
    buildTree.buildTree(jbodies, icells);
    upDownPass.upwardPass(icells);
    assert( icells.size() == jcells.size() );
    CellQueue Qi, Qj;
    Qi.push(icells.begin());
    Qj.push(jcells.begin());
    int ic=0;
    while (!Qi.empty()) {
      C_iter Ci=Qi.front(); Qi.pop();
      C_iter Cj=Qj.front(); Qj.pop();
      if (Ci->ICELL != Cj->ICELL) {
	std::cout << baseMPI.mpirank << " ICELL  : " << Ci->ICELL << " " << Cj->ICELL << std::endl;
	break;
      }
      if (Ci->NCHILD != Cj->NCHILD) {
	std::cout << baseMPI.mpirank << " NCHILD : " << Ci->NCHILD << " " << Cj->NCHILD << std::endl;
	break;
      }
      if (Ci->NCHILD == 0 && Cj->NCHILD == 0 && Ci->NBODY != Cj->NBODY) {
	std::cout << baseMPI.mpirank << " NBODY  : " << Ci->NBODY << " " << Cj->NBODY << std::endl;
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
	std::cout << baseMPI.mpirank << " " << Ci->ICELL << " " << sumi << " " << sumj << std::endl;
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
  for (int irank=0; irank<baseMPI.mpisize; irank++) {
    treeMPI.shiftBodies(jbodies);
    jcells.clear();
    localBounds = boundBox.getBounds(jbodies);
    jcells = buildTree.buildTree(jbodies, localBounds);
    upDownPass.upwardPass(jcells);
    traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
  }
#endif
  upDownPass.downwardPass(cells);

  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  logger.printTitle("MPI direct sum");
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  for (int i=0; i<baseMPI.mpisize; i++) {
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
    treeMPI.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
  }
  traversal.normalize(bodies);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
#if WRITE_TIME
  boundBox.writeTime(baseMPI.mpirank);
  buildTree.writeTime(baseMPI.mpirank);
  upDownPass.writeTime(baseMPI.mpirank);
  traversal.writeTime(baseMPI.mpirank);
  treeMPI.writeTime(baseMPI.mpirank);
#endif
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
  buildTree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();

#if VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->IBODY = 0;
  for (int irank=0; irank<baseMPI.mpisize; irank++) {
    treeMPI.gettreeMPI(jcells,(baseMPI.mpirank+irank)%baseMPI.mpisize);
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
  for (int i=1; i<baseMPI.mpisize; i++) {
    treeMPI.shiftBodies(jbodies);
    vtk.setGroupOfPoints(jbodies);
  }
  if (baseMPI.mpirank == 0) {
    vtk.plot();
  }
#endif
  return 0;
}
