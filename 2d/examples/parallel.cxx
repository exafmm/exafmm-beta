#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;
  Sort sort;

  const real_t eps2 = 0.0;
  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta,eps2);
  Traversal traversal(args.nspawn,args.images,eps2);
  LocalEssentialTree LET(args.images);
  args.numBodies /= LET.mpisize;
  args.verbose &= LET.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    LET.verbose = true;
  }
  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength, P);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies/1000, args.distribution, LET.mpirank, LET.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] += M_PI;
    B->X[0] *= 0.5;
  }
  Bounds localBounds = boundbox.getBounds(bodies);
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, LET.mpirank+LET.mpisize, LET.mpisize);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] -= M_PI;
    B->X[0] *= 0.5;
  }
  localBounds = boundbox.getBounds(jbodies,localBounds);
#endif
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
#if IneJ
  LET.partition(jbodies,globalBounds);
  jbodies = sort.sortBodies(jbodies);
  jbodies = LET.commBodies(jbodies);
#endif
  Cells cells = tree.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
#if IneJ
  Cells jcells = tree.buildTree(jbodies, localBounds);
  pass.upwardPass(jcells);
#endif

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree : Step 1
#if IneJ
  LET.setLET(jcells,localBounds,cycle);
#else
  LET.setLET(cells,localBounds,cycle);
#endif
  LET.commBodies();
  LET.commCells();
#if IneJ
  traversal.dualTreeTraversal(cells, jcells, cycle);
#else
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Bodies jbodies = bodies;
  Cells jcells;
#endif
  for (int irank=1; irank<LET.mpisize; irank++) {
    LET.getLET(jcells,(LET.mpirank+irank)%LET.mpisize);

#if 0 // Set to 1 for debugging full LET communication : Step 2 (LET must be set to full tree)
    LET.shiftBodies(jbodies); // This will overwrite recvBodies. (define recvBodies2 in partition.h to avoid this)
    Cells icells;
    tree.buildTree(jbodies, icells);
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
	std::cout << LET.mpirank << " ICELL  : " << Ci->ICELL << " " << Cj->ICELL << std::endl;
	break;
      }
      if (Ci->NCHILD != Cj->NCHILD) {
	std::cout << LET.mpirank << " NCHILD : " << Ci->NCHILD << " " << Cj->NCHILD << std::endl;
	break;
      }
      if (Ci->NCBODY != Cj->NCBODY) {
	std::cout << LET.mpirank << " NCBODY : " << Ci->NCBODY << " " << Cj->NCBODY << std::endl;
	break;
      }
      real_t sumi = 0, sumj = 0;
      if (Ci->NCBODY != 0) {
	for (int i=0; i<Ci->NCBODY; i++) {
	  B_iter Bi = Ci->BODY+i;
	  B_iter Bj = Cj->BODY+i;
	  sumi += Bi->X[0];
	  sumj += Bj->X[0];
	}
      }
      if (fabs(sumi-sumj)/fabs(sumi) > 1e-6) std::cout << LET.mpirank << " " << Ci->ICELL << " " << sumi << " " << sumj << std::endl;
      assert( fabs(sumi-sumj)/fabs(sumi) < 1e-6 );
      for (int i=0; i<Ci->NCHILD; i++) Qi.push(icells.begin()+Ci->CHILD+i);
      for (int i=0; i<Cj->NCHILD; i++) Qj.push(jcells.begin()+Cj->CHILD+i);
      ic++;
    }
    assert( ic == int(icells.size()) );
#endif
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
#else
  for (int irank=0; irank<LET.mpisize; irank++) {
    LET.shiftBodies(jbodies);
    jcells.clear();
    localBounds = boundbox.getBounds(jbodies);
    jcells = tree.buildTree(jbodies, localBounds);
    pass.upwardPass(jcells);
    traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
  }
#endif
  pass.downwardPass(cells);

#if 0
  LET.unpartition(bodies);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->ICELL = B->IBODY;
  }
  bodies = sort.sortBodies(bodies);
#endif
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  logger.printTitle("MPI direct sum");
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  for (int i=0; i<LET.mpisize; i++) {
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << LET.mpisize << std::endl;
    LET.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
  }
  traversal.normalize(bodies);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
  boundbox.writeTime(LET.mpirank);
  tree.writeTime(LET.mpirank);
  pass.writeTime(LET.mpirank);
  traversal.writeTime(LET.mpirank);
  LET.writeTime(LET.mpirank);
  LET.writeSendCount();
  boundbox.resetTimer();
  tree.resetTimer();
  pass.resetTimer();
  traversal.resetTimer();
  LET.resetTimer();
  logger.resetTimer();
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  data.evalError(bodies2, bodies, diff1, norm1);
  MPI_Reduce(&diff1, &diff2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  logger.printTitle("FMM vs. direct");
  logger.printError(diff2, norm2);
  tree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();
  MPI_Finalize();
  return 0;
}
