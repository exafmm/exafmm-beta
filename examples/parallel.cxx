#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;
  Sort sort;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta);
  Traversal traversal(args.nspawn,args.images);
  LocalEssentialTree LET(args.images);
  args.numBodies /= LET.MPISIZE;
  logger.verbose = LET.MPIRANK == 0;
  args.verbose &= logger.verbose;
  if (args.verbose) {
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    LET.verbose = true;
    logger.printTitle("Parameters");
  }
  if(LET.MPIRANK == 0) args.print(logger.stringLength,P);
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  if (args.verbose) logger.printTitle("Profiling");
  logger.startTimer("Total FMM");
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, args.chargeSign, LET.MPIRANK, LET.MPISIZE);
  Bounds localBounds = boundbox.getBounds(bodies);
#if IneJ
  Bodies jbodies = data.initBodies(args.numBodies, args.distribution, LET.MPIRANK+LET.MPISIZE, LET.MPISIZE);
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
  Box box = boundbox.bounds2box(localBounds);
  logger.startPAPI();
  Cells cells = tree.buildTree(bodies, box);
  pass.upwardPass(cells);
#if IneJ
  Cells jcells = tree.buildTree(jbodies, box);
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
  for (int irank=1; irank<LET.MPISIZE; irank++) {
    LET.getLET(jcells,(LET.MPIRANK+irank)%LET.MPISIZE);

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
	std::cout << LET.MPIRANK << " ICELL  : " << Ci->ICELL << " " << Cj->ICELL << std::endl;
	break;
      }
      if (Ci->NCHILD != Cj->NCHILD) {
	std::cout << LET.MPIRANK << " NCHILD : " << Ci->NCHILD << " " << Cj->NCHILD << std::endl;
	break;
      }
      if (Ci->NCBODY != Cj->NCBODY) {
	std::cout << LET.MPIRANK << " NCBODY : " << Ci->NCBODY << " " << Cj->NCBODY << std::endl;
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
      if (fabs(sumi-sumj)/fabs(sumi) > 1e-6) std::cout << LET.MPIRANK << " " << Ci->ICELL << " " << sumi << " " << sumj << std::endl;
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
  for (int irank=0; irank<LET.MPISIZE; irank++) {
    LET.shiftBodies(jbodies);
    jcells.clear();
    localBounds = boundbox.getBounds(jbodies);
    box = boundbox.bounds2box(localBounds);
    tree.buildTree(jbodies, jcells, box);
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
  if (args.verbose) logger.printTitle("MPI direct sum");
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies2);
  logger.startTimer("Total Direct");
  for (int i=0; i<LET.MPISIZE; i++) {
    LET.shiftBodies(jbodies);
    traversal.direct(bodies2, jbodies, cycle);
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << LET.MPISIZE << std::endl;
  }
  traversal.normalize(bodies2);
  if (args.verbose) logger.printTitle("Total runtime");
  if (logger.verbose) logger.printTime("Total FMM");
  logger.stopTimer("Total Direct",logger.verbose);
  boundbox.writeTime();
  tree.writeTime();
  pass.writeTime();
  traversal.writeTime();
  LET.writeTime();
  boundbox.resetTimer();
  tree.resetTimer();
  pass.resetTimer();
  traversal.resetTimer();
  LET.resetTimer();
  logger.resetTimer();
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  data.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
  MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(args.verbose) {
    logger.printTitle("FMM vs. direct");
    logger.printError(diff3, norm3, diff4, norm4);
    tree.printTreeData(cells);
    traversal.printTraversalData();
    logger.printPAPI();
  }

#if VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) B->ICELL = 0;
  for (int irank=0; irank<LET.MPISIZE; irank++) {
    LET.getLET(jcells,(LET.MPIRANK+irank)%LET.MPISIZE);
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
  for (int i=1; i<LET.MPISIZE; i++) {
    LET.shiftBodies(jbodies);
    vtk.setGroupOfPoints(jbodies);
  }
  if (LET.MPIRANK == 0) {
    vtk.plot();
  }
#endif
  return 0;
}
