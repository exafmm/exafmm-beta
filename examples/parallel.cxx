#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset data;
  Logger logger;
  Sort sort;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.NSPAWN);
  BuildTree tree(args.NCRIT,args.NSPAWN);
  UpDownPass pass(args.IMAGES,args.THETA);
  Traversal traversal(args.NSPAWN,args.IMAGES);
  LocalEssentialTree LET(args.IMAGES);
  logger.printNow = LET.MPIRANK == 0;
  boundbox.printNow = LET.MPIRANK == 0;
  tree.printNow = LET.MPIRANK == 0;
  pass.printNow = LET.MPIRANK == 0;
  traversal.printNow = LET.MPIRANK == 0;
  LET.printNow = LET.MPIRANK == 0;
#if AUTO
  traversal.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
    int numBodies = int(pow(10,(it+24)/8.0));
#else
  {
    int numBodies = args.numBodies / LET.MPISIZE;
#endif
    if (LET.printNow) std::cout << std::endl
      << "Num bodies           : " << numBodies << std::endl;
    if (LET.printNow) std::cout << "--- Profiling --------------------" << std::endl;
    bodies.resize(numBodies);
    data.initBodies(bodies, args.distribution, LET.MPIRANK, LET.MPISIZE);
    logger.startTimer("Total FMM");
    Bounds localBounds = boundbox.getBounds(bodies);
    Bounds globalBounds = LET.allreduceBounds(localBounds);
    localBounds = LET.partition(bodies,globalBounds);
    bodies = sort.sortBodies(bodies);
    bodies = LET.commBodies(bodies);
    Box box = boundbox.bounds2box(localBounds);
    tree.buildTree(bodies, cells, box);
    pass.upwardPass(cells);
    logger.startPAPI();

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree : Step 1
    LET.setLET(cells,localBounds,cycle);
    LET.commBodies();
    LET.commCells();
    traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
    jbodies = bodies;
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
      traversal.dualTreeTraversal(cells, jcells, cycle, args.mutual);
    }
#else
    jbodies = bodies;
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
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->ICELL = B->IBODY;                                      //  Do this to sort accroding to IPROC
    }                                                           // End loop over bodies
    Bodies buffer = bodies;                                     // Resize sort buffer
    bodies = sort.sortBodies(buffer);                                // Sort bodies in ascending order
#endif
    logger.stopPAPI();
    if (logger.printNow) std::cout << "--- Total runtime ----------------" << std::endl;
    logger.stopTimer("Total FMM",logger.printNow);
    if (logger.printNow) std::cout << "--- Round-robin MPI direct sum ---" << std::endl;
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
    jbodies = bodies;
    if (int(bodies.size()) > args.numTarget) data.sampleBodies(bodies, args.numTarget);
    Bodies bodies2 = bodies;
    data.initTarget(bodies2);
    logger.startTimer("Total Direct");
    for (int i=0; i<LET.MPISIZE; i++) {
      LET.shiftBodies(jbodies);
      pass.direct(bodies2, jbodies, cycle);
      if (logger.printNow) std::cout << "Direct loop          : " << i+1 << "/" << LET.MPISIZE << std::endl;
    }
    pass.normalize(bodies2);
    if (logger.printNow) std::cout << "----------------------------------" << std::endl;
    logger.stopTimer("Total Direct",logger.printNow);
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
    data.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(logger.printNow) {
      data.printError(diff3, norm3, diff4, norm4);
      tree.printTreeData(cells);
      traversal.printTraversalData();
    }
  }

#ifdef VTK
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); ++B) B->ICELL = 0;
  for (C_iter C=jcells.begin(); C!=jcells.end(); ++C) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }
  int Ncell = 0;
  vtkPlot vtk;
  if (LET.MPIRANK == 0) {
    vtk.setDomain(M_PI,0);
    vtk.setGroupOfPoints(jbodies,Ncell);
  }
  for (int i=1; i<LET.MPISIZE; i++) {
    LET.shiftBodies(jbodies);
    if (LET.MPIRANK == 0) {
      vtk.setGroupOfPoints(jbodies,Ncell);
    }
  }
  if (LET.MPIRANK == 0) {
    vtk.plot(Ncell);
  }
#endif
  return 0;
}
