#include "args.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"
#include "parallelfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args ARGS(argc, argv);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset DAT;
  Logger LOG;
  BuildTree BLD(ARGS.NCRIT,ARGS.NSPAWN);
  UpDownPass UDP(ARGS.IMAGES,ARGS.THETA);
  Traversal TRV(ARGS.NSPAWN,ARGS.IMAGES);
  ParallelFMM FMM(ARGS.NSPAWN,ARGS.IMAGES);
  LOG.printNow = FMM.MPIRANK == 0;
  BLD.printNow = FMM.MPIRANK == 0;
  UDP.printNow = FMM.MPIRANK == 0;
  TRV.printNow = FMM.MPIRANK == 0;
  FMM.printNow = FMM.MPIRANK == 0;
#if AUTO
  TRV.timeKernels();
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
    int numBodies = ARGS.numBodies / FMM.MPISIZE;
#endif
    if(FMM.printNow) std::cout << std::endl
      << "Num bodies           : " << numBodies << std::endl;
    if(FMM.printNow) std::cout << "--- Profiling --------------------" << std::endl;
    bodies.resize(numBodies);
    DAT.initBodies(bodies, ARGS.distribution, FMM.MPIRANK, FMM.MPISIZE);
    LOG.startTimer("Total FMM");

    FMM.partition(bodies);
    FMM.setLocal(bodies);
    BLD.buildTree(bodies, cells, FMM.localBox);
    UDP.upwardPass(cells);
    LOG.startPAPI();

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree : Step 1
    FMM.setLET(cells);
    FMM.commBodies();
    FMM.commCells();
    TRV.dualTreeTraversal(cells, cells, FMM.CYCLE, ARGS.mutual);
    jbodies = bodies;
    for( int irank=1; irank<FMM.MPISIZE; irank++ ) {
      FMM.getLET(jcells,(FMM.MPIRANK+irank)%FMM.MPISIZE);

#if 0 // Set to 1 for debugging full LET communication : Step 2 (LET must be set to full tree)
      FMM.shiftBodies(jbodies); // This will overwrite recvBodies. (define recvBodies2 in partition.h to avoid this)
      Cells icells;
      BLD.buildTree(jbodies, icells);
      UDP.upwardPass(icells);
      assert( icells.size() == jcells.size() );
      CellQueue Qi, Qj;
      Qi.push(icells.begin());
      Qj.push(jcells.begin());
      int ic=0;
      while( !Qi.empty() ) {
        C_iter Ci=Qi.front(); Qi.pop();
        C_iter Cj=Qj.front(); Qj.pop();
        if( Ci->ICELL != Cj->ICELL ) {
          std::cout << FMM.MPIRANK << " ICELL  : " << Ci->ICELL << " " << Cj->ICELL << std::endl;
          break;
        }
        if( Ci->NCHILD != Cj->NCHILD ) {
          std::cout << FMM.MPIRANK << " NCHILD : " << Ci->NCHILD << " " << Cj->NCHILD << std::endl;
          break;
        }
        if( Ci->NCBODY != Cj->NCBODY ) {
          std::cout << FMM.MPIRANK << " NCBODY : " << Ci->NCBODY << " " << Cj->NCBODY << std::endl;
          break;
        }
        real_t sumi = 0, sumj = 0;
        if( Ci->NCBODY != 0 ) {
          for( int i=0; i<Ci->NCBODY; i++ ) {
            B_iter Bi = Ci->BODY+i;
            B_iter Bj = Cj->BODY+i;
            sumi += Bi->X[0];
            sumj += Bj->X[0];
          }
        }
        if( fabs(sumi-sumj)/fabs(sumi) > 1e-6 ) std::cout << FMM.MPIRANK << " " << Ci->ICELL << " " << sumi << " " << sumj << std::endl;
        assert( fabs(sumi-sumj)/fabs(sumi) < 1e-6 );
        for( int i=0; i<Ci->NCHILD; i++ ) Qi.push(icells.begin()+Ci->CHILD+i);
        for( int i=0; i<Cj->NCHILD; i++ ) Qj.push(jcells.begin()+Cj->CHILD+i);
        ic++;
      }
      assert( ic == int(icells.size()) );
#endif
      TRV.dualTreeTraversal(cells, jcells, FMM.CYCLE, ARGS.mutual);
    }
#else
    jbodies = bodies;
    for( int irank=0; irank!=FMM.MPISIZE; irank++ ) {
      FMM.shiftBodies(jbodies);
      jcells.clear();
      FMM.setLocal(jbodies);
      BLD.buildTree(jbodies, jcells);
      UDP.upwardPass(jcells);
      TRV.dualTreeTraversal(cells, jcells, FMM.periodicCycle, ARGS.mutual);
    }
#endif

    UDP.downwardPass(cells);

    LOG.stopPAPI();
    if(LOG.printNow) std::cout << "--- Total runtime ----------------" << std::endl;
    LOG.stopTimer("Total FMM",LOG.printNow);
    if(LOG.printNow) std::cout << "--- Round-robin MPI direct sum ---" << std::endl;
    BLD.writeTime();
    UDP.writeTime();
    TRV.writeTime();
    FMM.writeTime();
    BLD.resetTimer();
    UDP.resetTimer();
    TRV.resetTimer();
    FMM.resetTimer();
    jbodies = bodies;
    if (int(bodies.size()) > ARGS.numTarget) DAT.sampleBodies(bodies, ARGS.numTarget);
    Bodies bodies2 = bodies;
    DAT.initTarget(bodies2);
    LOG.startTimer("Total Direct");
    for( int i=0; i!=FMM.MPISIZE; ++i ) {
      FMM.shiftBodies(jbodies);
      UDP.direct(bodies2, jbodies, FMM.CYCLE);
      if(LOG.printNow) std::cout << "Direct loop          : " << i+1 << "/" << FMM.MPISIZE << std::endl;
    }
    UDP.normalize(bodies2);
    if(LOG.printNow) std::cout << "----------------------------------" << std::endl;
    LOG.stopTimer("Total Direct",LOG.printNow);
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
    DAT.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(LOG.printNow) {
      DAT.printError(diff3, norm3, diff4, norm4);
      BLD.printTreeData(cells);
      TRV.printTraversalData();
    }
  }

#ifdef VTK
  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) B->ICELL = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }
  int Ncell = 0;
  vtkPlot vtk;
  if( FMM.MPIRANK == 0 ) {
    vtk.setDomain(M_PI,0);
    vtk.setGroupOfPoints(jbodies,Ncell);
  }
  for( int i=1; i!=FMM.MPISIZE; ++i ) {
    FMM.shiftBodies(jbodies);
    if( FMM.MPIRANK == 0 ) {
      vtk.setGroupOfPoints(jbodies,Ncell);
    }
  }
  if( FMM.MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
