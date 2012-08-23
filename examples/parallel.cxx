#include "dataset.h"
#include "parallelfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset DATA;
  ParallelFMM FMM;
  FMM.printNow = FMM.MPIRANK == 0;
#if AUTO
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
  numBodies = int(pow(10,(it+24)/8.0));
#else
  {
  numBodies = 1000000 / FMM.MPISIZE;
#endif
  if(FMM.printNow) std::cout << "N                    : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.initBodies(bodies,FMM.MPIRANK,FMM.MPISIZE);
  FMM.startTimer("FMM");

  FMM.partition(bodies);
  FMM.setBounds(bodies);
  FMM.buildTree(bodies,cells);
  FMM.upwardPass(cells);
  FMM.startPAPI();

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree : Step 1
  FMM.setLET(cells);
  FMM.commBodies();
  FMM.commCells();
  FMM.evaluate(cells,cells);
  jbodies = bodies;
  for( int irank=1; irank<FMM.MPISIZE; irank++ ) {
    FMM.getLET(jcells,(FMM.MPIRANK+irank)%FMM.MPISIZE);

#if 0 // Set to 1 for debugging full LET communication : Step 2 (LET must be set to full tree)
    FMM.shiftBodies(jbodies); // This will overwrite recvBodies. (define recvBodies2 in partition.h to avoid this)
    Cells icells;
    FMM.buildTree(jbodies,icells);
    FMM.upwardPass(icells);
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
      if( Ci->NCLEAF != Cj->NCLEAF ) {
        std::cout << FMM.MPIRANK << " NCLEAF : " << Ci->NCLEAF << " " << Cj->NCLEAF << std::endl;
        break;
      }
      real_t sumi = 0, sumj = 0;
      if( Ci->NCLEAF != 0 ) {
        for( int i=0; i<Ci->NCLEAF; i++ ) {
          B_iter Bi = Ci->LEAF+i;
          B_iter Bj = Cj->LEAF+i;
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
    FMM.evaluate(cells,jcells);
  }
#else
  jbodies = bodies;
  for( int irank=0; irank!=FMM.MPISIZE; irank++ ) {
    FMM.shiftBodies(jbodies);
    jcells.clear();
    FMM.setBounds(jbodies);
    FMM.buildTree(jbodies,jcells);
    FMM.upwardPass(jcells);
    FMM.evaluate(cells,jcells);
  }
#endif

  FMM.downwardPass(cells);

  FMM.stopPAPI();
  FMM.stopTimer("FMM",FMM.printNow);
  FMM.eraseTimer("FMM");
  FMM.writeTime();
  FMM.resetTimer();

  jbodies = bodies;
  if (bodies.size() > 100) bodies.resize(100);
  Bodies bodies2 = bodies;
  DATA.initTarget(bodies2);
  FMM.startTimer("Direct sum");
  for( int i=0; i!=FMM.MPISIZE; ++i ) {
    FMM.shiftBodies(jbodies);
    FMM.direct(bodies2,jbodies);
    if(FMM.printNow) std::cout << "Direct loop          : " << i+1 << "/" << FMM.MPISIZE << std::endl;
  }
  FMM.normalize(bodies2);
  FMM.stopTimer("Direct sum",FMM.printNow);
  FMM.eraseTimer("Direct sum");
  real_t diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = FMM.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(FMM.printNow) DATA.printError(diff3,norm3,diff4,norm4);
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
