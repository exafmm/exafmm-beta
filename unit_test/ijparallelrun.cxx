#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numTarget);
  Bodies jbodies(numBodies);
  Bodies jbodies2;
  Cells cells,jcells;
  Dataset D;
  D.kernelName = "Laplace";
  LocalEssentialTree T;
  T.setKernel(D.kernelName);
  T.initialize();
  if( MPIRANK == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,MPIRANK+1);
  D.random(jbodies,MPIRANK+MPISIZE+1);
  Bodies bodies2 = bodies;
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies2);
  T.stopTimer("Set domain   ",T.printNow);

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies2 = T.periodicBodies(jbodies);
    T.stopTimer("Set periodic ",T.printNow);
    T.eraseTimer("Set periodic ");
  } else {
    jbodies2 = jbodies;
  }

  T.startTimer("Direct sum   ");
  for( int i=0; i!=MPISIZE; ++i ) {
    T.shiftBodies(jbodies2);
    T.evalP2P(bodies2,jbodies2);
    if(T.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  T.stopTimer("Direct sum   ",T.printNow);
  T.eraseTimer("Direct sum   ");

  D.initTarget(bodies);

  T.octsection(bodies);
  T.octsection(jbodies);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
  T.topdown(jbodies,jcells);
#else
  T.bottomup(bodies,cells);
  T.bottomup(jbodies,jcells);
#endif

  T.commBodies(jcells);

  T.commCells(jbodies,jcells);

  T.startTimer("Downward     ");
  T.downward(cells,jcells,1);
  T.stopTimer("Downward     ",T.printNow);
  T.eraseTimer("Downward     ");

  T.startTimer("Unpartition  ");
  T.unpartition(bodies);
  T.stopTimer("Unpartition  ",T.printNow);
  T.eraseTimer("Unpartition  ");

  T.startTimer("Unsort bodies");
  std::sort(bodies.begin(),bodies.end());
  T.stopTimer("Unsort bodies",T.printNow);
  T.eraseTimer("Unsort bodies");
  if(T.printNow) T.writeTime();
  if(T.printNow) T.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = T.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(T.printNow) D.printError(diff3,norm3,diff4,norm4);
#ifdef DEBUG
  T.print(std::sqrt(potDiff/potNorm));
#endif

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
  if( MPIRANK == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(jbodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    T.shiftBodies(jbodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(jbodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
