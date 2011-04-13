#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 100000;
  std::string kernelName = "Laplace";
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  if( MPIRANK == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.random(bodies,MPIRANK+1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

  T.octsection(bodies);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  T.commBodies(cells);

  jbodies = bodies;
  Cells jcells = cells;
  T.commCells(jbodies,jcells);

  T.startTimer("Downward     ");
  T.downward(cells,jcells,1);
  T.stopTimer("Downward     ",T.printNow);
  T.eraseTimer("Downward     ");

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies = T.periodicBodies(bodies);
    T.stopTimer("Set periodic ",T.printNow);
    T.eraseTimer("Set periodic ");
  } else {
    jbodies = bodies;
  }
  T.startTimer("Direct sum   ");
  Bodies bodies2 = bodies;
  D.initTarget(bodies2);
  for( int i=0; i!=MPISIZE; ++i ) {
    T.shiftBodies(jbodies);
    T.evalP2P(bodies2,jbodies);
    if(T.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  T.stopTimer("Direct sum   ",T.printNow);
  T.resetTimer();

  D.initTarget(bodies);
  T.evalP2M(cells);
  T.evalM2M(cells);
  T.commBodies(cells);
  jbodies = bodies;
  jcells = cells;
  T.commCells(jbodies,jcells);
  T.downward(cells,jcells,1);
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
  for( B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.I = 1;
    body.X = C->X;
    body.SRC = 0;
    bodies.push_back(body);
  }

  int Ncell = 0;
  vtkPlot vtk;
  if( MPIRANK == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    T.shiftBodies(bodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  T.finalize();
}
