#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 1000;
  std::string kernelName = "Laplace";
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  if( MPIRANK == 0 ) T.printNow = true;

  for( int it=0; it!=9; ++it ) {
    numBodies = int(pow(10,(it+24)/8.0));
    if(T.printNow) std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    D.random(bodies,MPIRANK+1);
    T.startTimer("FMM          ");
    T.setGlobDomain(bodies);
    T.octsection(bodies);
    cells.clear();
#ifdef TOPDOWN
    T.topdown(bodies,cells);
#else
    T.bottomup(bodies,cells);
#endif
    T.commBodies(cells);
    jbodies = bodies;
    Cells jcells = cells;
    T.commCells(jbodies,jcells);

    T.downward(cells,jcells,1);
    T.stopTimer("FMM          ",T.printNow);
    T.eraseTimer("FMM          ");

    T.startTimer("Direct sum   ");
    Bodies bodies2 = bodies;
#if 1
    D.initTarget(bodies2);
    if( IMAGES != 0 ) {
      jbodies = T.periodicBodies(bodies2);
    } else {
      jbodies = bodies2;
    }
    for( int i=0; i!=MPISIZE; ++i ) {
      T.shiftBodies(jbodies);
      T.evalP2P(bodies2,jbodies);
      if(T.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
    }
    D.writeTarget(bodies2);
#else
    D.readTarget(bodies2);
#endif
    T.stopTimer("Direct sum   ",T.printNow);
    T.eraseTimer("Direct sum   ");
    if(T.printNow) T.writeTime();
    if(T.printNow) T.resetTimer();

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
  }

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
