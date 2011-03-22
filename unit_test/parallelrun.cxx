#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  Evaluator E;
  if( T.commRank() == 0 ) T.printNow = true;

  T.startTimer("Set bodies   ");
  D.sphere(bodies,T.commRank()+1);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setGlobDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

  T.bisection(bodies);

#ifdef TOPDOWN
  T.topdown(bodies,cells);
#else
  T.bottomup(bodies,cells);
#endif

  T.commBodies(cells);

  Bodies bodies2 = bodies;
  Cells jcells = cells;
  T.commCells(bodies2,jcells);

  T.downward(cells,jcells,1);

  if( IMAGES != 0 ) {
    T.startTimer("Set periodic ");
    jbodies = T.periodicBodies(bodies);
    T.stopTimer("Set periodic ",T.printNow);
  } else {
    jbodies = bodies;
  }

  T.startTimer("Direct sum   ");
  bodies2 = bodies;
  D.initTarget(bodies2);
  for( int i=0; i!=T.commSize(); ++i ) {
    T.shiftBodies(jbodies);
    E.evalP2P(bodies2,jbodies);
    if(T.printNow) std::cout << "Direct loop   : " << i+1 << "/" << T.commSize() << std::endl;
  }
  T.stopTimer("Direct sum   ",T.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = T.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(T.printNow) D.printError(diff1,norm1,diff2,norm2);
#ifdef DEBUG
  T.print(std::sqrt(potDiff/potNorm));
#endif

#ifdef VTK
  for( B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.I = 1;
    body.X = C->X;
    body.Q = 0;
    bodies.push_back(body);
  }

  int Ncell = 0;
  vtkPlot vtk;
  if( T.commRank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=T.commSize(); ++i ) {
    T.shiftBodies(bodies);
    if( T.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( T.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
