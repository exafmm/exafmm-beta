#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 1000000;
  std::string kernelName = "BiotSavart";
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
  T.resetTimer();

  T.setKernel("BiotSavart");
  D.initTarget(bodies);
  T.evalP2M(cells);
  T.evalM2M(cells);
  jcells = cells;
  if( MPISIZE != 1 ) {
    #pragma omp parallel sections num_threads(2)
    {
      #pragma omp section
      {
         T.downward(cells,jcells,1,false);
      }
      #pragma omp section
      {
        T.updateBodies();
      }
    }
    jbodies = bodies;
    jcells = cells;
    T.commCells(jbodies,jcells);
    T.eraseLocalTree(jcells);
  }
  T.downward(cells,jcells,1);
  T.setKernel("Stretching");
  D.initTarget(bodies);
  T.downward(cells,jcells,1);
  if(T.printNow) T.writeTime();
  if(T.printNow) T.writeTime();
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
