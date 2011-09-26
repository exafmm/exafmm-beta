#include "vortex.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  IMAGES = 3;
  THETA = 1/sqrtf(4);
  const int numGrid1D = 512;
  const int numSteps  = 1;
  const int numSkip   = 0;
  const int numSkip2  = 0;
  const float dt      = 5e-3;
  const float nu      = 1e-2;
  Bodies bodies, bodies2;
  Cells cells;
  Vortex T(numGrid1D);
  T.setKernel("Gaussian");
  T.initialize();
  bool printNow = (MPIRANK == 0);

  T.startTimer("Read data    ");
  bodies.resize(T.numBodies);
  bodies2.resize(T.numBodies);
  T.readData(bodies,bodies2,cells);
  T.stopTimer("Read data    ",printNow);
  T.eraseTimer("Read data    ");

  T.startTimer("Validate data");
  T.gridVelocity(bodies,bodies2,cells);
  T.initialError(bodies2);
  T.stopTimer("Validate data",printNow);
  T.eraseTimer("Validate data");

  for( int step=0; step!=numSteps; ++step ) {
    T.print("Step         : ",0);
    T.print(step,0);
    T.print("\n",0);
    if( step%(numSkip+1) == 0 ) {
      T.startTimer("Statistics   ");
      T.gridVelocity(bodies,bodies2,cells);
      T.statistics(bodies2,nu,dt);
      T.stopTimer("Statistics   ",printNow);
      T.eraseTimer("Statistics   ");
    }

    T.startTimer("BiotSavart   ");
    T.BiotSavart(bodies,cells);
    T.stopTimer("BiotSavart   ",printNow);
    T.eraseTimer("BiotSavart   ");

    T.startTimer("Stretching   ");
    T.Stretching(bodies,cells);
    T.stopTimer("Stretching   ",printNow);
    T.eraseTimer("Stretching   ");

    T.startTimer("Convect      ");
    T.update(bodies,nu,dt);
    T.stopTimer("Convect      ",printNow);
    T.eraseTimer("Convect      ");

    if( step%(numSkip2+1) == numSkip2 ) {
      T.startTimer("Reinitialize ");
      T.reinitialize(bodies,cells);
      T.stopTimer("Reinitialize ",printNow);
      T.eraseTimer("Reinitialize ");
    }
    T.writeTime();
    T.resetTimer();
  }

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->ICELL = 0;
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
    Body body;
    body.ICELL = 1;
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
