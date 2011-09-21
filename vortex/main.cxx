#include "vortex.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  IMAGES = 3;
  THETA = 1/sqrtf(4);
  const int numGrid1D = 128;
  const int numSteps  = 100;
  const int numSkip   = 9;
  const float dt      = 5e-3;
  const float nu      = 1e-2;
  Bodies bodies, bodies2;
  Cells cells,jcells;
  Vortex T(numGrid1D);
  T.setKernel("Gaussian");
  T.initialize();
  bool printNow = (MPIRANK == 0);

  T.startTimer("Read data    ");
  bodies.resize(T.numBodies);
  T.readData(bodies,cells,jcells);
  T.stopTimer("Read data    ",printNow);
  T.eraseTimer("Read data    ");

  T.startTimer("Validate data");
  bodies2 = bodies;
  T.gridVelocity(bodies,cells,jcells);
  T.initialError(bodies);
  bodies = bodies2;
  T.stopTimer("Validate data",printNow);
  T.eraseTimer("Validate data");

  for( int step=0; step!=numSteps; ++step ) {
    if( step%(numSkip+1) == 0 ) {
      T.startTimer("Statistics   ");
      bodies2 = bodies;
      T.gridVelocity(bodies,cells,jcells);
      T.statistics(bodies,nu,dt);
      bodies = bodies2;
      T.stopTimer("Statistics   ",printNow);
      T.eraseTimer("Statistics   ");
    }

    T.startTimer("BiotSavart   ");
    T.BiotSavart(bodies,cells,jcells);
    T.stopTimer("BiotSavart   ",printNow);
    T.eraseTimer("BiotSavart   ");

    T.startTimer("Stretching   ");
    T.Stretching(bodies,cells,jcells);
    T.stopTimer("Stretching   ",printNow);
    T.eraseTimer("Stretching   ");

    T.startTimer("Convect      ");
    T.update(bodies,nu,dt);
    T.stopTimer("Convect      ",printNow);
    T.eraseTimer("Convect      ");

    if( step%(numSkip+1) == numSkip ) {
      T.startTimer("Reinitialize ");
      T.reinitialize(bodies,bodies2,cells,jcells);
      T.stopTimer("Reinitialize ",printNow);
      T.eraseTimer("Reinitialize ");
    }
    T.writeTime();
    T.resetTimer();
  }

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->ICELL = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
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
