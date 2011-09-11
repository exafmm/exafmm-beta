#include "vortex.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  IMAGES = 3;
  THETA = 1/sqrtf(4);
  const int numGrid1D = 128;
  const int numSteps  = 10;
  const int numSkip   = 1;
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
  T.readData(bodies,cells);
  T.stopTimer("Read data    ",printNow);
  T.eraseTimer("Read data    ");

  T.startTimer("Validate data");
  bodies2 = bodies;
  T.gridVelocity(bodies,cells);
  T.initialError(bodies);
  bodies = bodies2;
  T.stopTimer("Validate data",printNow);
  T.eraseTimer("Validate data");

  for( int step=0; step!=numSteps; ++step ) {
    T.startTimer("Statistics   ");
    if( step%(numSkip+1) == 0 ) {
      bodies2 = bodies;
      T.gridVelocity(bodies,cells);
      T.statistics(bodies,nu,dt);
    }
    T.stopTimer("Statistics   ",printNow);
    T.eraseTimer("Statistics   ");

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

    T.startTimer("Reinitialize ");
    if( step%(numSkip+1) == numSkip ) T.reinitialize(bodies,bodies2);
    T.stopTimer("Reinitialize ",printNow);
    T.eraseTimer("Reinitialize ");
    T.writeTime();
    T.resetTimer();
  }

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  T.finalize();
}
