#include "vortex.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  IMAGES = 3;
  THETA = 1/sqrtf(4);
  const int numGrid1D = 128;
  const int numSteps  = 1;
  const int numSkip   = 0;
  const int numSkip2  = 0;
  const float dt      = 5e-3;
  const float nu      = 1e-2;
  Bodies bodies, bodies2;
  Cells cells;
  Vortex FMM(numGrid1D);
  FMM.initialize();
  bool printNow = (MPIRANK == 0);

  FMM.startTimer("Read data    ");
  bodies.resize(FMM.numBodies);
  bodies2.resize(FMM.numBodies);
  FMM.readData(bodies,bodies2,cells);
  FMM.stopTimer("Read data    ",printNow);
  FMM.eraseTimer("Read data    ");

  FMM.startTimer("Validate data");
  FMM.gridVelocity(bodies,bodies2,cells);
  FMM.initialError(bodies2);
  FMM.stopTimer("Validate data",printNow);
  FMM.eraseTimer("Validate data");

  for( int step=0; step!=numSteps; ++step ) {
    FMM.print("Step          : ",0);
    FMM.print(step,0);
    FMM.print("\n",0);
    if( step%(numSkip+1) == 0 ) {
      FMM.startTimer("Statistics   ");
      FMM.gridVelocity(bodies,bodies2,cells);
      FMM.statistics(bodies2,nu,dt);
      FMM.stopTimer("Statistics   ",printNow);
      FMM.eraseTimer("Statistics   ");
    }

    FMM.startTimer("BiotSavart   ");
    FMM.BiotSavart(bodies,cells);
    FMM.stopTimer("BiotSavart   ",printNow);
    FMM.eraseTimer("BiotSavart   ");

    FMM.startTimer("Stretching   ");
    FMM.Stretching(bodies,cells);
    FMM.stopTimer("Stretching   ",printNow);
    FMM.eraseTimer("Stretching   ");

    FMM.startTimer("Convect      ");
    FMM.update(bodies,nu,dt);
    FMM.stopTimer("Convect      ",printNow);
    FMM.eraseTimer("Convect      ");

    if( step%(numSkip2+1) == numSkip2 ) {
      FMM.startTimer("Reinitialize ");
      FMM.reinitialize(bodies,cells);
      FMM.stopTimer("Reinitialize ",printNow);
      FMM.eraseTimer("Reinitialize ");
    }
    if(printNow) FMM.writeTime();
    FMM.resetTimer();
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
    vtk.setDomain(FMM.getR0(),FMM.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    FMM.shiftBodies(bodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  FMM.finalize();
}
