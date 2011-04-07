#include "vortex.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numGrid1D = 128;
  const int numSteps  = 1;
  const int numSkip   = 0;
  const float dt      = 5e-4;
  const float nu      = 2e-2;
  Bodies bodies,jbodies;
  Cells cells,jcells;
  Vortex T(numGrid1D);
  T.setKernel("Gaussian");
  T.initialize();

  T.startTimer("Read data    ");
  bodies.resize(T.numBodies);
  T.readData(bodies);
  T.stopTimer("Read data    ",true);
  T.eraseTimer("Read data    ");

  T.startTimer("Validate data");
  T.initialError(bodies);
  T.stopTimer("Validate data",true);
  T.eraseTimer("Validate data");

  for( int step=0; step!=numSteps; ++step ) {
    T.startTimer("Statistics   ");
    T.statistics(bodies,step%(numSkip+1)==0);
    T.stopTimer("Statistics   ",true);
    T.eraseTimer("Statistics   ");

    T.startTimer("Stretching   ");
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    T.setKernel("Stretching");
    cells.clear();
    T.setDomain(bodies);
    T.bottomup(bodies,cells);
    jcells = cells;
    T.downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    T.stopTimer("Stretching   ",true);
    T.eraseTimer("Stretching   ");

    T.startTimer("Convect      ");
    T.convect(bodies,nu,dt);
    T.stopTimer("Convect      ",true);
    T.eraseTimer("Convect      ");

    T.startTimer("Reinitialize ");
    if(step%(numSkip+1)==numSkip) T.reinitialize(bodies);
    T.stopTimer("Reinitialize ",true);
    T.eraseTimer("Reinitialize ");

    T.startTimer("BiotSavart   ");
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    T.setKernel("BiotSavart");
    cells.clear();
    T.setDomain(bodies);
    T.bottomup(bodies,cells);
    jcells = cells;
    T.downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    T.stopTimer("BiotSavart   ",true);
    T.eraseTimer("BiotSavart   ");
  }
  T.startTimer("Statistics   ");
  T.statistics(bodies);
  T.stopTimer("Statistics   ",true);
  T.eraseTimer("Statistics   ");

#ifdef VTK
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(bodies,Ncell);
  vtk.plot(Ncell);
#endif
  T.finalize();
}
