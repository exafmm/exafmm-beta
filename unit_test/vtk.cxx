#include "dataset.h"
#include "tree.h"
#include "vtk.h"

int main() {
  const int numBodies = 10000;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset dataset;
  dataset.kernelName = "Laplace";
  TreeStructure FMM;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.sphere(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  vtkPlot vtk;
  vtk.setDomain(FMM.getR0(),FMM.getX0());
  vtk.setGroup(0,bodies.size()/2);
  for( B_iter B=bodies.begin(); B!=bodies.begin()+bodies.size()/2; ++B ) {
    vtk.setPoints(0,B->X);
  }
  vtk.setGroup(1,bodies.size()/2);
  for( B_iter B=bodies.begin()+bodies.size()/2; B!=bodies.end(); ++B ) {
    vtk.setPoints(1,B->X);
  }
  vtk.plot(2);
  FMM.finalize();
}
