#include "dataset.h"
#include "tree.h"
#include "vtk.h"

int main() {
  const int numBodies = 10000;
  std::string kernelName = "Laplace";
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset D;
  TreeStructure T;
  T.setKernel(kernelName);
  T.initialize();
  D.kernelName = kernelName;
  T.printNow = true;

  T.startTimer("Set bodies   ");
  D.sphere(bodies);
  T.stopTimer("Set bodies   ",T.printNow);

  T.startTimer("Set domain   ");
  T.setDomain(bodies);
  T.stopTimer("Set domain   ",T.printNow);

  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroup(0,bodies.size()/2);
  for( B_iter B=bodies.begin(); B!=bodies.begin()+bodies.size()/2; ++B ) {
    vtk.setPoints(0,B->X);
  }
  vtk.setGroup(1,bodies.size()/2);
  for( B_iter B=bodies.begin()+bodies.size()/2; B!=bodies.end(); ++B ) {
    vtk.setPoints(1,B->X);
  }
  vtk.plot(2);
  T.finalize();
}
