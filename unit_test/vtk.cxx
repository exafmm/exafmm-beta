#include "dataset.h"
#include "tree.h"
#include "vtk.h"

int main() {
  int const numBodies(5000);
  double tic,toc;
  Bodies bodies(numBodies);
  Dataset D;
  TreeStructure T;

  tic = get_time();
  D.sphere(bodies);
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain(bodies);
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroup(0,bodies.size()/2);
  for( B_iter B=bodies.begin(); B!=bodies.begin()+bodies.size()/2; ++B ) {
    vtk.setPoints(0,B->pos);
  }
  vtk.setGroup(1,bodies.size()/2);
  for( B_iter B=bodies.begin()+bodies.size()/2; B!=bodies.end(); ++B ) {
    vtk.setPoints(1,B->pos);
  }
  vtk.plot(2);
}
