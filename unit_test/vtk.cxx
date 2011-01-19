#include "dataset.h"
#include "tree.h"
#include "vtk.h"

int main()
{
  int const Nbody=5000;
  Bodies bodies(Nbody);
  TreeStructure T(bodies);
  Dataset D(bodies);

  D.random();
  T.setDomain();

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
