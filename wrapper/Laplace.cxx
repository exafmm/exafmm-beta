#include "construct.h"

void laplace(int numBodies, float *x, float *y, float *z, float *charge,
             float *pot, float *fx, float *fy, float *fz) {
  Bodies bodies(numBodies);
  Cells cells,jcells;
  TreeConstructor T;
  T.printNow = true;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->X[0]   = x[B-bodies.begin()];
    B->X[1]   = y[B-bodies.begin()];
    B->X[2]   = z[B-bodies.begin()];
    B->Q      = charge[B-bodies.begin()];
    B->pot    = pot[B-bodies.begin()];
    B->acc[0] = fx[B-bodies.begin()];
    B->acc[1] = fy[B-bodies.begin()];
    B->acc[2] = fz[B-bodies.begin()];
    B->IBODY  = B-bodies.begin();
  }

  T.setDomain(bodies);
  T.bottomup(bodies,cells);
  jcells = cells;
  T.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    x[B-bodies.begin()]      = B->X[0];
    y[B-bodies.begin()]      = B->X[1];
    z[B-bodies.begin()]      = B->X[2];
    charge[B-bodies.begin()] = B->Q;
    pot[B-bodies.begin()]    = B->pot;
    fx[B-bodies.begin()]     = B->acc[0];
    fy[B-bodies.begin()]     = B->acc[1];
    fz[B-bodies.begin()]     = B->acc[2];
  }
}
