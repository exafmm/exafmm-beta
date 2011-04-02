#include "construct.h"

void biotsavart(int numBodies, float *x, float *y, float *z, float *qx, float *qy, float *qz, float *s,
                float *u, float *v, float *w) {
  Bodies bodies(numBodies);
  Cells cells,jcells;
  TreeConstructor T;
  T.printNow = false;

  T.startTimer("BiotSavart   ");
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->X[0]   = x[B-bodies.begin()];
    B->X[1]   = y[B-bodies.begin()];
    B->X[2]   = z[B-bodies.begin()];
    B->Q[0]   = qx[B-bodies.begin()];
    B->Q[1]   = qy[B-bodies.begin()];
    B->Q[2]   = qz[B-bodies.begin()];
    B->S      = s[B-bodies.begin()];
    B->vel[0] = u[B-bodies.begin()];
    B->vel[1] = v[B-bodies.begin()];
    B->vel[2] = w[B-bodies.begin()];
    B->IBODY  = B-bodies.begin();
  }

  T.setDomain(bodies);
  T.topdown(bodies,cells);
  jcells = cells;
  T.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    x[B-bodies.begin()]  = B->X[0];
    y[B-bodies.begin()]  = B->X[1];
    z[B-bodies.begin()]  = B->X[2];
    qx[B-bodies.begin()] = B->Q[0];
    qy[B-bodies.begin()] = B->Q[1];
    qz[B-bodies.begin()] = B->Q[2];
    s[B-bodies.begin()]  = B->S;
    u[B-bodies.begin()]  = B->vel[0];
    v[B-bodies.begin()]  = B->vel[1];
    w[B-bodies.begin()]  = B->vel[2];
  }
  T.stopTimer("BiotSavart   ",true);
}
