#include "construct.h"

void gaussian(int numBodies, float *xj, float *yj, float *zj, float *charge, float *s,
              float *xi, float *yi, float *zi, float *val) {
  Bodies bodies(numBodies);
  Bodies jbodies = bodies;
  Cells cells,jcells;
  TreeConstructor T;
  T.printNow = true;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->X[0]  = xi[B-bodies.begin()];
    B->X[1]  = yi[B-bodies.begin()];
    B->X[2]  = zi[B-bodies.begin()];
    B->val   = val[B-bodies.begin()];
    B->IBODY = B-bodies.begin();
  }
  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
    B->X[0]  = xj[B-jbodies.begin()];
    B->X[1]  = yj[B-jbodies.begin()];
    B->X[2]  = zj[B-jbodies.begin()];
    B->Q     = charge[B-jbodies.begin()];
    B->S     = s[B-jbodies.begin()];
  }

  T.setDomain(bodies);
  T.bottomup(bodies,cells);
  T.bottomup(jbodies,jcells);
  T.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    xi[B-bodies.begin()]  = B->X[0];
    yi[B-bodies.begin()]  = B->X[1];
    zi[B-bodies.begin()]  = B->X[2];
    val[B-bodies.begin()] = B->val;
  }
}
