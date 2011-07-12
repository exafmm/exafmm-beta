#include "construct.h"

extern "C" void FMMcalccoulomb (double* x, int n, double* q, double rscale, int tblno,
  double xmax, int periodicflag, int natchangeflag, double* force) {
  int nimages = ((periodicflag & 0x1) == 0) ? 0 : 3;
  IMAGES = nimages;
  THETA = 1/sqrtf(3);
  vect shift = xmax/2;
  Bodies bodies(n);
  Cells cells,jcells;
  TreeConstructor T;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0]   = x[3*i+0];
    B->X[1]   = x[3*i+1];
    B->X[2]   = x[3*i+2];
    B->SRC[0] = q[i];
    switch (tblno) {
    case 0 :
      B->TRG[1] = -force[3*i+0];
      B->TRG[2] = -force[3*i+1];
      B->TRG[3] = -force[3*i+2];
      break;
    case 1 :
      B->TRG[0] = force[3*i+0];
      break;
    }
    B->IBODY = i;
  }

  T.setKernel("Laplace");
  T.setDomain(bodies,shift,xmax/2);
  T.bottomup(bodies,cells);
  jcells = cells;
  T.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i]     = B->SRC[0];
    switch (tblno) {
    case 0 :
      force[3*i+0] = -B->TRG[1];
      force[3*i+1] = -B->TRG[2];
      force[3*i+2] = -B->TRG[3];
      break;
    case 1 :
      force[3*i+0] = B->TRG[0];
      break;
    }
  }
}
