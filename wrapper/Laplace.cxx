#include "construct.h"

/*
extern "C" void MR3calccoulomb_fmm
(double* x, int n, double* q, double rscale, int tblno, double xmax, int periodicflag,
 int natchangeflag, double* force){

  int nimages = ((periodicflag & 0x1) == 0) ? 0 : 3;
  static int n_saved = 0;
  static double* x_saved[3] = {NULL,NULL,NULL};
  static double* f_saved[3] = {NULL,NULL,NULL};
  if (n > n_saved){
    n_saved = n * 2;
    if (x_saved[0] != NULL){
      for (int i = 0; i < 3; i++) free(x_saved[i]);
      for (int i = 0; i < 3; i++) free(f_saved[i]);
    }
    for (int i = 0; i < 3; i++) x_saved[i] = (double*)malloc(n_saved*sizeof(double));
    for (int i = 0; i < 3; i++) f_saved[i] = (double*)malloc(n_saved*sizeof(double));
  }
  double coef[3];
  for (int i = 0; i < 3; i++) coef[i] = 2*M_PI/xmax;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < 3; j++) x_saved[j][i] = x[3*i+j] * coef[j] - M_PI;
    for (int j = 0; j < 3; j++) f_saved[j][i] = 0.0;
  }
  laplace(n,x_saved[0],x_saved[1],x_saved[2],q,f_saved[0],f_saved[1],f_saved[2],nimages);
  double rcoef[3];
  for (int i = 0; i < 3; i++) rcoef[i] = coef[i] * coef[i];
  for (int i = 0; i < n; i++)
  for (int j = 0; j < 3; j++) force[3*i+j] = - f_saved[j][i] * rcoef[j];
}
*/

void laplace(int numBodies, float *x, float *y, float *z, float *charge,
             float *fx, float *fy, float *fz, int nimages) {
  IMAGES = nimages;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Cells cells,jcells;
  TreeConstructor T;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->X[0]   = x[B-bodies.begin()];
    B->X[1]   = y[B-bodies.begin()];
    B->X[2]   = z[B-bodies.begin()];
    B->SRC[0] = charge[B-bodies.begin()];
    B->TRG[1] = fx[B-bodies.begin()];
    B->TRG[2] = fy[B-bodies.begin()];
    B->TRG[3] = fz[B-bodies.begin()];
    B->IBODY = B-bodies.begin();
  }

  T.setKernel("Laplace");
  T.setDomain(bodies);
  T.bottomup(bodies,cells);
  jcells = cells;
  T.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    x[B-bodies.begin()]      = B->X[0];
    y[B-bodies.begin()]      = B->X[1];
    z[B-bodies.begin()]      = B->X[2];
    charge[B-bodies.begin()] = B->SRC[0];
    fx[B-bodies.begin()]     = B->TRG[1];
    fy[B-bodies.begin()]     = B->TRG[2];
    fz[B-bodies.begin()]     = B->TRG[3];
  }
}

