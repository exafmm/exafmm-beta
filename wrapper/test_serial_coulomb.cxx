#include <cmath>
#include <cstdlib>
#include <iostream>
#include "mr3.h"

const double R2MIN = 0.25;
const double R2MAX = 64;

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double rscale, int tblno, double size, int periodicflag);

int main() {
  const int N = 10000;
  const double size = 2;
  double *xi     = new double [3*N];
  double *qi     = new double [N];
  double *pi     = new double [3*N];
  double *fi     = new double [3*N];
  double *pd     = new double [3*N];
  double *fd     = new double [3*N];
  double *xj     = new double [3*N];
  double *qj     = new double [N];

  srand48(0);
  float average = 0;
  for( int i=0; i!=N; ++i ) {
    xi[3*i+0] = drand48() * size - size/2;
    xi[3*i+1] = drand48() * size - size/2;
    xi[3*i+2] = drand48() * size - size/2;
    qi[i] = drand48()*2.0-1.0;
    average += qi[i];
    pi[3*i+0] = pi[3*i+1] = pi[3*i+2] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[3*i+0] = pd[3*i+1] = pd[3*i+2] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
  }
  average /= N;
  for( int i=0; i!=N; ++i ) {
    qi[i] -= average;
  }
  average = 0;
  for( int i=0; i!=N; ++i ) {
    xj[3*i+0] = drand48() * size - size/2;
    xj[3*i+1] = drand48() * size - size/2;
    xj[3*i+2] = drand48() * size - size/2;
    qj[i] = drand48()*2.0-1.0;
  }
  average /= N;
  for( int i=0; i!=N; ++i ) {
    qj[i] -= average;
  }

  FMMcalccoulomb_ij(N, xi, qi, pi, N, xj, qj, 0.0, 1, size, 0);
  FMMcalccoulomb_ij(N, xi, qi, fi, N, xj, qj, 0.0, 0, size, 0);
#if 0
  MR3calccoulomb_ij(N, xi, qi, pd, N, xj, qj, 1.0, 1, size, 0);
  MR3calccoulomb_ij(N, xi, qi, fd, N, xj, qj, 1.0, 0, size, 0);
#else
  for( int i=0; i!=N; ++i ) {
    double P = 0, Fx = 0, Fy = 0, Fz = 0;
    for( int j=0; j!=N; ++j ) {
      double dx = xi[3*i+0] - xj[3*j+0];
      double dy = xi[3*i+1] - xj[3*j+1];
      double dz = xi[3*i+2] - xj[3*j+2];
      double R2 = dx * dx + dy * dy + dz * dz;
      double invR = 1 / std::sqrt(R2);
      if( R2 == 0 ) invR = 0;
      double invR3 = qj[j] * invR * invR * invR;
      P += qj[j] * invR;
      Fx += dx * invR3;
      Fy += dy * invR3;
      Fz += dz * invR3;
    }
    pd[3*i+0] += P;
    fd[3*i+0] += Fx;
    fd[3*i+1] += Fy;
    fd[3*i+2] += Fz;
  }
#endif
  double Pd = 0, Pn = 0, Fd = 0, Fn = 0;
  float sums = 0;
  for( int i=0; i!=N; ++i ) {
    Pd += (pi[3*i+0] - pd[3*i+0]) * (pi[3*i+0] - pd[3*i+0]);
    Pn += pd[3*i+0] * pd[3*i+0];
    Fd += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
        + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
        + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    Fn += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
    sums += qi[i]*fd[3*i+0];
  }
  std::cout << "Coulomb       potential : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Coulomb       force     : " << sqrtf(Fd/Fn) << std::endl;
  std::cout << sums << std::endl;

  delete[] xi;
  delete[] qi;
  delete[] pi;
  delete[] fi;
  delete[] pd;
  delete[] fd;
  delete[] xj;
  delete[] qj;
}
