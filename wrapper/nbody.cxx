#include <cmath>
#include <cstdlib>
#include <iostream>

extern "C" void FMMcalccoulomb(double* x, int n, double* q, double rscale, int tblno,
  double xmax, int periodicflag, int natchangeflag, double* force);

int main() {
  const int N = 10000;
  const double eps2 = 1e-4;
  const double xmax = 2;
  double *x = new double [3*N];
  double *q = new double [N];
  double *f = new double [3*N];

  for( int i=0; i!=N; ++i ) {
    x[3*i+0] = rand() / (1. + RAND_MAX) * xmax;
    x[3*i+1] = rand() / (1. + RAND_MAX) * xmax;
    x[3*i+2] = rand() / (1. + RAND_MAX) * xmax;
    q[i] = 1. / N;
    f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
  }

  FMMcalccoulomb(x, N, q, 0.0, 0, xmax, 0, 0, f);
  double Fd = 0, Fn = 0;
  for( int i=0; i!=N; ++i ) {
    double Fx = 0, Fy = 0, Fz = 0;
    for( int j=0; j!=N; ++j ) {
      double dx = x[3*i+0] - x[3*j+0];
      double dy = x[3*i+1] - x[3*j+1];
      double dz = x[3*i+2] - x[3*j+2];
      double invR = 1 / sqrtf(dx * dx + dy * dy + dz * dz + eps2);
      double invR3 = q[j] * invR * invR * invR;
      Fx += dx * invR3;
      Fy += dy * invR3;
      Fz += dz * invR3;
    }
    Fd += (f[3*i+0] - Fx) * (f[3*i+0] - Fx)
        + (f[3*i+1] - Fy) * (f[3*i+1] - Fy)
        + (f[3*i+2] - Fz) * (f[3*i+2] - Fz);
    Fn += Fx * Fx + Fy * Fy + Fz * Fz;
  }
  std::cout << "Error (acc)   : " << sqrtf(Fd/Fn) << std::endl;

  delete[] x;
  delete[] q;
  delete[] f;
}
