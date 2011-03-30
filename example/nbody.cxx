#include <cmath>
#include <cstdlib>
#include <iostream>

extern void laplace(int, float*, float*, float*, float*, float*, float*, float*, float*);

int main() {
  const int N = 10000;
  const float eps2 = 1e-4;
  float *x  = new float [N];
  float *y  = new float [N];
  float *z  = new float [N];
  float *m  = new float [N];
  float *p  = new float [N];
  float *fx = new float [N];
  float *fy = new float [N];
  float *fz = new float [N];

  for( int i=0; i!=N; ++i ) {
    x[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    y[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    z[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    m[i] = 1. / N;
    p[i] = -m[i] / sqrtf(eps2);
    fx[i] = fy[i] = fz[i] = 0;
  }

  laplace(N,x,y,z,m,p,fx,fy,fz);

  float Pd = 0, Pn = 0, Fd = 0, Fn = 0;
  for( int i=0; i!=N; ++i ) {
    float P = -m[i] / sqrtf(eps2);
    float Fx = 0, Fy = 0, Fz = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = x[i] - x[j];
      float dy = y[i] - y[j];
      float dz = z[i] - z[j];
      float invR = 1 / sqrtf(dx * dx + dy * dy + dz * dz + eps2);
      float invR3 = m[j] * invR * invR * invR;
      P += m[j] * invR;
      Fx -= dx * invR3;
      Fy -= dy * invR3;
      Fz -= dz * invR3;
    }
    Pd += (p[i] - P) * (p[i] - P);
    Pn += P * P;
    Fd += (fx[i] - Fx) * (fx[i] - Fx) + (fy[i] - Fy) * (fy[i] - Fy) + (fz[i] - Fz) * (fz[i] - Fz);
    Fn += Fx * Fx + Fy * Fy + Fz * Fz;
  }
  std::cout << "Error (pot)   : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Error (acc)   : " << sqrtf(Fd/Fn) << std::endl;

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] m;
  delete[] p;
  delete[] fx;
  delete[] fy;
  delete[] fz;
}
