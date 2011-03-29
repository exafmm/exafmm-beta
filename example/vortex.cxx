#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

extern void biotsavart(int, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
extern void stretching(int, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);

int main() {
  const int N = 10;
  const float eps2 = 1e-4;
  float *x  = new float [N];
  float *y  = new float [N];
  float *z  = new float [N];
  float *qx = new float [N];
  float *qy = new float [N];
  float *qz = new float [N];
  float *s  = new float [N];
  float *u  = new float [N];
  float *v  = new float [N];
  float *w  = new float [N];

  for( int i=0; i!=N; ++i ) {
    x[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    y[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    z[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    qx[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    qy[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    qz[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    s[i] = 2 * powf(N,-1./3);
    u[i] = v[i] = w[i] = 0;
  }

  biotsavart(N,x,y,z,qx,qy,qz,s,u,v,w);

  float Ud = 0, Un = 0;
  for( int i=0; i!=N; ++i ) {
    float U = 0, V = 0, W = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = x[i] - x[j];
      float dy = y[i] - y[j];
      float dz = z[i] - z[j];
      float S2 = 2 * s[j] * s[j];
      float R2 = dx * dx + dy * dy + dz * dz + eps2;
      float RS = R2 / S2;
      float cutoff = 0.25 / M_PI / R2 / sqrtf(R2) * (erf( sqrtf(RS) )
                   - sqrtf(4 / M_PI * RS) * exp(-RS));
      U += (dy * qz[j] - dz * qy[j]) * cutoff;
      V += (dz * qx[j] - dx * qz[j]) * cutoff;
      W += (dx * qy[j] - dy * qx[j]) * cutoff;
    }
    Ud += (u[i] - U) * (u[i] - U) + (v[i] - V) * (v[i] - V) + (w[i] - W) * (w[i] - W);
    Un += U * U + V * V + W * W;
    u[i] = v[i] = w[i] = 0;
  }
  std::cout << "Error (BS)    : " << sqrtf(Ud/Un) << std::endl;

  stretching(N,x,y,z,qx,qy,qz,s,u,v,w);

  float Qd = 0, Qn = 0;
  for( int i=0; i!=N; ++i ) {
    float U = 0, V = 0, W = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = x[i] - x[j];
      float dy = y[i] - y[j];
      float dz = z[i] - z[j];
      float S2 = 2 * s[j] * s[j];
      float R2 = dx * dx + dy * dy + dz * dz + eps2;
      float RS = R2 / S2;
      float cutoff = 0.25 / M_PI / R2 / sqrtf(R2) * (erf( sqrtf(RS) )
                   - sqrtf(4 / M_PI * RS) * exp(-RS));
      U += (qy[i] * qz[j] - qz[i] * qy[j]) * cutoff;
      V += (qz[i] * qx[j] - qx[i] * qz[j]) * cutoff;
      W += (qx[i] * qy[j] - qy[i] * qx[j]) * cutoff;
      cutoff = 0.25 / M_PI / R2 / R2 / sqrtf(R2) * (3 * erf( sqrtf(RS) )
             - (2 * RS + 3) * sqrtf(4 / M_PI * RS) * exp(-RS))
             * (qx[i] * dx + qy[i] * dy + qz[i] * dz);
      U += (qy[j] * dz - qz[j] * dy) * cutoff;
      V += (qz[j] * dx - qx[j] * dz) * cutoff;
      W += (qx[j] * dy - qy[j] * dx) * cutoff;
    }
    Qd += (u[i] - U) * (u[i] - U) + (v[i] - V) * (v[i] - V) + (w[i] - W) * (w[i] - W);
    Qn += U * U + V * V + W * W;
  }
  std::cout << "Error (ST)    : " << sqrtf(Qd/Qn) << std::endl;

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] qx;
  delete[] qy;
  delete[] qz;
  delete[] s;
  delete[] u;
  delete[] v;
  delete[] w;
}
