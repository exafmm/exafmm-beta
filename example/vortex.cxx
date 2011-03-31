#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

extern void biotsavart(int, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
extern void stretching(int, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*);
extern void gaussian(int, float*, float*, float*, float*, float*, float*, float*, float*, float*);

int main() {
  const int N = 1000;
  const float eps2 = 1e-4;
  float *xi = new float [N];
  float *yi = new float [N];
  float *zi = new float [N];
  float *xj = new float [N];
  float *yj = new float [N];
  float *zj = new float [N];
  float *qx = new float [N];
  float *qy = new float [N];
  float *qz = new float [N];
  float *s  = new float [N];
  float *u  = new float [N];
  float *v  = new float [N];
  float *w  = new float [N];

  for( int i=0; i!=N; ++i ) {
    xi[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    yi[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    zi[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    xj[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    yj[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    zj[i] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;
    qx[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    qy[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    qz[i] = (rand() / (1. + RAND_MAX) * 2 - 1) / N;
    s[i] = 2 * powf(N,-1./3);
    u[i] = v[i] = w[i] = 0;
  }

  biotsavart(N,xi,yi,zi,qx,qy,qz,s,u,v,w);

  float Ud = 0, Un = 0;
  for( int i=0; i!=N; ++i ) {
    float U = 0, V = 0, W = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = xi[i] - xi[j];
      float dy = yi[i] - yi[j];
      float dz = zi[i] - zi[j];
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

  stretching(N,xi,yi,zi,qx,qy,qz,s,u,v,w);

  float Qd = 0, Qn = 0;
  for( int i=0; i!=N; ++i ) {
    float U = 0, V = 0, W = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = xi[i] - xi[j];
      float dy = yi[i] - yi[j];
      float dz = zi[i] - zi[j];
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
    u[i] = v[i] = w[i] = 0;
  }
  std::cout << "Error (ST)    : " << sqrtf(Qd/Qn) << std::endl;

  gaussian(N,xj,yj,zj,qx,s,xi,yi,zi,v);

  float Vd = 0, Vn = 0;
  for( int i=0; i!=N; ++i ) {
    float V = 0;
    for( int j=0; j!=N; ++j ) {
      float dx = xi[i] - xj[j];
      float dy = yi[i] - yj[j];
      float dz = zi[i] - zj[j];
      float S2 = 2 * s[j] * s[j];
      float R2 = dx * dx + dy * dy + dz * dz + eps2;
      V += qx[j] / (M_PI * S2) / std::sqrt(M_PI * S2) * exp(-R2 / S2);
    }
    Vd += (v[i] - V) * (v[i] - V);
    Vn += V * V;
  }
  std::cout << "Error (GA)    : " << sqrtf(Vd/Vn) << std::endl;

  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] xj;
  delete[] yj;
  delete[] zj;
  delete[] qx;
  delete[] qy;
  delete[] qz;
  delete[] s;
  delete[] u;
  delete[] v;
  delete[] w;
}
