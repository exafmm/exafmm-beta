#ifndef pregpu_h
#define pregpu_h
#include <omp.h>

static size_t keysDevcSize = 0;                                 // Size of offsets for rangeHost
static size_t rangeDevcSize = 0;                                // Size of offsets for sourceHost
static size_t sourceDevcSize = 0;                               // Size of sources
static size_t targetDevcSize = 0;                               // Size of targets
static int    *keysDevc;                                        // Keys on device
static int    *rangeDevc;                                       // Ranges on device
static float  *sourceDevc;                                      // Sources on device
static float  *targetDevc;                                      // Targets on device
#pragma omp threadprivate(keysDevcSize,rangeDevcSize,sourceDevcSize,targetDevcSize)
#pragma omp threadprivate(keysDevc,rangeDevc,sourceDevc,targetDevc)
__device__ __constant__ float constDevc[1];                     // Constants on device

namespace {
__device__ void cart2sph(float& r, float& theta, float& phi, float dx, float dy, float dz) {
  r = sqrtf(dx * dx + dy * dy + dz * dz)+EPS;
  theta = acosf(dz / r);
  if( fabs(dx) + fabs(dy) < EPS ) {
    phi = 0;
  } else if( fabs(dx) < EPS ) {
    phi = dy / fabs(dy) * M_PI * 0.5;
  } else if( dx > 0 ) {
    phi = atanf(dy / dx);
  } else {
    phi = atanf(dy / dx) + M_PI;
  }
}

__device__ void sph2cart(float r, float theta, float phi, float *spherical, float *cartesian) {
  cartesian[0] = sinf(theta) * cosf(phi) * spherical[0]
               + cosf(theta) * cosf(phi) / r * spherical[1]
               - sinf(phi) / r / sinf(theta) * spherical[2];
  cartesian[1] = sinf(theta) * sinf(phi) * spherical[0]
               + cosf(theta) * sinf(phi) / r * spherical[1]
               + cosf(phi) / r / sinf(theta) * spherical[2];
  cartesian[2] = cosf(theta) * spherical[0]
               - sinf(theta) / r * spherical[1];
}

__device__ void evalMultipole(float *YnmShrd, float rho, float alpha, float *factShrd) {
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
  float fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<P; ++m ){
    float p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    YnmShrd[npn] = rhom * p / factShrd[2*m];
    YnmShrd[nmn] = YnmShrd[npn];
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= -rho;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ){
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      YnmShrd[npm] = rhon * p / factShrd[n+m];
      YnmShrd[nmm] = YnmShrd[npm];
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= -rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
}

__device__ void evalLocal(float *YnmShrd, float rho, float alpha, float *factShrd) {
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
  float rho_1 = 1 / rho;
  for( int l=threadIdx.x; l<(2*P+1)*P; l+=THREADS ){
    float fact = 1;
    float pn = 1;
    float rhom = rho_1;
    int nn = floor(sqrtf(2*l+0.25)-0.5);
    int mm = 0;
    float Ynm;
    for( int i=0; i<=nn; ++i ) mm += i;
    mm = l - mm;
    int n;
    for( int m=0; m<mm; ++m ){
      rhom *= rho_1;
      pn = -pn * fact * s;
      fact += 2;
    }
    int m = mm;
    float p = pn;
    if( mm == nn ) Ynm = rhom * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho_1;
    float rhon = rhom;
    for( n=m+1; n<nn; ++n ){
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= rho_1;
    }
    if( n <= nn ) Ynm = rhon * p * factShrd[n-m];
    YnmShrd[l] = Ynm;
  }
  __syncthreads();
}
}

#endif
