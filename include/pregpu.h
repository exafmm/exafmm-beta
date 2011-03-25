#ifndef pregpu_h
#define pregpu_h
#include "cuprintf.h"

int   *keysDevc;                                                // Keys on device
int   *rangeDevc;                                               // Ranges on device
float *sourceDevc;                                              // Sources on device
float *targetDevc;                                              // Targets on device
__device__ __constant__ float constDevc[1];                     // Constants on device

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
  float fact = 1;
  float pn = 1;
  float rhom = 1.0 / rho;
  for( int m=0; m<2*P; ++m ){
    float p = pn;
    int i = m * (m + 1) /2 + m;
    YnmShrd[i] = rhom * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom /= rho;
    float rhon = rhom;
    for( int n=m+1; n<2*P; ++n ){
      i = n * (n + 1) / 2 + m;
      YnmShrd[i] = rhon * p * factShrd[n-m];
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon /= rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
}

void Kernel::initialize() {
  precalculate();
  startTimer("Init GPU     ");                                  // Start timer
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
  cudaPrintfInit();
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel::allocGPU() {
  cudaThreadSynchronize();
  startTimer("cudaMalloc   ");
  cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(float) );
  cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(float) );
  cudaThreadSynchronize();
  stopTimer("cudaMalloc   ");
  startTimer("cudaMemcpy   ");
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(float));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

void Kernel::deallocGPU() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(float),cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
  startTimer("cudaFree     ");
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
  cudaThreadSynchronize();
  stopTimer("cudaFree     ");
}

void Kernel::finalize() {
  postcalculate();
  cudaPrintfDisplay(stdout, true);
  cudaPrintfEnd();
}

#endif
