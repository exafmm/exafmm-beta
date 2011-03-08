#include <sys/time.h>
#include "kernel.h"
#include "gpu.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const real EPS = 1e-6;

double get_gpu_time() {
  cudaThreadSynchronize();
  struct timeval tv;                                            // Time value
  gettimeofday(&tv, NULL);                                      // Get time of day in seconds and microseconds
  return double(tv.tv_sec+tv.tv_usec*1e-6);                     // Combine seconds and microseconds and return
}

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
    fact = fact + 2;
  }
}

__device__ void evalLocal(float *YnmShrd, float rho, float alpha, float *factShrd) {
  float x = cosf(alpha);
  float s = sqrt(1 - x * x);
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
    fact = fact + 2;
  }
}

__device__ void P2M_core(float *target, float rho, float alpha, float beta, float source) {
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[NCOEF];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int nn = floor(sqrt(2*threadIdx.x+0.25)-0.5);
  int mm = 0;
  for( int i=0; i<=nn; ++i ) mm += i;
  mm = threadIdx.x - mm;
  if( threadIdx.x >= NCOEF ) nn = mm = 0;
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
  fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<=mm; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    YnmShrd[i] = rhom * p * rsqrtf(factShrd[2*m]);
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho;
    float rhon = rhom;
    for( int n=m+1; n<=nn; ++n ) {
      i = n * (n + 1) / 2 + m;
      YnmShrd[i] = rhon * p * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= rho;
    }
    pn = -pn * fact * s;
    fact = fact + 2;
  }
  int i = nn * (nn + 1) / 2 + mm;
  float ere = cosf(-mm * beta);
  float eim = sinf(-mm * beta);
  target[0] += source * YnmShrd[i] * ere;
  target[1] += source * YnmShrd[i] * eim;
}

__global__ void P2M_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float targetShrd[3];
  __shared__ float sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS;
  targetShrd[0] = targetGlob[2*itarget+0];
  targetShrd[1] = targetGlob[2*itarget+1];
  targetShrd[2] = targetGlob[2*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    int size  = rangeGlob[keys+2*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
      __syncthreads();
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = sourceShrd[4*i+0] - targetShrd[0];
        d.y = sourceShrd[4*i+1] - targetShrd[1];
        d.z = sourceShrd[4*i+2] - targetShrd[2];
        float rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        P2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
    sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
    sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
    sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = sourceShrd[4*i+0] - targetShrd[0];
      d.y = sourceShrd[4*i+1] - targetShrd[1];
      d.z = sourceShrd[4*i+2] - targetShrd[2];
      float rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      P2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void M2M_core(float *target, float beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrt(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NCOEF ) j = k = 0;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<=j; ++n ) {
    for( int m=-n; m<=min(k-1,n); ++m ) {
      if( j-n >= k-m ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
        float ere = cosf(-m * beta);
        float eim = sinf(-m * beta);
        float ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        float cnm = ODDEVEN((m-abs(m))/2+j);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] -= sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] += sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
    for( int m=k; m<=n; ++m ) {
      if( j-n >= m-k ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
        float ere = cosf(-m * beta);
        float eim = sinf(-m * beta);
        float ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        float cnm = ODDEVEN(k+j+m);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] += sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
  }
}

__global__ void M2M_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[P*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    M2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void M2L_core(float *target, float  beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrt(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NCOEF ) j = k = 0;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=-n; m<0; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = anm * ajk * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 - m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] += sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] -= sourceShrd[2*i+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = ODDEVEN((abs(k - m) - k - m) / 2);
      cnm *= anm * ajk * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 + m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] -= sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] += sourceShrd[2*i+1] * CnmReal;
    }
  }
}

__global__ void M2L_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[4*NCOEF];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalLocal(YnmShrd,rho,alpha,factShrd);
    M2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void M2P_core(float &pot, float r, float theta, float phi, float *factShrd, float *sourceShrd) {
  float x = cosf(theta);
  float s = sqrtf(1 - x * x);
  float fact = 1;
  float pn = 1;
  float rhom = 1.0 / r;
  for( int m=0; m<P; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    float ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    float eim = sinf(m * phi);
    float anm = rhom * rsqrt(factShrd[2*m]);
    float Ynm = anm * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    float realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    pot += 2 * Ynm * realj;
    rhom /= r;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrt(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      pot += 2 * Ynm * realj;
      rhon /= r;
    }
    pn = -pn * fact * s;
    fact = fact+2;
  }
}

__global__ void M2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[4];
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  target[0] = targetGlob[4*itarget+0];
  target[1] = targetGlob[4*itarget+1];
  target[2] = targetGlob[4*itarget+2];
  target[3] = targetGlob[4*itarget+3];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    float3 d;
    d.x = target[0] - sourceGlob[begin+0];
    d.y = target[1] - sourceGlob[begin+1];
    d.z = target[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    M2P_core(target[3],r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

__global__ void P2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[4];
  __shared__ float sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  target[0] = targetGlob[4*itarget+0];
  target[1] = targetGlob[4*itarget+1];
  target[2] = targetGlob[4*itarget+2];
  target[3] = targetGlob[4*itarget+3];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    int size  = rangeGlob[keys+2*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
      __syncthreads();
#pragma unroll 64
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = target[0] - sourceShrd[4*i+0];
        d.y = target[1] - sourceShrd[4*i+1];
        d.z = target[2] - sourceShrd[4*i+2];
        target[3] += sourceShrd[4*i+3] * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
    sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
    sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
    sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = target[0] - sourceShrd[4*i+0];
      d.y = target[1] - sourceShrd[4*i+1];
      d.z = target[2] - sourceShrd[4*i+2];
      target[3] += sourceShrd[4*i+3] * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
    }
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

__device__ void L2L_core(float *target, float beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrt(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NCOEF ) j = k = 0;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=j+k-n; m<0; ++m ) {
      int nms = n * (n + 1) / 2 - m;
      int jnkm = (n - j) * (n - j) + n - j + m - k;
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = ODDEVEN(k-n) * ajk / anm * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      target[0] += sourceShrd[2*nms+0] * CnmReal;
      target[0] += sourceShrd[2*nms+1] * CnmImag;
      target[1] += sourceShrd[2*nms+0] * CnmImag;
      target[1] -= sourceShrd[2*nms+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      if( n-j >= abs(m-k) ) {
        int nms = n * (n + 1) / 2 + m;
        int jnkm = (n - j) * (n - j) + n - j + m - k;
        float ere = cosf((m - k) * beta);
        float eim = sinf((m - k) * beta);
        float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
        float cnm = ODDEVEN((m-k-abs(m-k)) / 2 - n);
        cnm *= ajk / anm * YnmShrd[jnkm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*nms+0] * CnmReal;
        target[0] -= sourceShrd[2*nms+1] * CnmImag;
        target[1] += sourceShrd[2*nms+0] * CnmImag;
        target[1] += sourceShrd[2*nms+1] * CnmReal;
      }
    }
  }
}

__global__ void L2L_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[P*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    L2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void L2P_core(float &pot, float r, float theta, float phi, float *factShrd, float *sourceShrd) {
  float x = cosf(theta);
  float s = sqrtf(1 - x * x);
  float fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<P; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    float ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    float eim = sinf(m * phi);
    float anm = rhom * rsqrt(factShrd[2*m]);
    float Ynm = anm * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    float realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    pot += 2 * Ynm * realj;
    rhom *= r;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrt(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      pot += 2 * Ynm * realj;
      rhon *= r;
    }
    pn = -pn * fact * s;
    fact = fact + 2;
  }
}

__global__ void L2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[4];
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  target[0] = targetGlob[4*itarget+0];
  target[1] = targetGlob[4*itarget+1];
  target[2] = targetGlob[4*itarget+2];
  target[3] = targetGlob[4*itarget+3];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    float3 d;
    d.x = target[0] - sourceGlob[begin+0];
    d.y = target[1] - sourceGlob[begin+1];
    d.z = target[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    L2P_core(target[3],r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

void Kernel::initialize() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
}

void Kernel::allocGPU() {
  cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(float) );
  cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(float) );
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(float),cudaMemcpyHostToDevice);
}

void Kernel::deallocGPU() {
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(float),cudaMemcpyDeviceToHost);
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
}

void Kernel::P2M() {
  double tic,toc;
  bool p = true;
  if( MPIRANK != 0 ) p = false;
  if(p) std::cout << "-----------" << std::endl;
  tic = get_gpu_time();
  allocGPU();
  toc = get_gpu_time();
  if(p) std::cout << "Alloc GPU     : " << toc-tic << std::endl;
  tic = get_gpu_time();
  int numBlocks = keysHost.size();
  P2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  toc = get_gpu_time();
  if(p) std::cout << "GPU kernel    : " << toc-tic << std::endl;
  tic = get_gpu_time();
  deallocGPU();
  toc = get_gpu_time();
  if(p) std::cout << "Dealloc GPU   : " << toc-tic << std::endl;
  if(p) std::cout << "-----------" << std::endl;
}

void Kernel::M2M() {
  allocGPU();
  int numBlocks = keysHost.size();
  M2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::M2L() {
  allocGPU();
  int numBlocks = keysHost.size();
  M2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::M2P() {
  allocGPU();
  int numBlocks = keysHost.size();
  M2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::P2P() {
  allocGPU();
  int numBlocks = keysHost.size();
  P2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::L2L() {
  allocGPU();
  int numBlocks = keysHost.size();
  L2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::L2P() {
  allocGPU();
  int numBlocks = keysHost.size();
  L2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  deallocGPU();
}

void Kernel::finalize() {}
