/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#include "cuprintf.h"
#define KERNEL
#include "kernel.h"
#undef KERNEL
#include "laplace.h"
#include "pregpu.h"

template<>
void Kernel<Laplace>::initialize() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(DEVICE);                                        // Set GPU device
  cudaPrintfInit();                                             // Initialize cuPrintf
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

template<>
void Kernel<Laplace>::M2M_CPU() {
  const complex I(0.,1.);                                       // Imaginary unit
  vect dist = CI->X - CJ->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,-beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex M = 0;
      for( int n=0; n<=j; ++n ) {
        for( int m=-n; m<=std::min(k-1,n); ++m ) {
          if( j-n >= k-m ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            const int nm    = n * n + n + m;
            M += CJ->M[3*jnkms] * std::pow(I,double(m-abs(m))) * Ynm[nm]
               * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            M += std::conj(CJ->M[3*jnkms]) * Ynm[nm]
               * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
      }
      CI->M[3*jks] += M;
    }
  }
}

template<>
void Kernel<Laplace>::finalize() {
  if( MPIRANK == 0 ) cudaPrintfDisplay(stdout, true);           // Print cuPrintf buffer to display
  cudaPrintfEnd();                                              // Finalize cuPrintf
}

__device__ void LaplaceP2M_core(gpureal *target, gpureal rho, gpureal alpha, gpureal beta, gpureal source) {
  __shared__ gpureal factShrd[2*P];
  gpureal Ynm;
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int nn = floorf(sqrtf(2*threadIdx.x+0.25)-0.5);
  int mm = 0;
  for( int i=0; i<=nn; ++i ) mm += i;
  mm = threadIdx.x - mm;
  if( threadIdx.x >= NTERM ) nn = mm = 0;
  gpureal x = cosf(alpha);
  gpureal s = sqrtf(1 - x * x);
  fact = 1;
  gpureal pn = 1;
  gpureal rhom = 1;
  for( int m=0; m<mm; ++m ) {
    rhom *= rho;
    pn = -pn * fact * s;
    fact += 2;
  }
  int m=mm;
  gpureal p = pn;
  if(mm==nn) Ynm = rhom * p * rsqrtf(factShrd[2*m]);
  gpureal p1 = p;
  p = x * (2 * m + 1) * p;
  rhom *= rho;
  gpureal rhon = rhom;
  for( int n=m+1; n<=nn; ++n ) {
    if(n==nn){
      Ynm = rhon * p * rsqrtf(factShrd[n+m] / factShrd[n-m]);
    }
    gpureal p2 = p1;
    p1 = p;
    p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
    rhon *= rho;
  }
  gpureal ere = cosf(-mm * beta);
  gpureal eim = sinf(-mm * beta);
  target[0] += source * Ynm * ere;
  target[1] += source * Ynm * eim;
}

__global__ void LaplaceP2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[2] = {0, 0};
  __shared__ gpureal targetShrd[3];
  __shared__ gpureal sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS;
  targetShrd[0] = targetGlob[6*itarget+0];
  targetShrd[1] = targetGlob[6*itarget+1];
  targetShrd[2] = targetGlob[6*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    int size  = rangeGlob[keys+3*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[7*isource+3];
      __syncthreads();
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = sourceShrd[4*i+0] - targetShrd[0];
        d.y = sourceShrd[4*i+1] - targetShrd[1];
        d.z = sourceShrd[4*i+2] - targetShrd[2];
        gpureal rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        LaplaceP2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    if( threadIdx.x < size - iblok * THREADS ) {
      sourceShrd[4*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[7*isource+3];
    }
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = sourceShrd[4*i+0] - targetShrd[0];
      d.y = sourceShrd[4*i+1] - targetShrd[1];
      d.z = sourceShrd[4*i+2] - targetShrd[2];
      gpureal rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      LaplaceP2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
}

__device__ void LaplaceM2M_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
  int j = floorf(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<=j; ++n ) {
    for( int m=-n; m<=min(k-1,n); ++m ) {
      if( j-n >= k-m ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
        gpureal ere = cosf(-m * beta);
        gpureal eim = sinf(-m * beta);
        gpureal ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        gpureal cnm = ODDEVEN((m-abs(m))/2+j);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        gpureal CnmReal = cnm * ere;
        gpureal CnmImag = cnm * eim;
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
        gpureal ere = cosf(-m * beta);
        gpureal eim = sinf(-m * beta);
        gpureal ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        gpureal cnm = ODDEVEN(k+j+m);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        gpureal CnmReal = cnm * ere;
        gpureal CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] += sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
  }
}

__global__ void LaplaceM2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[2] = {0, 0};
  __shared__ gpureal sourceShrd[2*THREADS];
  __shared__ gpureal factShrd[2*P];
  __shared__ gpureal YnmShrd[P*P];
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetGlob[6*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[6*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[6*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
    }
    __syncthreads();
    gpureal rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    LaplaceM2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
}

__device__ void LaplaceM2L_core(gpureal *target, gpureal  beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
  int j = floorf(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=-n; m<0; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
      gpureal ere = cosf((m - k) * beta);
      gpureal eim = sinf((m - k) * beta);
      gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      gpureal cnm = anm * ajk * YnmShrd[jnkm];
      gpureal CnmReal = cnm * ere;
      gpureal CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 - m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] += sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] -= sourceShrd[2*i+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
      gpureal ere = cosf((m - k) * beta);
      gpureal eim = sinf((m - k) * beta);
      gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      gpureal cnm = ODDEVEN((abs(k - m) - k - m) / 2);
      cnm *= anm * ajk * YnmShrd[jnkm];
      gpureal CnmReal = cnm * ere;
      gpureal CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 + m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] -= sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] += sourceShrd[2*i+1] * CnmReal;
    }
  }
}

__global__ void LaplaceM2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal target[2] = {0, 0};
  __shared__ gpureal sourceShrd[2*THREADS];
  __shared__ gpureal factShrd[2*P];
  __shared__ gpureal YnmShrd[4*NTERM];
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
    }
    __syncthreads();
    int I = 0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz, ++I ) {
          if( Iperiodic & (1 << I) ) {
            float3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            d.x += targetGlob[6*itarget+0] - sourceGlob[begin+0];
            d.y += targetGlob[6*itarget+1] - sourceGlob[begin+1];
            d.z += targetGlob[6*itarget+2] - sourceGlob[begin+2];
            gpureal rho,alpha,beta;
            cart2sph(rho,alpha,beta,d.x,d.y,d.z);
            evalLocal(YnmShrd,rho,alpha,factShrd);
            LaplaceM2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
          }
        }
      }
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
}

__device__ void LaplaceM2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd) {
  gpureal x = cosf(theta);
  gpureal y = sinf(theta);
  if( fabsf(y) < EPS ) y = 1 / EPS;
  gpureal s = sqrtf(1 - x * x);
  gpureal spherical[3] = {0, 0, 0};
  gpureal cartesian[3] = {0, 0, 0};
  gpureal fact = 1;
  gpureal pn = 1;
  gpureal rhom = 1.0 / r;
  for( int m=0; m<P; ++m ) {
    gpureal p = pn;
    int i = m * (m + 1) / 2 + m;
    gpureal ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    gpureal eim = sinf(m * phi);
    gpureal anm = rhom * rsqrtf(factShrd[2*m]);
    gpureal Ynm = anm * p;
    gpureal p1 = p;
    p = x * (2 * m + 1) * p;
    gpureal YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    gpureal realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    gpureal imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom /= r;
    gpureal rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      gpureal p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
      target[0] += 2 * Ynm * realj;
      spherical[0] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      rhon /= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,spherical,cartesian);
  target[1] += cartesian[0];
  target[2] += cartesian[1];
  target[3] += cartesian[2];
}

__global__ void LaplaceM2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal targetX[3];
  gpureal target[4] = {0, 0, 0, 0};
  __shared__ gpureal sourceShrd[2*THREADS];
  __shared__ gpureal factShrd[2*P];
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
    }
    __syncthreads();
    int I = 0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz, ++I ) {
          if( Iperiodic & (1 << I) ) {
            float3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            d.x += targetX[0] - sourceGlob[begin+0];
            d.y += targetX[1] - sourceGlob[begin+1];
            d.z += targetX[2] - sourceGlob[begin+2];
            gpureal r,theta,phi;
            cart2sph(r,theta,phi,d.x,d.y,d.z);
            LaplaceM2P_core(target,r,theta,phi,factShrd,sourceShrd);
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
}

__device__ inline void LaplaceP2P_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[4*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[4*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[4*i+2];
  gpureal invR = rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
  gpureal invR3 = sourceShrd[4*i+3] * invR * invR * invR;
  target[0] += sourceShrd[4*i+3] * invR;
  target[1] -= d.x * invR3;
  target[2] -= d.y * invR3;
  target[3] -= d.z * invR3;
}

__global__ void LaplaceP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal targetX[3];
  gpureal target[4] = {0, 0, 0, 0};
  __shared__ gpureal sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int size      = rangeGlob[keys+3*ilist+2];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[7*isource+3];
      __syncthreads();
      int I = 0;
      for( int ix=-1; ix<=1; ++ix ) {
        for( int iy=-1; iy<=1; ++iy ) {
          for( int iz=-1; iz<=1; ++iz, ++I ) {
            if( Iperiodic & (1 << I) ) {
              float3 d;
              d.x = ix * D0;
              d.y = iy * D0;
              d.z = iz * D0;
#pragma unroll 64
              for( int i=0; i<THREADS; ++i ) {
                LaplaceP2P_core(target,targetX,sourceShrd,d,i);
              }
            }
          }
        }
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    if( threadIdx.x < size - iblok * THREADS ) {
      sourceShrd[4*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[7*isource+3];
    }
    __syncthreads();
    int I = 0;
    int icounter=0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz, ++I ) {
          if( Iperiodic & (1 << I) ) {
            icounter++;
            float3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            for( int i=0; i<size-iblok*THREADS; ++i ) {
              LaplaceP2P_core(target,targetX,sourceShrd,d,i);
            }
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
}

__device__ void LaplaceL2L_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
  int j = floorf(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  gpureal ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=j+k-n; m<0; ++m ) {
      int nms = n * (n + 1) / 2 - m;
      int jnkm = (n - j) * (n - j) + n - j + m - k;
      gpureal ere = cosf((m - k) * beta);
      gpureal eim = sinf((m - k) * beta);
      gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      gpureal cnm = ODDEVEN(k-n) * ajk / anm * YnmShrd[jnkm];
      gpureal CnmReal = cnm * ere;
      gpureal CnmImag = cnm * eim;
      target[0] += sourceShrd[2*nms+0] * CnmReal;
      target[0] += sourceShrd[2*nms+1] * CnmImag;
      target[1] += sourceShrd[2*nms+0] * CnmImag;
      target[1] -= sourceShrd[2*nms+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      if( n-j >= abs(m-k) ) {
        int nms = n * (n + 1) / 2 + m;
        int jnkm = (n - j) * (n - j) + n - j + m - k;
        gpureal ere = cosf((m - k) * beta);
        gpureal eim = sinf((m - k) * beta);
        gpureal anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
        gpureal cnm = ODDEVEN((m-k-abs(m-k)) / 2 - n);
        cnm *= ajk / anm * YnmShrd[jnkm];
        gpureal CnmReal = cnm * ere;
        gpureal CnmImag = cnm * eim;
        target[0] += sourceShrd[2*nms+0] * CnmReal;
        target[0] -= sourceShrd[2*nms+1] * CnmImag;
        target[1] += sourceShrd[2*nms+0] * CnmImag;
        target[1] += sourceShrd[2*nms+1] * CnmReal;
      }
    }
  }
}

__global__ void LaplaceL2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[2] = {0, 0};
  __shared__ gpureal sourceShrd[2*THREADS];
  __shared__ gpureal factShrd[2*P];
  __shared__ gpureal YnmShrd[P*P];
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetGlob[6*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[6*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[6*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
    }
    __syncthreads();
    gpureal rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    LaplaceL2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
}

__device__ void LaplaceL2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd) {
  gpureal x = cosf(theta);
  gpureal y = sinf(theta);
  if( fabsf(y) < EPS ) y = 1 / EPS;
  gpureal s = sqrtf(1 - x * x);
  gpureal spherical[3] = {0, 0, 0};
  gpureal cartesian[3] = {0, 0, 0};
  gpureal fact = 1;
  gpureal pn = 1;
  gpureal rhom = 1;
  for( int m=0; m<P; ++m ) {
    gpureal p = pn;
    int i = m * (m + 1) / 2 + m;
    gpureal ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    gpureal eim = sinf(m * phi);
    gpureal anm = rhom * rsqrtf(factShrd[2*m]);
    gpureal Ynm = anm * p;
    gpureal p1 = p;
    p = x * (2 * m + 1) * p;
    gpureal YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    gpureal realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    gpureal imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] += 2 * m / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom *= r;
    gpureal rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      gpureal p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
      target[0] += 2 * Ynm * realj;
      spherical[0] += 2 * n / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      rhon *= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,spherical,cartesian);
  target[1] += cartesian[0];
  target[2] += cartesian[1];
  target[3] += cartesian[2];
}

__global__ void LaplaceL2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal targetX[3];
  gpureal target[4] = {0, 0, 0, 0};
  __shared__ gpureal sourceShrd[2*THREADS];
  __shared__ gpureal factShrd[2*P];
  gpureal fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  __syncthreads();
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetX[0] - sourceGlob[begin+0];
    d.y = targetX[1] - sourceGlob[begin+1];
    d.z = targetX[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
    }
    __syncthreads();
    gpureal r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    LaplaceL2P_core(target,r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
}

#include "gpu.h"

CALL_GPU(Laplace,P2M,P2M GPUkernel);
CALL_GPU(Laplace,M2M,M2M GPUkernel);
CALL_GPU(Laplace,M2L,M2L GPUkernel);
CALL_GPU(Laplace,M2P,M2P GPUkernel);
CALL_GPU(Laplace,P2P,P2P GPUkernel);
CALL_GPU(Laplace,L2L,L2L GPUkernel);
CALL_GPU(Laplace,L2P,L2P GPUkernel);
