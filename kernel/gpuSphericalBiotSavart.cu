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
#define KERNEL
#include "kernel.h"
#undef KERNEL
#include "biotsavart.h"
#include "gpu.h"

template<>
void Kernel<BiotSavart>::initialize() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(DEVICE);                                        // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

template<>
void Kernel<BiotSavart>::M2M_CPU(C_iter CI, C_iter CJ) {
  const complex I(0.,1.);                                   // Imaginary unit
  vect dist = CI->X - CJ->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,-beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex M[3] = {0., 0., 0.};
      for( int n=0; n<=j; ++n ) {
        for( int m=-n; m<=std::min(k-1,n); ++m ) {
          if( j-n >= k-m ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            const int nm    = n * n + n + m;
            for( int d=0; d!=3; ++d ) {
              M[d] += CJ->M[3*jnkms+d] * std::pow(I,double(m-abs(m))) * Ynm[nm]
                    * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            for( int d=0; d!=3; ++d ) {
              M[d] += std::conj(CJ->M[3*jnkms+d]) * Ynm[nm]
                    * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
      }
      for( int d=0; d!=3; ++d ) {
        CI->M[3*jks+d] += M[d];
      }
    }
  }
}

template<>
void Kernel<BiotSavart>::finalize() {}

template<>
void Kernel<BiotSavart>::allocate() {
  cudaThreadSynchronize();
  startTimer("cudaMalloc   ");
  if( keysHost.size() > keysDevcSize ) {
    if( keysDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(keysDevc));
    CUDA_SAFE_CALL(cudaMalloc( (void**) &keysDevc, keysHost.size()*sizeof(int) ));
    keysDevcSize = keysHost.size();
  }
  if( rangeHost.size() > rangeDevcSize ) {
    if( rangeDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(rangeDevc));
    CUDA_SAFE_CALL(cudaMalloc( (void**) &rangeDevc, rangeHost.size()*sizeof(int) ));
    rangeDevcSize = rangeHost.size();
  }
  if( sourceHost.size() > sourceDevcSize ) {
    if( sourceDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(sourceDevc));
    CUDA_SAFE_CALL(cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) ));
    sourceDevcSize = sourceHost.size();
  }
  if( targetHost.size() > targetDevcSize ) {
    if( targetDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(targetDevc));
    CUDA_SAFE_CALL(cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) ));
    targetDevcSize = targetHost.size();
  }
  cudaThreadSynchronize();
  stopTimer("cudaMalloc   ");
}

template<>
void Kernel<BiotSavart>::hostToDevice() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  CUDA_SAFE_CALL(cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(gpureal)));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

template<>
void Kernel<BiotSavart>::deviceToHost() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

__device__ void BiotSavartP2M_core(gpureal *target, gpureal rho, gpureal alpha, gpureal beta,
                                   gpureal *sourceShrd, int ithread) {
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
  target[0] += sourceShrd[6*ithread+3] * Ynm * ere;
  target[1] += sourceShrd[6*ithread+3] * Ynm * eim;
  target[2] += sourceShrd[6*ithread+4] * Ynm * ere;
  target[3] += sourceShrd[6*ithread+4] * Ynm * eim;
  target[4] += sourceShrd[6*ithread+5] * Ynm * ere;
  target[5] += sourceShrd[6*ithread+5] * Ynm * eim;
}

__global__ void BiotSavartP2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ gpureal targetShrd[3];
  __shared__ gpureal sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[7*isource+5];
      __syncthreads();
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = sourceShrd[6*i+0] - targetShrd[0];
        d.y = sourceShrd[6*i+1] - targetShrd[1];
        d.z = sourceShrd[6*i+2] - targetShrd[2];
        gpureal rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        BiotSavartP2M_core(target,rho,alpha,beta,sourceShrd,i);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    if( threadIdx.x < size - iblok * THREADS ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[7*isource+5];
    }
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = sourceShrd[6*i+0] - targetShrd[0];
      d.y = sourceShrd[6*i+1] - targetShrd[1];
      d.z = sourceShrd[6*i+2] - targetShrd[2];
      gpureal rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      BiotSavartP2M_core(target,rho,alpha,beta,sourceShrd,i);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

template<>
void Kernel<BiotSavart>::P2M() {
  cudaThreadSynchronize();
  startTimer("P2M GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartP2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("P2M GPUkernel");\
}

__device__ void BiotSavartM2M_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
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
        target[0] += sourceShrd[6*jnkms+0] * CnmReal;
        target[0] -= sourceShrd[6*jnkms+1] * CnmImag;
        target[1] += sourceShrd[6*jnkms+0] * CnmImag;
        target[1] += sourceShrd[6*jnkms+1] * CnmReal;
        target[2] += sourceShrd[6*jnkms+2] * CnmReal;
        target[2] -= sourceShrd[6*jnkms+3] * CnmImag;
        target[3] += sourceShrd[6*jnkms+2] * CnmImag;
        target[3] += sourceShrd[6*jnkms+3] * CnmReal;
        target[4] += sourceShrd[6*jnkms+4] * CnmReal;
        target[4] -= sourceShrd[6*jnkms+5] * CnmImag;
        target[5] += sourceShrd[6*jnkms+4] * CnmImag;
        target[5] += sourceShrd[6*jnkms+5] * CnmReal;
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
        target[0] += sourceShrd[6*jnkms+0] * CnmReal;
        target[0] += sourceShrd[6*jnkms+1] * CnmImag;
        target[1] += sourceShrd[6*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[6*jnkms+1] * CnmReal;
        target[2] += sourceShrd[6*jnkms+2] * CnmReal;
        target[2] += sourceShrd[6*jnkms+3] * CnmImag;
        target[3] += sourceShrd[6*jnkms+2] * CnmImag;
        target[3] -= sourceShrd[6*jnkms+3] * CnmReal;
        target[4] += sourceShrd[6*jnkms+4] * CnmReal;
        target[4] += sourceShrd[6*jnkms+5] * CnmImag;
        target[5] += sourceShrd[6*jnkms+4] * CnmImag;
        target[5] -= sourceShrd[6*jnkms+5] * CnmReal;
      }
    }
  }
}

__global__ void BiotSavartM2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ gpureal sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    gpureal rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    BiotSavartM2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

template<>
void Kernel<BiotSavart>::M2M() {
  cudaThreadSynchronize();
  startTimer("M2M GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartM2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("M2M GPUkernel");\
}

__device__ void BiotSavartM2L_core(gpureal *target, gpureal  beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
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
      target[0] += sourceShrd[6*i+0] * CnmReal;
      target[0] += sourceShrd[6*i+1] * CnmImag;
      target[1] += sourceShrd[6*i+0] * CnmImag;
      target[1] -= sourceShrd[6*i+1] * CnmReal;
      target[2] += sourceShrd[6*i+2] * CnmReal;
      target[2] += sourceShrd[6*i+3] * CnmImag;
      target[3] += sourceShrd[6*i+2] * CnmImag;
      target[3] -= sourceShrd[6*i+3] * CnmReal;
      target[4] += sourceShrd[6*i+4] * CnmReal;
      target[4] += sourceShrd[6*i+5] * CnmImag;
      target[5] += sourceShrd[6*i+4] * CnmImag;
      target[5] -= sourceShrd[6*i+5] * CnmReal;
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
      target[0] += sourceShrd[6*i+0] * CnmReal;
      target[0] -= sourceShrd[6*i+1] * CnmImag;
      target[1] += sourceShrd[6*i+0] * CnmImag;
      target[1] += sourceShrd[6*i+1] * CnmReal;
      target[2] += sourceShrd[6*i+2] * CnmReal;
      target[2] -= sourceShrd[6*i+3] * CnmImag;
      target[3] += sourceShrd[6*i+2] * CnmImag;
      target[3] += sourceShrd[6*i+3] * CnmReal;
      target[4] += sourceShrd[6*i+4] * CnmReal;
      target[4] -= sourceShrd[6*i+5] * CnmImag;
      target[5] += sourceShrd[6*i+4] * CnmImag;
      target[5] += sourceShrd[6*i+5] * CnmReal;
    }
  }
}

__global__ void BiotSavartM2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ gpureal sourceShrd[6*THREADS];
  __shared__ gpureal factShrd[2*P];
  __shared__ gpureal YnmShrd[4*NTERM];
  gpureal fact = 1e-6;
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
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
            BiotSavartM2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
          }
        }
      }
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0] * 1e-6;
  targetGlob[6*itarget+1] = target[1] * 1e-6;
  targetGlob[6*itarget+2] = target[2] * 1e-6;
  targetGlob[6*itarget+3] = target[3] * 1e-6;
  targetGlob[6*itarget+4] = target[4] * 1e-6;
  targetGlob[6*itarget+5] = target[5] * 1e-6;
}

template<>
void Kernel<BiotSavart>::M2L() {
  cudaThreadSynchronize();
  startTimer("M2L GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartM2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("M2L GPUkernel");\
}

__device__ void BiotSavartM2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd) {
  gpureal x = cosf(theta);
  gpureal y = sinf(theta);
  if( fabsf(y) < EPS ) y = 1 / EPS;
  gpureal s = sqrtf(1 - x * x);
  gpureal spherical[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  gpureal cartesian[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  gpureal fact = 1;
  gpureal pn = 1;
  gpureal rhom = 1 / r;
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
    gpureal realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
    gpureal imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
    spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
    imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
    spherical[3] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[4] += 2 * YnmTheta * realj;
    spherical[5] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
    imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
    spherical[6] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[7] += 2 * YnmTheta * realj;
    spherical[8] -= 2 * m * Ynm * imagj;
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
      realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
      imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
      spherical[0] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
      imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
      spherical[3] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[4] += 2 * YnmTheta * realj;
      spherical[5] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
      imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
      spherical[6] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[7] += 2 * YnmTheta * realj;
      spherical[8] -= 2 * m * Ynm * imagj;
      rhon /= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,&spherical[0],&cartesian[0]);
  sph2cart(r,theta,phi,&spherical[3],&cartesian[3]);
  sph2cart(r,theta,phi,&spherical[6],&cartesian[6]);
  target[0] += 0.25 / M_PI * (cartesian[5] - cartesian[7]);
  target[1] += 0.25 / M_PI * (cartesian[6] - cartesian[2]);
  target[2] += 0.25 / M_PI * (cartesian[1] - cartesian[3]);
}

__global__ void BiotSavartM2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal targetX[3];
  gpureal target[3] = {0, 0, 0};
  __shared__ gpureal sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
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
            BiotSavartM2P_core(target,r,theta,phi,factShrd,sourceShrd);
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

template<>
void Kernel<BiotSavart>::M2P() {
  cudaThreadSynchronize();
  startTimer("M2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartM2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("M2P GPUkernel");\
}

__device__ inline void BiotSavartP2P_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[7*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[7*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[7*i+2];
#if 0
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z + EPS2;
  gpureal RS = R2 * sourceShrd[7*i+6];
  gpureal invR = rsqrtf(R2);
  if( R2 == 0 ) invR = 0;
  gpureal invR3 = invR * invR * invR;
  gpureal cutoff = 0.25 / M_PI * invR3 * (erff( sqrtf(RS) )
               - sqrtf(4 / M_PI * RS) * expf(-RS));
  target[0] += (d.y * sourceShrd[7*i+5] - d.z * sourceShrd[7*i+4]) * cutoff;
  target[1] += (d.z * sourceShrd[7*i+3] - d.x * sourceShrd[7*i+5]) * cutoff;
  target[2] += (d.x * sourceShrd[7*i+4] - d.y * sourceShrd[7*i+3]) * cutoff;
#else
  const gpureal SQRT4PI = M_2_SQRTPI;
  const gpureal FOURPI = 0.25 * M_1_PI;
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z + EPS2;
  gpureal invR = rsqrtf(R2);
  if( R2 == 0 ) invR = 0;
  gpureal invR3 = invR * invR * invR;
  gpureal RS = R2 * sourceShrd[7*i+6];
  gpureal SQRT_RS = sqrtf(RS);
  gpureal z = SQRT_RS,t,r;
  (t)=1.0f/(1.0f+0.5f*(z));
  (r)=(t)*expf(-(z)*(z)-1.26551223f+(t)*(1.00002368f+(t)*(0.37409196f+(t)*(0.09678418f+
      (t)*(-0.18628806f+(t)*(0.27886807f+(t)*(-1.13520398f+(t)*(1.48851587f+
      (t)*(-0.82215223f+(t)*0.17087277f)))))))));
  gpureal cutoff = FOURPI * invR3 * ( 1.0f - r
               - SQRT4PI * SQRT_RS * expf(-RS));
  target[0] += (d.y * sourceShrd[7*i+5] - d.z * sourceShrd[7*i+4]) * cutoff;
  target[1] += (d.z * sourceShrd[7*i+3] - d.x * sourceShrd[7*i+5]) * cutoff;
  target[2] += (d.x * sourceShrd[7*i+4] - d.y * sourceShrd[7*i+3]) * cutoff;
#endif
}

__global__ void BiotSavartP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal targetX[3];
  gpureal target[3] = {0, 0, 0};
  __shared__ gpureal sourceShrd[7*THREADS];
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
      sourceShrd[7*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[7*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[7*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[7*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[7*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[7*threadIdx.x+5] = sourceGlob[7*isource+5];
      sourceShrd[7*threadIdx.x+6] = 0.5f / (sourceGlob[7*isource+6] * sourceGlob[7*isource+6]);
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
                BiotSavartP2P_core(target,targetX,sourceShrd,d,i);
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
      sourceShrd[7*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[7*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[7*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[7*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[7*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[7*threadIdx.x+5] = sourceGlob[7*isource+5];
      sourceShrd[7*threadIdx.x+6] = 0.5f / (sourceGlob[7*isource+6] * sourceGlob[7*isource+6]);
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
              BiotSavartP2P_core(target,targetX,sourceShrd,d,i);
            }
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

template<>
void Kernel<BiotSavart>::P2P() {
  cudaThreadSynchronize();
  startTimer("P2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartP2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("P2P GPUkernel");\
}

__device__ void BiotSavartL2L_core(gpureal *target, gpureal beta, gpureal *factShrd, gpureal *YnmShrd, gpureal *sourceShrd) {
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
      target[0] += sourceShrd[6*nms+0] * CnmReal;
      target[0] += sourceShrd[6*nms+1] * CnmImag;
      target[1] += sourceShrd[6*nms+0] * CnmImag;
      target[1] -= sourceShrd[6*nms+1] * CnmReal;
      target[2] += sourceShrd[6*nms+2] * CnmReal;
      target[2] += sourceShrd[6*nms+3] * CnmImag;
      target[3] += sourceShrd[6*nms+2] * CnmImag;
      target[3] -= sourceShrd[6*nms+3] * CnmReal;
      target[4] += sourceShrd[6*nms+4] * CnmReal;
      target[4] += sourceShrd[6*nms+5] * CnmImag;
      target[5] += sourceShrd[6*nms+4] * CnmImag;
      target[5] -= sourceShrd[6*nms+5] * CnmReal;
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
        target[0] += sourceShrd[6*nms+0] * CnmReal;
        target[0] -= sourceShrd[6*nms+1] * CnmImag;
        target[1] += sourceShrd[6*nms+0] * CnmImag;
        target[1] += sourceShrd[6*nms+1] * CnmReal;
        target[2] += sourceShrd[6*nms+2] * CnmReal;
        target[2] -= sourceShrd[6*nms+3] * CnmImag;
        target[3] += sourceShrd[6*nms+2] * CnmImag;
        target[3] += sourceShrd[6*nms+3] * CnmReal;
        target[4] += sourceShrd[6*nms+4] * CnmReal;
        target[4] -= sourceShrd[6*nms+5] * CnmImag;
        target[5] += sourceShrd[6*nms+4] * CnmImag;
        target[5] += sourceShrd[6*nms+5] * CnmReal;
      }
    }
  }
}

__global__ void BiotSavartL2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ gpureal sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    gpureal rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    BiotSavartL2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

template<>
void Kernel<BiotSavart>::L2L() {
  cudaThreadSynchronize();
  startTimer("L2L GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartL2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("L2L GPUkernel");\
}

__device__ void BiotSavartL2P_core(gpureal *target, gpureal r, gpureal theta, gpureal phi, gpureal *factShrd, gpureal *sourceShrd) {
  gpureal x = cosf(theta);
  gpureal y = sinf(theta);
  if( fabsf(y) < EPS ) y = 1 / EPS;
  gpureal s = sqrtf(1 - x * x);
  gpureal spherical[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  gpureal cartesian[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    gpureal realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
    gpureal imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
    spherical[0] += 2 * m / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
    imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
    spherical[3] += 2 * m / r * Ynm * realj;
    spherical[4] += 2 * YnmTheta * realj;
    spherical[5] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
    imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
    spherical[6] += 2 * m / r * Ynm * realj;
    spherical[7] += 2 * YnmTheta * realj;
    spherical[8] -= 2 * m * Ynm * imagj;
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
      realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
      imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
      spherical[0] += 2 * n / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
      imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
      spherical[3] += 2 * n / r * Ynm * realj;
      spherical[4] += 2 * YnmTheta * realj;
      spherical[5] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
      imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
      spherical[6] += 2 * n / r * Ynm * realj;
      spherical[7] += 2 * YnmTheta * realj;
      spherical[8] -= 2 * m * Ynm * imagj;
      rhon *= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,&spherical[0],&cartesian[0]);
  sph2cart(r,theta,phi,&spherical[3],&cartesian[3]);
  sph2cart(r,theta,phi,&spherical[6],&cartesian[6]);
  target[0] += 0.25 / M_PI * (cartesian[5] - cartesian[7]);
  target[1] += 0.25 / M_PI * (cartesian[6] - cartesian[2]);
  target[2] += 0.25 / M_PI * (cartesian[1] - cartesian[3]);
}

__global__ void BiotSavartL2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal targetX[3];
  gpureal target[3] = {0, 0, 0};
  __shared__ gpureal sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    gpureal r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    BiotSavartL2P_core(target,r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

template<>
void Kernel<BiotSavart>::L2P() {
  cudaThreadSynchronize();
  startTimer("L2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    BiotSavartL2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("L2P GPUkernel");\
}
