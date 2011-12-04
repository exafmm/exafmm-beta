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
#include "gaussian.h"
#include "gpu.h"

template<>
void Kernel<Gaussian>::initialize() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(DEVICE);                                        // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

template<>
void Kernel<Gaussian>::M2M_CPU() {}

template<>
void Kernel<Gaussian>::finalize() {}

template<>
void Kernel<Gaussian>::allocate() {
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
void Kernel<Gaussian>::hostToDevice() {
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
void Kernel<Gaussian>::deviceToHost() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

template<>
void Kernel<Gaussian>::P2M() {}

template<>
void Kernel<Gaussian>::M2M() {}

template<>
void Kernel<Gaussian>::M2L() {}

template<>
void Kernel<Gaussian>::M2P() {}

__device__ inline void GaussianP2P_core(gpureal &target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[5*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[5*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[5*i+2];
  gpureal S2 = 2 * sourceShrd[5*i+4] * sourceShrd[5*i+4];
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z;
  target += sourceShrd[5*i+3] / (M_PI * S2) * rsqrtf(M_PI * S2) * expf(-R2 / S2);
}

__global__ void GaussianP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal targetX[3];
  gpureal target = 0;
  __shared__ gpureal sourceShrd[5*THREADS];
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
      sourceShrd[5*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[5*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[5*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[5*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[5*threadIdx.x+4] = sourceGlob[7*isource+6];
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
                GaussianP2P_core(target,targetX,sourceShrd,d,i);
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
      sourceShrd[5*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[5*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[5*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[5*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[5*threadIdx.x+4] = sourceGlob[7*isource+6];
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
              GaussianP2P_core(target,targetX,sourceShrd,d,i);
            }
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target;
}

template<>
void Kernel<Gaussian>::P2P() {
  cudaThreadSynchronize();
  startTimer("P2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    GaussianP2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("P2P GPUkernel");\
}

template<>
void Kernel<Gaussian>::L2L() {}

__global__ void GaussianL2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = 0;
}

template<>
void Kernel<Gaussian>::L2P() {
  cudaThreadSynchronize();
  startTimer("L2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    GaussianL2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("L2P GPUkernel");\
}
