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
#include <cutil.h>
__device__ __constant__ gpureal constDevc[514];                 // Constants on device

template<>
void Kernel<VanDerWaals>::initialize() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(DEVICE);                                        // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

template<>
void Kernel<VanDerWaals>::finalize() {
}

template<>
void Kernel<VanDerWaals>::allocate() {
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
void Kernel<VanDerWaals>::hostToDevice() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  constHost.push_back(ATOMS);
  for( int i=0; i!=ATOMS*ATOMS; ++i ) {
    constHost.push_back(RSCALE[i]);
    constHost.push_back(GSCALE[i]);
  }
  assert( constHost.size() == 514 );
  CUDA_SAFE_CALL(cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(gpureal)));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

template<>
void Kernel<VanDerWaals>::deviceToHost() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

__device__ inline void VanDerWaalsP2P_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd,
                                          float3 d, int i, int atoms) {
  d.x += targetX[0];
  d.x -= sourceShrd[4*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[4*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[4*i+2];
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z;
  if( R2 == 0 ) return;
  int atypei = targetX[3];
  int atypej = sourceShrd[4*i+3];
  gpureal rscale = constDevc[2*(atypei*atoms+atypej)+2];
  gpureal gscale = constDevc[2*(atypei*atoms+atypej)+3];
  gpureal R2s = R2 * rscale;
  if( R2MIN > R2s || R2s >= R2MAX ) return;
  gpureal invR = rsqrtf(R2s);
  gpureal invR2 = invR * invR;
  gpureal invR6 = invR2 * invR2 * invR2;
  gpureal dtmp = gscale * invR6 * invR2 * (2.0f * invR6 - 1.0f);
  target[0] += gscale * invR6 * (invR6 - 1.0f);
  target[1] -= d.x * dtmp;
  target[2] -= d.y * dtmp;
  target[3] -= d.z * dtmp;
}

__global__ void VanDerWaalsP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  int atoms = constDevc[1];
  gpureal targetX[4];
  gpureal target[4] = {0, 0, 0, 0};
  __shared__ gpureal sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[4*itarget+0];
  targetX[1] = targetGlob[4*itarget+1];
  targetX[2] = targetGlob[4*itarget+2];
  targetX[3] = targetGlob[4*itarget+3];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int size      = rangeGlob[keys+3*ilist+2];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
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
                VanDerWaalsP2P_core(target,targetX,sourceShrd,d,i,atoms);
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
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
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
              VanDerWaalsP2P_core(target,targetX,sourceShrd,d,i,atoms);
            }
          }
        }
      }
    }
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

template<>
void Kernel<VanDerWaals>::P2P() {
  cudaThreadSynchronize();
  startTimer("P2P GPUkernel");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    VanDerWaalsP2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("P2P GPUkernel");\
}

template<>
void Kernel<VanDerWaals>::P2M() {}

template<>
void Kernel<VanDerWaals>::M2M() {}

template<>
void Kernel<VanDerWaals>::M2L() {}

template<>
void Kernel<VanDerWaals>::M2P() {}

template<>
void Kernel<VanDerWaals>::L2L() {}

template<>
void Kernel<VanDerWaals>::L2P() {}
