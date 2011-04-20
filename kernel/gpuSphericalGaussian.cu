#include "kernel.h"
#include "gaussian.h"
#include "pregpu.h"

void Kernel::GaussianInit() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(DEVICE);                                        // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

__global__ void GaussianP2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {}

__global__ void GaussianM2M_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {}

void Kernel::GaussianM2M_CPU() {}

__global__ void GaussianM2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {}

__global__ void GaussianM2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {}

__device__ inline void GaussianP2P_core(gpureal &target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[5*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[5*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[5*i+2];
  gpureal S2 = 2 * sourceShrd[5*i+4] * sourceShrd[5*i+4];
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z + EPS2;
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

__global__ void GaussianL2L_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {}

__global__ void GaussianL2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = 0;
}

void Kernel::GaussianFinal() {}

#include "gpu.h"

CALL_GPU(GaussianP2M,P2M GPUkernel);
CALL_GPU(GaussianM2M,M2M GPUkernel);
CALL_GPU(GaussianM2L,M2L GPUkernel);
CALL_GPU(GaussianM2P,M2P GPUkernel);
CALL_GPU(GaussianP2P,P2P GPUkernel);
CALL_GPU(GaussianL2L,L2L GPUkernel);
CALL_GPU(GaussianL2P,L2P GPUkernel);
