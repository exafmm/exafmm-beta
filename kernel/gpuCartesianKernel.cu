#include "kernel.h"
#include "gpu.h"

void Kernel::initialize() {
//  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
  cudaSetDevice(1);
  cudaThreadSynchronize();                                      // Sync GPU threads
}

void Kernel::P2M() {
  for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {
    vect dist = CJ->X - B->pos;
    CJ->M[0] += B->scal;
    CJ->M[1] += B->scal * dist[0];
    CJ->M[2] += B->scal * dist[1];
    CJ->M[3] += B->scal * dist[2];
    CJ->M[4] += B->scal * dist[0] * dist[0] / 2;
    CJ->M[5] += B->scal * dist[1] * dist[1] / 2;
    CJ->M[6] += B->scal * dist[2] * dist[2] / 2;
    CJ->M[7] += B->scal * dist[0] * dist[1];
    CJ->M[8] += B->scal * dist[1] * dist[2];
    CJ->M[9] += B->scal * dist[2] * dist[0];
  }
}

void Kernel::M2M() {
  vect dist = CI->X - CJ->X;
  CI->M[0] += CJ->M[0];
  CI->M[1] += CJ->M[1] +  dist[0] * CJ->M[0];
  CI->M[2] += CJ->M[2] +  dist[1] * CJ->M[0];
  CI->M[3] += CJ->M[3] +  dist[2] * CJ->M[0];
  CI->M[4] += CJ->M[4] +  dist[0] * CJ->M[1] + dist[0] * dist[0]  * CJ->M[0] / 2;
  CI->M[5] += CJ->M[5] +  dist[1] * CJ->M[2] + dist[1] * dist[1]  * CJ->M[0] / 2;
  CI->M[6] += CJ->M[6] +  dist[2] * CJ->M[3] + dist[2] * dist[2]  * CJ->M[0] / 2;
  CI->M[7] += CJ->M[7] + (dist[0] * CJ->M[2] + dist[1] * CJ->M[1] + dist[0] * dist[1] * CJ->M[0]) / 2;
  CI->M[8] += CJ->M[8] + (dist[1] * CJ->M[3] + dist[2] * CJ->M[2] + dist[1] * dist[2] * CJ->M[0]) / 2;
  CI->M[9] += CJ->M[9] + (dist[2] * CJ->M[1] + dist[0] * CJ->M[3] + dist[2] * dist[0] * CJ->M[0]) / 2;
}

void Kernel::M2L() {
  vect dist = CI->X - CJ->X;
  real R = std::sqrt(norm(dist));
  real R3 = R * R * R;
  real R5 = R3 * R * R;
  CI->L[0] += CJ->M[0] / R;
  CI->L[0] += CJ->M[1] * (-dist[0] / R3);
  CI->L[0] += CJ->M[2] * (-dist[1] / R3);
  CI->L[0] += CJ->M[3] * (-dist[2] / R3);
  CI->L[0] += CJ->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[0] += CJ->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[0] += CJ->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[0] += CJ->M[7] * (3 * dist[0] * dist[1] / R5);
  CI->L[0] += CJ->M[8] * (3 * dist[1] * dist[2] / R5);
  CI->L[0] += CJ->M[9] * (3 * dist[2] * dist[0] / R5);
  CI->L[1] += CJ->M[0] * (-dist[0] / R3);
  CI->L[1] += CJ->M[1] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[1] += CJ->M[2] * (3 * dist[0] * dist[1] / R5);
  CI->L[1] += CJ->M[3] * (3 * dist[0] * dist[2] / R5);
  CI->L[2] += CJ->M[0] * (-dist[1] / R3);
  CI->L[2] += CJ->M[1] * (3 * dist[1] * dist[0] / R5);
  CI->L[2] += CJ->M[2] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[2] += CJ->M[3] * (3 * dist[1] * dist[2] / R5);
  CI->L[3] += CJ->M[0] * (-dist[2] / R3);
  CI->L[3] += CJ->M[1] * (3 * dist[2] * dist[0] / R5);
  CI->L[3] += CJ->M[2] * (3 * dist[2] * dist[1] / R5);
  CI->L[3] += CJ->M[3] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[4] += CJ->M[0] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[5] += CJ->M[0] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[6] += CJ->M[0] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[7] += CJ->M[0] * (3 * dist[0] * dist[1] / R5);
  CI->L[8] += CJ->M[0] * (3 * dist[1] * dist[2] / R5);
  CI->L[9] += CJ->M[0] * (3 * dist[2] * dist[0] / R5);
}

void Kernel::M2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CJ->X;
    real R = std::sqrt(norm(dist));
    real R3 = R * R * R;
    real R5 = R3 * R * R;
    B->pot += CJ->M[0] / R;
    B->pot += CJ->M[1] * (-dist[0] / R3);
    B->pot += CJ->M[2] * (-dist[1] / R3);
    B->pot += CJ->M[3] * (-dist[2] / R3);
    B->pot += CJ->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    B->pot += CJ->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    B->pot += CJ->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    B->pot += CJ->M[7] * (3 * dist[0] * dist[1] / R5);
    B->pot += CJ->M[8] * (3 * dist[1] * dist[2] / R5);
    B->pot += CJ->M[9] * (3 * dist[2] * dist[0] / R5);
  }
}

__global__ void P2P_GPU(int *keysGlob, int *rangeGlob, float4 *targetGlob, float4 *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float3 d;
  __shared__ float4 sourceShrd[THREADS];
  float4 target = targetGlob[blockIdx.x * THREADS + threadIdx.x];
  for( int ilist=0; ilist!=numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    int size  = rangeGlob[keys+2*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      __syncthreads();
      sourceShrd[threadIdx.x] = sourceGlob[begin + iblok * THREADS + threadIdx.x];
      __syncthreads();
#pragma unroll 32
      for( int i=0; i<THREADS; i++ ) {
        d.x = target.x - sourceShrd[i].x;
        d.y = target.y - sourceShrd[i].y;
        d.z = target.z - sourceShrd[i].z;
        target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
      }
    }
    int iblok = (size-1)/THREADS;
    __syncthreads();
    sourceShrd[threadIdx.x] = sourceGlob[begin + iblok * THREADS + threadIdx.x];
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; i++ ) {
      d.x = target.x - sourceShrd[i].x;
      d.y = target.y - sourceShrd[i].y;
      d.z = target.z - sourceShrd[i].z;
      target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
    }
  }
  targetGlob[blockIdx.x * THREADS + threadIdx.x] = target;
}

void Kernel::P2P() {
#if 0                                                           // Emulated GPU kernel
  for( int blockIdx=0; blockIdx!=numBlocks; ++blockIdx ) {
    int keys = keysHost[blockIdx];
    int numList = rangeHost[keys];
    float3 d;
    float4 sourceShrd[THREADS];
    float4 target;
    for( int threadIdx=0; threadIdx!=THREADS; ++threadIdx ) {
      target.x = targetHost[4*(blockIdx * THREADS + threadIdx)+0];
      target.y = targetHost[4*(blockIdx * THREADS + threadIdx)+1];
      target.z = targetHost[4*(blockIdx * THREADS + threadIdx)+2];
      target.w = targetHost[4*(blockIdx * THREADS + threadIdx)+3];
      for( int ilist=0; ilist!=numList; ++ilist ) {
        int begin = rangeHost[keys+2*ilist+1];
        int size  = rangeHost[keys+2*ilist+2];
        for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
          for( int i=0; i<THREADS; i++ ) {
            sourceShrd[i].x = sourceHost[4*(begin + iblok * THREADS + i)+0];
            sourceShrd[i].y = sourceHost[4*(begin + iblok * THREADS + i)+1];
            sourceShrd[i].z = sourceHost[4*(begin + iblok * THREADS + i)+2];
            sourceShrd[i].w = sourceHost[4*(begin + iblok * THREADS + i)+3];
          }
          for( int i=0; i<THREADS; i++ ) {
            d.x = target.x - sourceShrd[i].x;
            d.y = target.y - sourceShrd[i].y;
            d.z = target.z - sourceShrd[i].z;
            target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
          }
        }
        int iblok = (size-1)/THREADS;
        for( int i=0; i<THREADS; i++ ) {
          sourceShrd[i].x = sourceHost[4*(begin + iblok * THREADS + i)+0];
          sourceShrd[i].y = sourceHost[4*(begin + iblok * THREADS + i)+1];
          sourceShrd[i].z = sourceHost[4*(begin + iblok * THREADS + i)+2];
          sourceShrd[i].w = sourceHost[4*(begin + iblok * THREADS + i)+3];
        }
        for( int i=0; i<size-iblok*THREADS; i++ ) {
          d.x = target.x - sourceShrd[i].x;
          d.y = target.y - sourceShrd[i].y;
          d.z = target.z - sourceShrd[i].z;
          target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
        }
      }
      targetHost[4*(blockIdx * THREADS + threadIdx)+0] = target.x;
      targetHost[4*(blockIdx * THREADS + threadIdx)+1] = target.y;
      targetHost[4*(blockIdx * THREADS + threadIdx)+2] = target.z;
      targetHost[4*(blockIdx * THREADS + threadIdx)+3] = target.w;
    }
  }
#else
  int numBlocks = keysHost.size();
  int numRanges = rangeHost.size();
  int numTarget = targetHost.size() / 4;
  int numSource = sourceHost.size() / 4;
  cudaMalloc( (void**) &keysDevc,   numBlocks*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  numRanges*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, numTarget*sizeof(float4) );
  cudaMalloc( (void**) &sourceDevc, numSource*sizeof(float4) );
  cudaMemcpy(keysDevc,  &keysHost[0],  numBlocks*sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], numRanges*sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],numTarget*sizeof(float4),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],numSource*sizeof(float4),cudaMemcpyHostToDevice);
  P2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaMemcpy(&targetHost[0],targetDevc,numTarget*sizeof(float4),cudaMemcpyDeviceToHost);
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
#endif
}

void Kernel::L2L() {
  vect dist = CI->X - CJ->X;
  for( int i=0; i<10; ++i )
    CI->L[i] += CJ->L[i];
  CI->L[0] += CJ->L[1] * dist[0];
  CI->L[0] += CJ->L[2] * dist[1];
  CI->L[0] += CJ->L[3] * dist[2];
  CI->L[0] += CJ->L[4] * dist[0] * dist[0] / 2;
  CI->L[0] += CJ->L[5] * dist[1] * dist[1] / 2;
  CI->L[0] += CJ->L[6] * dist[2] * dist[2] / 2;
  CI->L[0] += CJ->L[7] * dist[0] * dist[1];
  CI->L[0] += CJ->L[8] * dist[1] * dist[2];
  CI->L[0] += CJ->L[9] * dist[2] * dist[0];
  CI->L[1] += CJ->L[4] * dist[0] * dist[0] / 2;
  CI->L[1] += CJ->L[7] * dist[0] * dist[1];
  CI->L[1] += CJ->L[9] * dist[0] * dist[2];
  CI->L[2] += CJ->L[7] * dist[1] * dist[0];
  CI->L[2] += CJ->L[5] * dist[1] * dist[1] / 2;
  CI->L[2] += CJ->L[8] * dist[1] * dist[2];
  CI->L[3] += CJ->L[9] * dist[2] * dist[0];
  CI->L[3] += CJ->L[8] * dist[2] * dist[1];
  CI->L[3] += CJ->L[6] * dist[2] * dist[2] / 2;
}

void Kernel::L2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CI->X;
    B->pot += CI->L[0];
    B->pot += CI->L[1] * dist[0];
    B->pot += CI->L[2] * dist[1];
    B->pot += CI->L[3] * dist[2];
    B->pot += CI->L[4] * dist[0] * dist[0] / 2;
    B->pot += CI->L[5] * dist[1] * dist[1] / 2;
    B->pot += CI->L[6] * dist[2] * dist[2] / 2;
    B->pot += CI->L[7] * dist[0] * dist[1];
    B->pot += CI->L[8] * dist[1] * dist[2];
    B->pot += CI->L[9] * dist[2] * dist[0];
  }
}

void Kernel::finalize() {}
