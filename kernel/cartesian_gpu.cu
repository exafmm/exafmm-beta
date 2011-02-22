#include "kernel.h"
#include "gpu.h"

__global__ void P2P_GPU(float4 *sourceGlob, float *targetGlob) {
  int N = deviceConstant[0];
  float3 d;
  __shared__ float4 sourceShrd[THREADS];
  float4 target = sourceGlob[blockIdx.x * THREADS + threadIdx.x];
  target.w *= -rsqrtf(EPS2);
  for( int iblok=0; iblok<(N-1)/THREADS; iblok++) {
    __syncthreads();
    sourceShrd[threadIdx.x] = sourceGlob[iblok * THREADS + threadIdx.x];
    __syncthreads();
#pragma unroll 32
    for( int i=0; i<THREADS; i++ ) {
      d.x = target.x - sourceShrd[i].x;
      d.y = target.y - sourceShrd[i].y;
      d.z = target.z - sourceShrd[i].z;
      target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
    }
  }
  int iblok = (N-1)/THREADS;
  __syncthreads();
  sourceShrd[threadIdx.x] = sourceGlob[iblok * THREADS + threadIdx.x];
  __syncthreads();
  for( int i=0; i<N - (iblok * THREADS); i++ ) {
    d.x = target.x - sourceShrd[i].x;
    d.y = target.y - sourceShrd[i].y;
    d.z = target.z - sourceShrd[i].z;
    target.w += sourceShrd[i].w * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
  }
  targetGlob[blockIdx.x * THREADS + threadIdx.x] = target.w;
}

void Kernel::setup() {
  cudaSetDevice(MPIRANK % GPUS);
  cudaThreadSynchronize();
}

void Kernel::P2P(B_iter B0, B_iter BN) {
  hostConstant[0] = BN - B0;
  int Nround = ((BN - B0) / THREADS + 1) * THREADS;
  float4 *sourceHost, *sourceDevc;
  float  *targetHost, *targetDevc;

  sourceHost = (float4*)     malloc( Nround*sizeof(float4) );
  targetHost = (float *)     malloc( Nround*sizeof(float ) );
  cudaMalloc(  (void**) &sourceDevc, Nround*sizeof(float4) );
  cudaMalloc(  (void**) &targetDevc, Nround*sizeof(float ) );

  for( B_iter B=B0; B!=BN; ++B ) {
    sourceHost[B-B0].x = B->pos[0];
    sourceHost[B-B0].y = B->pos[1];
    sourceHost[B-B0].z = B->pos[2];
    sourceHost[B-B0].w = B->scal;
  }

  cudaMemcpyToSymbol(deviceConstant,hostConstant,4*sizeof(int));
  cudaMemcpy(sourceDevc,sourceHost,Nround*sizeof(float4),cudaMemcpyHostToDevice);
  P2P_GPU<<< Nround/THREADS, THREADS >>>(sourceDevc,targetDevc);
  cudaMemcpy(targetHost,targetDevc,Nround*sizeof(float ),cudaMemcpyDeviceToHost);

  for( B_iter B=B0; B!=BN; ++B ) {
    B->pot = targetHost[B-B0];
  }
}

void Kernel::P2P(B_iter Bi0, B_iter BiN, B_iter Bj0, B_iter BjN) {
  for( B_iter Bi=Bi0; Bi!=BiN; ++Bi ) {
    for( B_iter Bj=Bj0; Bj!=BjN; ++Bj ) {
      vect dist = Bi->pos - Bj->pos;
      real r = std::sqrt(norm(dist) + EPS2);
      Bi->pot += Bj->scal / r;
    }
  }
}

void Kernel::P2M(C_iter C) {
  C->M = 0;
  for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {
    vect dist = C->X - B->pos;
    C->M[0] += B->scal;
    C->M[1] += B->scal * dist[0];
    C->M[2] += B->scal * dist[1];
    C->M[3] += B->scal * dist[2];
    C->M[4] += B->scal * dist[0] * dist[0] / 2;
    C->M[5] += B->scal * dist[1] * dist[1] / 2;
    C->M[6] += B->scal * dist[2] * dist[2] / 2;
    C->M[7] += B->scal * dist[0] * dist[1];
    C->M[8] += B->scal * dist[1] * dist[2];
    C->M[9] += B->scal * dist[2] * dist[0];
  }
}

void Kernel::M2M(C_iter CI, C_iter CJ) {
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

void Kernel::M2L(C_iter CI, C_iter CJ) {
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

void Kernel::L2L(C_iter CI, C_iter CJ) {
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

void Kernel::L2P(C_iter C) {
  for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {
    vect dist = B->pos - C->X;
    B->pot += C->L[0];
    B->pot += C->L[1] * dist[0];
    B->pot += C->L[2] * dist[1];
    B->pot += C->L[3] * dist[2];
    B->pot += C->L[4] * dist[0] * dist[0] / 2;
    B->pot += C->L[5] * dist[1] * dist[1] / 2;
    B->pot += C->L[6] * dist[2] * dist[2] / 2;
    B->pot += C->L[7] * dist[0] * dist[1];
    B->pot += C->L[8] * dist[1] * dist[2];
    B->pot += C->L[9] * dist[2] * dist[0];
  }
}

void Kernel::M2P(C_iter CI, C_iter CJ) {
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
