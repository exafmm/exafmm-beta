#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include "matrixMul.h"

#define CHECK_BANK_CONFLICTS 0

#if CHECK_BANK_CONFLICTS
#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
#define AS(i, j) CUT_BANK_CHECKER(((float*)&Cs[0][0]), (BLOCK_SIZE * i + j))
#else
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#define CS(i, j) Cs[i][j]
#endif

//#if CHECK_BANK_CONFLICTS
//#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
//#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
//#else
//#define AS(i, j) As[i][j]
//#define BS(i, j) Bs[i][j]
//#endif

extern "C"
__global__ void
gpuewwvpot_kernel( float* C, float* B, float* A, int num_atm)
{
  // A : k vertor
  // B : x and q

  int by = blockIdx.x;
  int ty = threadIdx.x;
  int ty4, i4, jty4, bty3;
  int bty4 = THD * by * 4 + ty*4;
  int jEnd = (num_atm+THD-1) / THD;

  float kr;
  float sin_theta, cos_theta;
  float qsin = 0.e0;
  float qcos = 0.e0;

  ty4 = ty * 4;

  __shared__ float As[1][THD * 4];
  AS(0, ty4)   = A[bty4];
  AS(0, ty4+1) = A[bty4+1];
  AS(0, ty4+2) = A[bty4+2];
  AS(0, ty4+3) = A[bty4+3];

  for (int j = 0; j < jEnd; j++){

    jty4 = j * THD * 4 + ty4;

    __shared__ float Bs[1][THD * 4];
    BS(0,ty4)   = B[jty4];
    BS(0,ty4+1) = B[jty4+1];
    BS(0,ty4+2) = B[jty4+2];
    BS(0,ty4+3) = B[jty4+3];

    __syncthreads();
    for (int i = 0; i < THD; i++){

      i4 = i * 4;

      kr = AS(0,ty4)   * BS(0,i4) 
	 + AS(0,ty4+1) * BS(0,i4+1) 
	 + AS(0,ty4+2) * BS(0,i4+2);
      sin_theta = sin(kr);
      cos_theta = cos(kr);
      qsin += sin_theta * BS(0,i4+3);
      qcos += cos_theta * BS(0,i4+3);

    }
    __syncthreads();
  }

  bty3 = THD * by * 3 + ty * 3;

  C[bty3]   = AS(0,ty4+3) * qsin;
  C[bty3+1] = AS(0,ty4+3) * qcos;
  C[bty3+2] = AS(0,ty4+3) * ( qsin * qsin + qcos * qcos );
  //C[0]   = 1.e0;
  //C[bty3+1] = 1.e0;
  //C[bty3+2] = 1.e0;
}


extern "C"
__global__ void
gpuewwvpot_kernel2( float* C, float* B, float* A, float* D, int num_k)
{
  // A : x and q
  // B : k vector

  int by = blockIdx.x;
  int ty = threadIdx.x;
  int i4, i3, jty4, jty3;
  int jEnd = (num_k+THD-1) / THD;
  int ty4 = ty * 4;
  int ty3 = ty * 3;
  int bty4 = THD * by * 4 + ty * 4;
  int bty3 = THD * by * 3 + ty * 3;

  float kr, tmp;
  float fx = 0.e0;
  float fy = 0.e0;
  float fz = 0.e0;

  __shared__ float As[1][THD * 4];
  AS(0, ty4)   = A[bty4];
  AS(0, ty4+1) = A[bty4+1];
  AS(0, ty4+2) = A[bty4+2];
  AS(0, ty4+3) = A[bty4+3];

  for (int j = 0; j < jEnd; j++){

    jty4 = j * THD * 4 + ty4;

    __shared__ float Bs[1][THD * 4];
    BS(0,ty4)   = B[jty4];
    BS(0,ty4+1) = B[jty4+1];
    BS(0,ty4+2) = B[jty4+2];
    BS(0,ty4+3) = B[jty4+3];

    jty3 = j * THD * 3 + ty3;

    __shared__ float Cs[1][THD * 3];
    CS(0,ty3)   = C[jty3]; // qsin
    CS(0,ty3+1) = C[jty3+1]; // qcos
    CS(0,ty3+2) = C[jty3+2];

    __syncthreads();
    for (int i = 0; i < THD; i++){

      i4 = i * 4;
      i3 = i * 3;
      kr = AS(0,ty4)   * BS(0,i4)
	+  AS(0,ty4+1) * BS(0,i4+1) 
	+  AS(0,ty4+2) * BS(0,i4+2);
      //      tmp = CS(0,i3+1) * sin(kr) - CS(0,i3) * cos(kr);
      //      fx += tmp * BS(0,i4);
      //      fy += tmp * BS(0,i4+1);
      //      fz += tmp * BS(0,i4+2);
      tmp = CS(0,i3+1) * cos(kr) + CS(0,i3) * sin(kr);
      fx += tmp;
    }
    __syncthreads();
  }

  D[bty3]   = AS(0,ty4+3) * fx;
  D[bty3+1] = AS(0,ty4+3) * fy;
  D[bty3+2] = AS(0,ty4+3) * fz;
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
