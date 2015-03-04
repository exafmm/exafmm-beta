#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include "matrixMul.h"

#define CHECK_BANK_CONFLICTS 0
#if CHECK_BANK_CONFLICTS
#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
#else
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#endif

extern "C"
__global__ void
coulombpot_kernel( float* C, float* A, int num, float xmax, float alpha)
{
  int by = blockIdx.x;
  int ty = threadIdx.x;
  int aBegin = THD * by * 4 + ty * 4;
  int jEnd = (num+THD-1) / THD;
  int ty4, i4, jty4;

  float alpha2 = alpha  * alpha;
  float l = xmax;
  float al = 1.e0 / xmax;
  float exclude_radius2 = 0.00001e0;
  //float cutoff_radius2 = 9.8e0;
  float ar2, sqdn;
  float dx, dy, dz;

  float Csub = 0.e0;

  ty4 = ty * 4;

  __shared__ float As[1][THD * 4];
  AS(0, ty4)   = A[aBegin];
  AS(0, ty4+1) = A[aBegin+1];
  AS(0, ty4+2) = A[aBegin+2];
  AS(0, ty4+3) = A[aBegin+3];

  for (int j = 0; j < jEnd; j++){

    jty4 = j * THD * 4 + ty4;

    __shared__ float Bs[1][THD * 4];
    BS(0,ty4)   = A[jty4];
    BS(0,ty4+1) = A[jty4+1];
    BS(0,ty4+2) = A[jty4+2];
    BS(0,ty4+3) = A[jty4+3];

    __syncthreads();
    for (int i = 0; i < THD; i++){

      i4 = i * 4;

      dx = AS(0,ty4)   - BS(0,i4);
      dy = AS(0,ty4+1) - BS(0,i4+1);
      dz = AS(0,ty4+2) - BS(0,i4+2);

      dx = dx - rintf(dx * al) * l;
      dy = dy - rintf(dy * al) * l;
      dz = dz - rintf(dz * al) * l;

//      dn2 = dx * dx + dy * dy + dz * dz;
      ar2 = alpha2 * (dx * dx + dy * dy + dz * dz);
      //if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
      if (ar2 > exclude_radius2){
//	sqdn = sqrt(dn2) * alpha;
//	Csub += erfc(sqdn) / sqdn * BS(0, i4+3);
	sqdn = sqrt(ar2);
	Csub += BS(0, i4+3)/sqdn;
      }
    }
    __syncthreads();
  }

//  C[THD * by * 3 + ty * 3] = Csub * AS(0,ty4+3) * alpha;
  C[THD * by * 3 + ty * 3] = Csub;
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
