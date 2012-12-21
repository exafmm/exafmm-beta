#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include "matrixMul.h"

#define CHECK_BANK_CONFLICTS 0
#if CHECK_BANK_CONFLICTS
#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
#define DS(i, j) CUT_BANK_CHECKER(((float*)&Ds[0][0]), (BLOCK_SIZE * i + j))
#define ES(i, j) CUT_BANK_CHECKER(((float*)&Es[0][0]), (BLOCK_SIZE * i + j))
#define FS(i, j) CUT_BANK_CHECKER(((float*)&Fs[0][0]), (BLOCK_SIZE * i + j))
#else
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#define DS(i, j) Ds[i][j]
#define ES(i, j) Es[i][j]
#define FS(i, j) Fs[i][j]
#endif

extern "C"
__global__ void
gpuvdwforce_kernel( float* C, float* A, float* B, int* D, int num, float xmax, int nat,
		    float r2min, float r2max)
{
  int by = blockIdx.x;
  int ty = threadIdx.x;
  int aBegin = THD * by * 3 + ty*3;
  int jEnd = (num+THD-1) / THD;
  //  int jEnd = num / THD;
  int ty3, i3, jty3;
  int itype;

  float l = xmax;
  float al = 1.e0 / xmax;
  //  float exclude_radius2 = 0.00001e0;
  //  float cutoff_radius2 = 64.e0;
  float dn2, dn6;
  float dx, dy, dz, tmp, r2inv;

  float Csub[3];
  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;

  ty3 = ty * 3;

  __shared__ float Fs[1][THD];
  FS(0, ty)   = B[ty];

  __shared__ float As[1][THD * 3];
  AS(0, ty3)   = A[aBegin];
  AS(0, ty3+1) = A[aBegin+1];
  AS(0, ty3+2) = A[aBegin+2];
  __shared__ int Ds[1][THD];
  //  DS(0, ty)   = D[THD * by + ty] - 1;
  DS(0, ty)   = D[THD * by + ty];

  for (int j = 0; j < jEnd; j++){

    jty3 = j * THD * 3 + ty*3;

    __shared__ float Bs[1][THD * 3];
    BS(0,ty3)   = A[jty3];
    BS(0,ty3+1) = A[jty3+1];
    BS(0,ty3+2) = A[jty3+2];
    __shared__ int Es[1][THD];
    //    ES(0, ty)   = D[THD * j + ty] - 1;
    ES(0, ty)   = D[THD * j + ty];

    __syncthreads();
    for (int i = 0; i < THD; i++){

      i3 = i * 3;

      dx = AS(0,ty3)   - BS(0,i3);
      dy = AS(0,ty3+1) - BS(0,i3+1);
      dz = AS(0,ty3+2) - BS(0,i3+2);

      dx = dx - rintf(dx * al) * l;
      dy = dy - rintf(dy * al) * l;
      dz = dz - rintf(dz * al) * l;

      itype = nat * DS(0,ty) + ES(0,i);

      //dn2 = (dx * dx + dy * dy + dz * dz) * 1.e0;
      dn2 = (dx * dx + dy * dy + dz * dz) * FS(0,itype*2+1);
      //      if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
      //      if (dn2 > exclude_radius2){
      if (dn2 >= r2min && dn2<r2max && dn2!=0.0f){
	r2inv=1.0f/dn2;
	dn6 = r2inv*r2inv*r2inv;
	//Csub += FS(0,itype*2) * dn6 * (dn6 - 1.e0);
	tmp = FS(0,itype*2) * dn6 * r2inv * (2.e0 * dn6 - 1.e0);
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
      }
    }
    __syncthreads();
  }

  C[aBegin]   = Csub[0];
  C[aBegin+1] = Csub[1];
  C[aBegin+2] = Csub[2];
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
