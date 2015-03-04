#ifndef _MATRIXMUL_KERNEL_H_
#define _MATRIXMUL_KERNEL_H_

#include <stdio.h>
#include "matrixMul.h"

#define CHECK_BANK_CONFLICTS 1
#if CHECK_BANK_CONFLICTS
#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
#else
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#endif


#define VMP 1
#define UNROLL
#define NO_SHARED_FOR_I

#ifdef NO_SHARED_FOR_I
#define KERNEL_CORE_IJ \
	dx  = Asreg[0] - BS(0,i4);\
	dy  = Asreg[1] - BS(0,i4+1);\
	dz  = Asreg[2] - BS(0,i4+2);\
	ar2 = dx * dx + dy * dy + dz * dz + eps2;\
	sqdn=rsqrtf(ar2);\
	tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;\
        i4+=4;\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;
#else
#define KERNEL_CORE_IJ \
	dx  = AS(0,ty4)   - BS(0,i4);\
	dy  = AS(0,ty4+1) - BS(0,i4+1);\
	dz  = AS(0,ty4+2) - BS(0,i4+2);\
	ar2 = dx * dx + dy * dy + dz * dz + eps2;\
	sqdn=rsqrtf(ar2);\
	tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;\
        i4+=4;\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;
#endif

#define KERNEL_CORE_IJ_VMP2(vmpid) \
	dxv[vmpid]  = AS(0,ty4)   - BS(0,i4v[vmpid]);\
	dyv[vmpid]  = AS(0,ty4+1) - BS(0,i4v[vmpid]+1);\
	dzv[vmpid]  = AS(0,ty4+2) - BS(0,i4v[vmpid]+2);\
	ar2v[vmpid] = dxv[vmpid] * dxv[vmpid] + dyv[vmpid] * dyv[vmpid] + dzv[vmpid] * dzv[vmpid] + eps2;\
	sqdnv[vmpid]=rsqrtf(ar2v[vmpid]);\
	tmpv[vmpid]=BS(0,i4v[vmpid]+3)*sqdnv[vmpid]*sqdnv[vmpid]*sqdnv[vmpid];\
        i4v[vmpid]+=4*VMP;\
	Csub[0] += tmpv[vmpid] * dxv[vmpid];\
	Csub[1] += tmpv[vmpid] * dyv[vmpid];\
	Csub[2] += tmpv[vmpid] * dzv[vmpid];

#define KERNEL_CORE_IJ_VMP22(vmpid) \
	dxv[vmpid]  = AS(0,ty4)   - BS(0,i4v[vmpid]);\
	dyv[vmpid]  = AS(0,ty4+1) - BS(0,i4v[vmpid]+1);\
	dzv[vmpid]  = AS(0,ty4+2) - BS(0,i4v[vmpid]+2);\
	dxv[vmpid+1]  = AS(0,ty4)   - BS(0,i4v[vmpid+1]);\
	dyv[vmpid+1]  = AS(0,ty4+1) - BS(0,i4v[vmpid+1]+1);\
	dzv[vmpid+1]  = AS(0,ty4+2) - BS(0,i4v[vmpid+1]+2);\
        ar2v[vmpid]  = dxv[vmpid] * dxv[vmpid] + eps2;\
        ar2v[vmpid+1]  = dxv[vmpid+1] * dxv[vmpid+1] + eps2;\
        ar2v[vmpid] += dyv[vmpid] * dyv[vmpid];\
        ar2v[vmpid+1] += dyv[vmpid+1] * dyv[vmpid+1];\
        ar2v[vmpid] += dzv[vmpid] * dzv[vmpid];\
        ar2v[vmpid+1] += dzv[vmpid+1] * dzv[vmpid+1];\
        sqdnv[vmpid] = ar2v[vmpid];\
        sqdnv[vmpid+1] = ar2v[vmpid+1];\
	tmpv[vmpid]=BS(0,i4v[vmpid]+3)*sqdnv[vmpid]*sqdnv[vmpid]*sqdnv[vmpid];\
	tmpv[vmpid+1]=BS(0,i4v[vmpid+1]+3)*sqdnv[vmpid+1]*sqdnv[vmpid+1]*sqdnv[vmpid+1];\
        i4v[vmpid]+=4*VMP;\
        i4v[vmpid+1]+=4*VMP;\
	Csub[0] += tmpv[vmpid] * dxv[vmpid];\
	Csub[0] += tmpv[vmpid+1] * dxv[vmpid+1];\
	Csub[1] += tmpv[vmpid] * dyv[vmpid];\
	Csub[1] += tmpv[vmpid+1] * dyv[vmpid+1];\
	Csub[2] += tmpv[vmpid] * dzv[vmpid];\
	Csub[2] += tmpv[vmpid+1] * dzv[vmpid+1];

#if 0
	sqdnv[vmpid]=rsqrtf(ar2v[vmpid]);\
	sqdnv[vmpid+1]=rsqrtf(ar2v[vmpid+1]);\

#endif


extern "C"
__global__ void
coulombforce_kernel( float* C, float* A, int num, float xmax, float alpha,
		     float eps2)
{
  int by = blockIdx.x;
#ifdef HALF_THREADS_PER_BLOCK   
  int ty = threadIdx.x*2;
#else
  int ty = threadIdx.x;
#endif
  int jEnd = (num+THD-1) / THD;
  int ty4, i4, jty4, byt3;

  float alpha2 = alpha  * alpha;
  float l = xmax;
  float al = 1.e0 / xmax;
  float exclude_radius2 = 0.00001e0;
  //  float exclude_radius2 = 0.01e0;
  //float cutoff_radius2 = 9.8e0 * alpha2;
  float ar2, sqdn, tmp, r2;
  float dx, dy, dz;
#ifdef THREAD_CONTIGUOUS
  int aBegin = THD * by * 4 + ty;
  __shared__ float As[1][THD*4];
  __shared__ float Bs[1][THD*4];
#elif defined(HALF_THREADS_PER_BLOCK)
  int aBegin = THD * by * 4 + ty * 4;
  __shared__ float As[1][THD * 4];
  __shared__ float Bs[1][THD * 4];
#else
  int aBegin = THD * by * 4 + ty * 4;
  __shared__ float As[1][THD * 4];
  __shared__ float Bs[1][THD * 4];
#endif
#ifdef HALF_THREADS_PER_BLOCK   
  float dx0,dx1,dy0,dy1,dz0,dz1;
  float ar20,ar21,sqdn0,sqdn1,tmp0,tmp1;
  float Csub[6];
  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;
  Csub[3] = 0.e0;
  Csub[4] = 0.e0;
  Csub[5] = 0.e0;
#else
  float Csub[3];
  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;
#endif

  ty4 = ty * 4;

  if(eps2==0.0f){
#ifdef THREAD_CONTIGUOUS
    AS(0, ty4)       = A[aBegin];
    AS(0, ty4+THD)   = A[aBegin+THD];
    AS(0, ty4+THD*2) = A[aBegin+THD*2];
    AS(0, ty4+THD*3) = A[aBegin+THD*3];
#else
    AS(0, ty4)   = A[aBegin];
    AS(0, ty4+1) = A[aBegin+1];
    AS(0, ty4+2) = A[aBegin+2];
    AS(0, ty4+3) = A[aBegin+3];
#endif

    for (int j = 0; j < jEnd; j++){
#ifdef THREAD_CONTIGUOUS
      jty4 = j * THD * 4 + ty;
      BS(0,ty4)       = A[jty4];
      BS(0,ty4+THD)   = A[jty4+THD];
      BS(0,ty4+THD*2) = A[jty4+THD*2];
      BS(0,ty4+THD*3) = A[jty4+THD*3];
#else
      jty4 = j * THD * 4 + ty4;
      BS(0,ty4)   = A[jty4];
      BS(0,ty4+1) = A[jty4+1];
      BS(0,ty4+2) = A[jty4+2];
      BS(0,ty4+3) = A[jty4+3];
#endif
      __syncthreads();
      for (int i = 0; i < THD; i++){
	i4 = i * 4;
#ifdef THREAD_CONTIGUOUS
	dx = AS(0,ty4)       - BS(0,i4);
	dy = AS(0,ty4+THD)   - BS(0,i4+THD);
	dz = AS(0,ty4+THD*2) - BS(0,i4+THD*2);
#else
	dx = AS(0,ty4)   - BS(0,i4);
	dy = AS(0,ty4+1) - BS(0,i4+1);
	dz = AS(0,ty4+2) - BS(0,i4+2);
#endif
	dx = dx - rintf(dx * al) * l;
	dy = dy - rintf(dy * al) * l;
	dz = dz - rintf(dz * al) * l;
	ar2 = alpha2 * (dx * dx + dy * dy + dz * dz);
	//if (ar2 > exclude_radius2 && ar2 < cutoff_radius2){
	if(ar2!=0.0f){
	  sqdn=rsqrtf(ar2);
#ifdef THREAD_CONTIGUOUS
	  tmp=BS(0,i4+THD*3)*sqdn*sqdn*sqdn;
#else
	  tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;
#endif
	  Csub[0] += tmp * dx;
	  Csub[1] += tmp * dy;
	  Csub[2] += tmp * dz;
	}
      }
      __syncthreads();
    }
  }
  else{
#ifdef THREAD_CONTIGUOUS
    AS(0, ty4)       = A[aBegin];
    AS(0, ty4+THD)   = A[aBegin+THD];
    AS(0, ty4+THD*2) = A[aBegin+THD*2];
    AS(0, ty4+THD*3) = A[aBegin+THD*3];
#elif defined(HALF_THREADS_PER_BLOCK)
    AS(0, ty4)   = A[aBegin];
    AS(0, ty4+1) = A[aBegin+1];
    AS(0, ty4+2) = A[aBegin+2];
    AS(0, ty4+3) = A[aBegin+3];
    AS(0, ty4+4) = A[aBegin+4];
    AS(0, ty4+5) = A[aBegin+5];
    AS(0, ty4+6) = A[aBegin+6];
    AS(0, ty4+7) = A[aBegin+7];
#else
    AS(0, ty4)   = A[aBegin];
    AS(0, ty4+1) = A[aBegin+1];
    AS(0, ty4+2) = A[aBegin+2];
    AS(0, ty4+3) = A[aBegin+3];
#endif
    for (int j = 0; j < jEnd; j++){
#ifdef THREAD_CONTIGUOUS
      jty4 = j * THD * 4 + ty;
      BS(0,ty4)       = A[jty4];
      BS(0,ty4+THD)   = A[jty4+THD];
      BS(0,ty4+THD*2) = A[jty4+THD*2];
      BS(0,ty4+THD*3) = A[jty4+THD*3];
#elif defined(HALF_THREADS_PER_BLOCK)
      jty4 = j * THD * 4 + ty4;
      BS(0,ty4)   = A[jty4];
      BS(0,ty4+1) = A[jty4+1];
      BS(0,ty4+2) = A[jty4+2];
      BS(0,ty4+3) = A[jty4+3];
      BS(0,ty4+4) = A[jty4+4];
      BS(0,ty4+5) = A[jty4+5];
      BS(0,ty4+6) = A[jty4+6];
      BS(0,ty4+7) = A[jty4+7];
#else
      jty4 = j * THD * 4 + ty4;
      BS(0,ty4)   = A[jty4];
      BS(0,ty4+1) = A[jty4+1];
      BS(0,ty4+2) = A[jty4+2];
      BS(0,ty4+3) = A[jty4+3];
#endif
      __syncthreads();
      for (int i = 0; i < THD; i++){
	i4 = i * 4;
#ifdef THREAD_CONTIGUOUS
	dx = AS(0,ty4)       - BS(0,i4);
	dy = AS(0,ty4+THD)   - BS(0,i4+THD);
	dz = AS(0,ty4+THD*2) - BS(0,i4+THD*2);
#elif defined(HALF_THREADS_PER_BLOCK)
	dx0 = AS(0,ty4)   - BS(0,i4);
	dy0 = AS(0,ty4+1) - BS(0,i4+1);
	dz0 = AS(0,ty4+2) - BS(0,i4+2);
	dx1 = AS(0,ty4+4) - BS(0,i4);
	dy1 = AS(0,ty4+5) - BS(0,i4+1);
	dz1 = AS(0,ty4+6) - BS(0,i4+2);
#else
	dx = AS(0,ty4)   - BS(0,i4);
	dy = AS(0,ty4+1) - BS(0,i4+1);
	dz = AS(0,ty4+2) - BS(0,i4+2);
#endif
#ifdef HALF_THREADS_PER_BLOCK   
	ar20 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0 + eps2;
	ar21 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1 + eps2;
	sqdn0=rsqrtf(ar20);
	sqdn1=rsqrtf(ar21);
#else
	ar2 = dx * dx + dy * dy + dz * dz + eps2;
	sqdn=rsqrtf(ar2);
#endif
#ifdef THREAD_CONTIGUOUS
	tmp=BS(0,i4+THD*3)*sqdn*sqdn*sqdn;
#elif defined(HALF_THREADS_PER_BLOCK)
	tmp0=BS(0,i4+3)*sqdn0*sqdn0*sqdn0;
	tmp1=BS(0,i4+3)*sqdn1*sqdn1*sqdn1;
#else
	tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;
#endif
#ifdef HALF_THREADS_PER_BLOCK   
	Csub[0] += tmp0 * dx0;
	Csub[1] += tmp0 * dy0;
	Csub[2] += tmp0 * dz0;
	Csub[3] += tmp1 * dx1;
	Csub[4] += tmp1 * dy1;
	Csub[5] += tmp1 * dz1;
#else
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
#endif
      }
      __syncthreads();
    }
  }
  byt3 = THD * by * 3 + ty * 3;
//  C[byt3]   = Csub[0] * tmp;
//  C[byt3+1] = Csub[1] * tmp;
//  C[byt3+2] = Csub[2] * tmp;
  C[byt3]   = Csub[0];
  C[byt3+1] = Csub[1];
  C[byt3+2] = Csub[2];
#ifdef HALF_THREADS_PER_BLOCK   
  C[byt3+3] = Csub[3];
  C[byt3+4] = Csub[4];
  C[byt3+5] = Csub[5];
#endif
}


extern "C"
__global__ void
coulombforce_ij_kernel( float* C, float *A, float *B, int numj, float xmax, float alpha,
			float eps2)
{
  //#define USE_TWOSETS_J
#if VMP!=1
  float dxv[VMP],dyv[VMP],dzv[VMP],ar2v[VMP],sqdnv[VMP],tmpv[VMP];
  float Csubv[VMP][3];
  int i4v[VMP];
#endif
  int by = blockIdx.x;
  int ty = threadIdx.x;
  int jEnd = (numj+THD-1) / THD;
  int ty4, i4, jty4, byt3;

  float alpha2 = alpha  * alpha;
  float l = xmax;
  float al = 1.e0 / xmax;
  float exclude_radius2 = 0.00001e0;
  float ar2, sqdn, tmp, r2, ar2x, ar2y, ar2z;
  float dx, dy, dz;
  int aBegin = THD * by * 4 + ty * 4;
  __shared__ float As[1][THD * 4];
#ifdef NO_SHARED_FOR_I
  float Asreg[4];
#endif
#ifdef USE_TWOSETS_J
  __shared__ float Bs[2][THD * 4];
  int msel=0;
#else
  __shared__ float Bs[1][THD * 4];
#endif
  float Csub[3];
  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;

  ty4 = ty * 4;

  if(eps2==0.0f){
    AS(0, ty4)   = A[aBegin];
    AS(0, ty4+1) = A[aBegin+1];
    AS(0, ty4+2) = A[aBegin+2];
    AS(0, ty4+3) = A[aBegin+3];
    for (int j = 0; j < jEnd; j++){
      jty4 = j * THD * 4 + ty4;
      BS(0,ty4)   = B[jty4];
      BS(0,ty4+1) = B[jty4+1];
      BS(0,ty4+2) = B[jty4+2];
      BS(0,ty4+3) = B[jty4+3];
      __syncthreads();
      for (int i = 0; i < THD; i++){
	i4 = i * 4;
	dx = AS(0,ty4)   - BS(0,i4);
	dy = AS(0,ty4+1) - BS(0,i4+1);
	dz = AS(0,ty4+2) - BS(0,i4+2);
	dx = dx - rintf(dx * al) * l;
	dy = dy - rintf(dy * al) * l;
	dz = dz - rintf(dz * al) * l;
	ar2 = alpha2 * (dx * dx + dy * dy + dz * dz);
	if(ar2!=0.0f){
	  sqdn=rsqrtf(ar2);
	  tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;
	  Csub[0] += tmp * dx;
	  Csub[1] += tmp * dy;
	  Csub[2] += tmp * dz;
	}
      }
      __syncthreads();
    }
  }
  else{
#ifdef NO_SHARED_FOR_I
    Asreg[0] = A[aBegin];
    Asreg[1] = A[aBegin+1];
    Asreg[2] = A[aBegin+2];
    Asreg[3] = A[aBegin+3];
#else
    AS(0, ty4)   = A[aBegin];
    AS(0, ty4+1) = A[aBegin+1];
    AS(0, ty4+2) = A[aBegin+2];
    AS(0, ty4+3) = A[aBegin+3];
#endif
#ifdef USE_TWOSETS_J
    // set first j-particle data
    BS(msel,ty4)   = B[ty4];
    BS(msel,ty4+1) = B[ty4+1];
    BS(msel,ty4+2) = B[ty4+2];
    BS(msel,ty4+3) = B[ty4+3];
    __syncthreads();
    for (int j = 0; j < jEnd; j++){
      jty4 = (j+1) * THD * 4 + ty4;
      BS(1-msel,ty4)   = B[jty4];
      BS(1-msel,ty4+1) = B[jty4+1];
      BS(1-msel,ty4+2) = B[jty4+2];
      BS(1-msel,ty4+3) = B[jty4+3];
      for (int i = 0; i < THD; i++){
	i4 = i * 4;
	dx = AS(0,ty4)   - BS(msel,i4);
	dy = AS(0,ty4+1) - BS(msel,i4+1);
	dz = AS(0,ty4+2) - BS(msel,i4+2);
	ar2 = dx * dx + dy * dy + dz * dz + eps2;
	sqdn=rsqrtf(ar2);
	tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
      }
      msel=1-msel;
      __syncthreads();
    }
#else
    for (int j = 0; j < jEnd; j++){
      jty4 = j * THD * 4 + ty4;
      BS(0,ty4)   = B[jty4];
      BS(0,ty4+1) = B[jty4+1];
      BS(0,ty4+2) = B[jty4+2];
      BS(0,ty4+3) = B[jty4+3];
      __syncthreads();
#if THD==64 && defined(UNROLL)
      i4=0;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
#elif THD==128 && VMP==1 && defined(UNROLL) // unroling
      i4=0;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
#elif THD==256 && VMP==1 && defined(UNROLL) // unroling
      i4=0;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
      KERNEL_CORE_IJ;
#elif VMP==2 // not work 
      for(int i=0;i<VMP;i++) i4v[i]=i*4;
      for (int i = 0; i < THD/VMP; i++){
	//	KERNEL_CORE_IJ_VMP2(0);
	//	KERNEL_CORE_IJ_VMP2(1);
	KERNEL_CORE_IJ_VMP22(0);
      }
#else
      for (int i = 0; i < THD; i++){
	i4  = i * 4;
	dx  = AS(0,ty4)   - BS(0,i4);
	dy  = AS(0,ty4+1) - BS(0,i4+1);
	dz  = AS(0,ty4+2) - BS(0,i4+2);
	ar2 = dx * dx + dy * dy + dz * dz + eps2;
	sqdn=rsqrtf(ar2);
	tmp=BS(0,i4+3)*sqdn*sqdn*sqdn;
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
      }
#endif
      __syncthreads();
    }
#endif
  }
  byt3 = THD * by * 3 + ty * 3;
  C[byt3]   = Csub[0];
  C[byt3+1] = Csub[1];
  C[byt3+2] = Csub[2];
}

#endif // #ifndef _MATRIXMUL_KERNEL_H_
