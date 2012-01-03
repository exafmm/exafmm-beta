#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cutil.h>

#include "mr3.h"

#if gpu
__device__ __constant__ VG_MATRIX      d_matrix[ATYPE2];

__device__ __inline__ 
void coulombforce_inter(int xj[3], float qj, int xi[3], float fi[3], float rscale2f, float al2)
{
  int k;
  float dn2,inr,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= rscale2f;
  inr   = rsqrtf(dn2);
  dphir = qj * inr * inr * inr;
  if(dn2==0.0f) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir * dr[k];
}

__global__ 
void coulombforce_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
			 float rscale2f, float xmax, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) coulombforce_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) coulombforce_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__ 
void coulombpot_inter(int xj[3], float qj, int xi[3], float fi[3], float rscale2f, float al2)
{
  int k;
  float dn2,inr,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= rscale2f;
  inr   = rsqrtf(dn2);
  dphir = qj * inr;
  if(dn2==0.0f) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir;
}

__global__ 
void coulombpot_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
		       float rscale2f, float xmax, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) coulombpot_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) coulombpot_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__ 
void realforce_inter(int xj[3], float qj, int xi[3], float fi[3], float rscale2f, float al2,
		     float r2min, float r2max)
{
  int k;
  float dn2,r,inr,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= rscale2f;
  inr   = rsqrtf(dn2);
  r     = inr * dn2;
  dphir = qj * ((float)(M_2_SQRTPI)*expf(-dn2) + erfcf(r)*inr)*inr*inr;
  if(dn2<r2min || dn2>=r2max) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir * dr[k];
}

__global__ 
void realforce_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
		      float rscale2f, float xmax, float r2min, float r2max, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) realforce_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2,r2min,r2max);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) realforce_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2,r2min,r2max);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__ 
void realpot_inter(int xj[3], float qj, int xi[3], float fi[3], float rscale2f, float al2,
		   float r2min, float r2max)
{
  int k;
  float dn2,r,inr,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= rscale2f;
  inr   = rsqrtf(dn2);
  r     = inr * dn2;
  dphir = qj * erfcf(r) * inr;
  if(dn2<r2min || dn2>=r2max) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir;
}

__global__ 
void realpot_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
		    float rscale2f, float xmax, float r2min, float r2max, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) realpot_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2,r2min,r2max);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) realpot_inter(s_xj[js].r,s_xj[js].qatype.q,xi,fi,rscale2f,al2,r2min,r2max);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__
void vdwforce_inter(int xj[3], int xi[3], float fi[3], int t, float al2,
                    float r2min, float r2max)
{
  int k;
  float dn2,inr2,dn6,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= d_matrix[t].rscale;
  inr2  = 1.0f/dn2;
  dn6   = inr2*inr2*inr2;
  dphir = d_matrix[t].gscale * dn6 * inr2 * (2.0f * dn6 - 1.0f);
  if(dn2<r2min || dn2>=r2max) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir * dr[k];
}

__global__
void vdwforce_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
                     int nat, float xmax, float r2min, float r2max, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js,atypei;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  atypei = xivec[i].qatype.atype * nat;
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) vdwforce_inter(s_xj[js].r,xi,fi,atypei+s_xj[js].qatype.atype,al2,r2min,r2max);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) vdwforce_inter(s_xj[js].r,xi,fi,atypei+s_xj[js].qatype.atype,al2,r2min,r2max);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__
void vdwpot_inter(int xj[3], int xi[3], float fi[3], int t, float al2,
                  float r2min, float r2max)
{
  int k;
  float dn2,inr2,dn6,dr[3],dphir;

  dn2 = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] - xj[k];
    dr[k] *= al2;
    dn2   += dr[k] * dr[k];
  }
  dn2  *= d_matrix[t].rscale;
  inr2  = 1.0f/dn2;
  dn6   = inr2*inr2*inr2;
  dphir = d_matrix[t].gscale * dn6 * (dn6 - 1.0f);
  if(dn2<r2min || dn2>=r2max) dphir = 0.0f;
  for(k=0; k<3; k++) fi[k] += dphir;
}

__global__
void vdwpot_kernel(int ni, VG_XVEC *xivec, int nj, VG_XVEC *xjvec,
                   int nat, float xmax, float r2min, float r2max, float *fvec)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float fi[3],al2;
  int js,atypei;
  __shared__ VG_XVEC s_xj[NLOAD];
  int xi[3];

  al2=scalbnf(xmax,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = xivec[i].r[k];
  atypei = xivec[i].qatype.atype * nat;
  for (j = 0; j < nj - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = xjvec[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) vdwpot_inter(s_xj[js].r,xi,fi,atypei+s_xj[js].qatype.atype,al2,r2min,r2max);
  }
  __syncthreads();
  if(tid < nj - j) s_xj[tid] = xjvec[j + tid];
  __syncthreads();
  for (js = 0; js < nj - j; js++) vdwpot_inter(s_xj[js].r,xi,fi,atypei+s_xj[js].qatype.atype,al2,r2min,r2max);
  if(i<ni) for(k=0; k<3; k++) fvec[i*3+k] = fi[k];
}


__device__ __inline__ 
void ewald_dft_inter(int xj[3], float qj, float ki[3], float factor1,
		     float bsbci[3], float al2)
{
  int k;
  float th,dr[3],s,c;

  th = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xj[k] * al2;
    th    += dr[k] * ki[k];
  }
  th  *= (float)(2.0 * M_PI);
  s    = qj * sinf(th);
  c    = qj * cosf(th);
  bsbci[0] += s * factor1;
  bsbci[1] += c * factor1;
  //  bsbci[2] += (s * s + c * c) * 0.5f * factor1;
}

__global__ 
void ewald_dft_kernel(VG_XVEC *x, int n, VG_KVEC *kvec, int knum, float *bsbc)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float bsbci[3],al2;
  int js;
  __shared__ VG_XVEC s_xj[NLOAD];
  float ki[3],factor1;

  al2=scalbnf(1.0f,-32);
  for(k=0; k<3; k++) bsbci[k] = 0.0f;
  for(k=0; k<3; k++) ki[k] = kvec[i].k[k];
  factor1                  = kvec[i].factor1;
  for (j = 0; j < n - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_xj[tid] = x[j + tid];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) ewald_dft_inter(s_xj[js].r,s_xj[js].qatype.q,ki,factor1,bsbci,al2);
  }
  __syncthreads();
  if(tid < n - j) s_xj[tid] = x[j + tid];
  __syncthreads();
  for (js = 0; js < n - j; js++) ewald_dft_inter(s_xj[js].r,s_xj[js].qatype.q,ki,factor1,bsbci,al2);
  if(i<knum) for(k=0; k<3; k++) bsbc[i*3+k] = bsbci[k];
}


__device__ __inline__ 
void ewaldforce_idft_inter(float kj[3], float bsbc[3], int xi[3], float fi[3], float al2)
{
  int k;
  float th,dr[3],s,c;

  th = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] * al2;
    th    += dr[k] * kj[k];
  }
  th  *= (float)(2.0 * M_PI);
  s    = sinf(th);
  c    = cosf(th);
  for(k=0; k<3; k++){
    fi[k] += (bsbc[1] * s - bsbc[0] * c) * kj[k];
  }
}

__global__ 
void ewaldforce_idft_kernel(VG_XVEC *x, int n, VG_KVEC *kvec, int knum, 
			    float *bsbc, float *force)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float al2;
  int js,xi[3];
  __shared__ VG_KVEC s_kj[NLOAD];
  __shared__ float   s_bsbcj[NLOAD][3];
  float fi[3];

  al2=scalbnf(1.0f,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = x[i].r[k];
  for (j = 0; j < knum - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_kj[tid] = kvec[j + tid];
    if(tid < NLOAD) for(k=0;k<2;k++) s_bsbcj[tid][k] = bsbc[(j + tid)*3 + k];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) ewaldforce_idft_inter(s_kj[js].k,s_bsbcj[js],xi,fi,al2);
  }
  __syncthreads();
  if(tid < knum - j) s_kj[tid] = kvec[j + tid];
  if(tid < knum - j) for(k=0;k<2;k++) s_bsbcj[tid][k] = bsbc[(j + tid)*3 + k];
  __syncthreads();
  for (js = 0; js < knum - j; js++) ewaldforce_idft_inter(s_kj[js].k,s_bsbcj[js],xi,fi,al2);
  if(i<n) for(k=0; k<3; k++) force[i*3+k] = fi[k];
}


__device__ __inline__ 
void ewaldpot_idft_inter(float kj[3], float bsbc[3], int xi[3], float fi[3], float al2)
{
  int k;
  float th,dr[3],s,c;

  th = 0.0f;
  for(k=0; k<3; k++){
    dr[k]  = xi[k] * al2;
    th    += dr[k] * kj[k];
  }
  th  *= (float)(2.0 * M_PI);
  s    = sinf(th);
  c    = cosf(th);
  fi[0] += bsbc[1] * c + bsbc[0] * s;
}

__global__ 
void ewaldpot_idft_kernel(VG_XVEC *x, int n, VG_KVEC *kvec, int knum, 
			  float *bsbc, float *force)
{
  int tid = threadIdx.x;
  int i = blockIdx.x * NTHRE + tid;
  int j,k;
  float al2;
  int js,xi[3];
  __shared__ VG_KVEC s_kj[NLOAD];
  __shared__ float   s_bsbcj[NLOAD][3];
  float fi[3];

  al2=scalbnf(1.0f,-32);
  for(k=0; k<3; k++) fi[k] = 0.0f;
  for(k=0; k<3; k++) xi[k] = x[i].r[k];
  for (j = 0; j < knum - NLOAD; j+=NLOAD){
    __syncthreads();
    if(tid < NLOAD) s_kj[tid] = kvec[j + tid];
    if(tid < NLOAD) for(k=0;k<2;k++) s_bsbcj[tid][k] = bsbc[(j + tid)*3 + k];
    __syncthreads();
#pragma unroll 16
    for (js = 0; js < NLOAD; js++) ewaldpot_idft_inter(s_kj[js].k,s_bsbcj[js],xi,fi,al2);
  }
  __syncthreads();
  if(tid < knum - j) s_kj[tid] = kvec[j + tid];
  if(tid < knum - j) for(k=0;k<2;k++) s_bsbcj[tid][k] = bsbc[(j + tid)*3 + k];
  __syncthreads();
  for (js = 0; js < knum - j; js++) ewaldpot_idft_inter(s_kj[js].k,s_bsbcj[js],xi,fi,al2);
  if(i<n) for(k=0; k<3; k++) force[i*3+k] = fi[k];
}

static void malloc_x(VG_XVEC **d_x, VG_XVEC **xf, int nalloc, int nthre)
{
  static VG_XVEC *d_x_static=NULL,*xf_static=NULL;
  static int nalloc_bak=0;
  if(nalloc>nalloc_bak){
    if(nalloc<NMAX) nalloc=NMAX;
    CUDA_SAFE_CALL(cudaFree(d_x_static));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_x_static,sizeof(VG_XVEC)*(nalloc+nthre)));
    free(xf_static);
    if((xf_static=(VG_XVEC *)malloc(sizeof(VG_XVEC)*(nalloc+nthre)))==NULL){
      fprintf(stderr,"** error : can't malloc xf_static **\n");
      exit(1);
    }
    nalloc_bak=nalloc;
  }
  *d_x=d_x_static;
  *xf=xf_static;
}


static void malloc_x2(VG_XVEC **d_x, VG_XVEC **xf, int nalloc, int nthre)
{
  static VG_XVEC *d_x_static=NULL,*xf_static=NULL;
  static int nalloc_bak=0;
  if(nalloc>nalloc_bak){
    if(nalloc<NMAX) nalloc=NMAX;
    CUDA_SAFE_CALL(cudaFree(d_x_static));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_x_static,sizeof(VG_XVEC)*(nalloc+nthre)));
    free(xf_static);
    if((xf_static=(VG_XVEC *)malloc(sizeof(VG_XVEC)*(nalloc+nthre)))==NULL){
      fprintf(stderr,"** error : can't malloc xf_static **\n");
      exit(1);
    }
    nalloc_bak=nalloc;
  }
  *d_x=d_x_static;
  *xf=xf_static;
}


static void malloc_f(float **d_force, float **forcef, int nalloc)
{
  static float *d_force_static=NULL,*forcef_static=NULL;
  static int nalloc_bak=0;
  if(nalloc>nalloc_bak){
    if(nalloc<NMAX) nalloc=NMAX;
    CUDA_SAFE_CALL(cudaFree(d_force_static));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_force_static,sizeof(float)*nalloc*3));
    free(forcef_static);
    if((forcef_static=(float *)malloc(sizeof(float)*nalloc*3))==NULL){
      fprintf(stderr,"** error : can't malloc forcef_static **\n");
      exit(1);
    }
    bzero(forcef_static,sizeof(float)*nalloc*3);
    nalloc_bak=nalloc;
  }
  *d_force=d_force_static;
  *forcef=forcef_static;
}


static void malloc_k(VG_KVEC **d_k, VG_KVEC **kf, int kalloc, int nthre)
{
  static VG_KVEC *d_k_static=NULL,*kf_static=NULL;
  static int kalloc_bak=0;
  if(kalloc>kalloc_bak){
    if(kalloc<KMAX) kalloc=KMAX;
    CUDA_SAFE_CALL(cudaFree(d_k_static));
    CUDA_SAFE_CALL(cudaMalloc((void**)&d_k_static,sizeof(VG_KVEC)*(kalloc+nthre)));
    free(kf_static);
    if((kf_static=(VG_KVEC *)malloc(sizeof(VG_KVEC)*(kalloc+nthre)))==NULL){
      fprintf(stderr,"** error : can't malloc kf_static **\n");
      exit(1);
    }
    kalloc_bak=kalloc;
  }
  *d_k=d_k_static;
  *kf=kf_static;
}


static void malloc_and_make_index(int n, int nat, int atype[], int **index_ret)
{
  int i,at,na[ATYPE],offset[ATYPE],*index;
  if((index=(int *)malloc(sizeof(int)*n))==NULL){
    fprintf(stderr,"** error : can't malloc index **\n");
    exit(1);
  }
  for(at=0;at<nat;at++) na[at]=0;
  for(i=0;i<n;i++) na[atype[i]]++;
  offset[0]=0;
  for(at=1;at<nat;at++) offset[at]=offset[at-1]+na[at-1];
  for(at=0;at<nat;at++) na[at]=0;
  for(i=0;i<n;i++){
    at=atype[i];
    index[i]=offset[at]+na[at];// i:original, index[i]:new
    na[at]++;
  }
  *index_ret=index;
}


static void free_index(int *index)
{
  free(index);
}


static void send_matrix(int nat, double *gscale, double *rscale)
{
  int i,j;
  VG_MATRIX *matrix;

  // send matrix
  if((matrix=(VG_MATRIX *)malloc(sizeof(VG_MATRIX)*nat*nat))==NULL){
    fprintf(stderr,"** error : can't malloc matrix **\n");
    exit(1);
  }
  for(i=0;i<nat;i++){
    for(j=0;j<nat;j++){
      matrix[i*nat+j].gscale=(float)(gscale[i*nat+j]);
      matrix[i*nat+j].rscale=(float)(rscale[i*nat+j]);
    }
  }
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_matrix,matrix,sizeof(VG_MATRIX)*nat*nat));
  free(matrix);
}


static void copy_and_send_x(int n, double *x, int *atype, int *index,
                            double xmax_1[3], VG_XVEC *vec, VG_XVEC *d_x)
{
  int i,j;
  DI2 di2;

  for(i=0;i<n;i++){
    int idx=index[i];
    for(j=0;j<3;j++){
      di2.d=x[i*3+j]*xmax_1[j]+0x180000;
      vec[idx].r[j]=di2.fi2.fi0.i;
    }
    vec[idx].qatype.atype=atype[i];
  }
  CUDA_SAFE_CALL(cudaMemcpy(d_x,vec,sizeof(VG_XVEC)*n,cudaMemcpyHostToDevice));
}


static void copy_and_send_xq(int n, double *x, double *q,
                             double xmax_1[3], VG_XVEC *vec, VG_XVEC *d_x)
{
  int i,j;
  DI2 di2;

  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      di2.d=x[i*3+j]*xmax_1[j]+0x180000;
      vec[i].r[j]=di2.fi2.fi0.i;
    }
    vec[i].qatype.q=q[i];
  }
  CUDA_SAFE_CALL(cudaMemcpy(d_x,vec,sizeof(VG_XVEC)*n,cudaMemcpyHostToDevice));
}


static void copy_and_send_k(int knum, int *k, double alpha, double epsilon,
                            double cellsize_1[3],
                            VG_KVEC *vec, VG_KVEC *d_k)
{
  int i,j;
  double ktmp,r2,kvtmp,eps1=1.0/epsilon,vol1;
  double alpha4=1.0/(4.0*alpha*alpha);

  for(j=0,vol1=1.0;j<3;j++) vol1*=cellsize_1[j];
  for(i=0;i<knum;i++){
    for(j=0,r2=0.0;j<3;j++){
      ktmp=(double)k[i*3+j];
      vec[i].k[j]=(float)ktmp;
      kvtmp=2.0*M_PI*ktmp*cellsize_1[j];
      r2+=kvtmp*kvtmp;
    }
    vec[i].factor1=2.0*eps1*vol1*exp(-r2*alpha4)/r2;
  }
  CUDA_SAFE_CALL(cudaMemcpy(d_k,vec,sizeof(VG_KVEC)*knum,cudaMemcpyHostToDevice));
}


static void get_result_q(int n, float *d_force, float *forcef,
                         double *q, double rfac, double *force)
{
  // copy GPU result to host, and convert it to double
  int i,j;
  double factor;

  CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
  for(i=0;i<n;i++){
    factor=q[i]*rfac;
    for(j=0;j<3;j++) force[i*3+j]+=forcef[i*3+j]*factor;
  }
}


static void get_result_q3(int n, float *d_force, float *forcef,
                          double *q, double rfac[3], double *force)
{
  // copy GPU result to host, and convert it to double
  int i,j;
  double factor;

  CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
  for(i=0;i<n;i++){
    for(j=0;j<3;j++){
      factor=q[i]*rfac[j];
      force[i*3+j]+=forcef[i*3+j]*factor;
    }
  }
}


static void get_result(int n, float *d_force, float *forcef, double *force)
{
  // copy GPU result to host, and convert it to double
  int i,j;

  CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
  for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=forcef[i*3+j];
}


static void get_result_index(int n, int *index, float *d_force, float *forcef, double *force)
{
  // copy GPU result to host, and convert it to double
  int i,j;

  CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
  for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]+=forcef[index[i]*3+j];
}


void MR3calccoulomb_ij(int ni, double xi[], double qi[], double force[],
                       int nj, double xj[], double qj[],
                       double rscale,
                       int tblno, double xmax, int periodicflag)
{
  VG_XVEC *d_xi,*xif;
  VG_XVEC *d_xj,*xjf;
  float *d_force,*forcef;
  if((periodicflag & 1)==0) xmax*=2.0;
  double xmax_1[3]={1.0/xmax,1.0/xmax,1.0/xmax};
  float r2min=MD_REAL_R2MIN,r2max=MD_REAL_R2MAX;
  float rscale2f=(float)(rscale*rscale);

  malloc_x(&d_xi,&xif,ni,NTHRE);
  malloc_x2(&d_xj,&xjf,nj,NTHRE);
  malloc_f(&d_force,&forcef,ni);
  copy_and_send_xq(ni,xi,qi,xmax_1,xif,d_xi);
  copy_and_send_xq(nj,xj,qj,xmax_1,xjf,d_xj);
  switch(tblno){
  case 0:
    coulombforce_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                       rscale2f,(float)xmax,d_force);
    break;
  case 1:
    coulombpot_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                     rscale2f,(float)xmax,d_force);
    break;
  case 6:
    realforce_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                    rscale2f,(float)xmax,r2min,r2max,d_force);
    break;
  case 7:
    realpot_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                       rscale2f,(float)xmax,r2min,r2max,d_force);
    break;
  default:
    fprintf(stderr,"** error : not supported tblno = %d **\n",tblno);
    exit(1);
    break;
  }
  CUT_CHECK_ERROR("Kernel execution failed");
  if(tblno==0 || tblno==1) get_result(ni,d_force,forcef,force);
  else if(tblno==6)        get_result_q(ni,d_force,forcef,qi,rscale*rscale*rscale,force);
  else if(tblno==7)        get_result_q(ni,d_force,forcef,qi,rscale,force);
}


void MR3calcvdw_ij(int ni, double xi[], int atypei[], double force[],
                   int nj, double xj[], int atypej[],
                   int nat, double gscale[], double rscale[],
                   int tblno, double xmax, int periodicflag)
{
  VG_XVEC *d_xi,*xif;
  VG_XVEC *d_xj,*xjf;
  float *d_force,*forcef;
  int *indexi,*indexj;
  if((periodicflag & 1)==0) xmax*=2.0;
  double xmax_1[3]={1.0/xmax,1.0/xmax,1.0/xmax};
  float r2min=MD_LJ_R2MIN,r2max=MD_LJ_R2MAX;

  if(nat>ATYPE){
    fprintf(stderr,"** error : nat is too large **\n");
    exit(1);
  }
  malloc_x(&d_xi,&xif,ni,NTHRE);
  malloc_x2(&d_xj,&xjf,nj,NTHRE);
  malloc_f(&d_force,&forcef,ni);
  malloc_and_make_index(ni,nat,atypei,&indexi);
  malloc_and_make_index(nj,nat,atypej,&indexj);
  send_matrix(nat,gscale,rscale);
  copy_and_send_x(ni,xi,atypei,indexi,xmax_1,xif,d_xi);
  copy_and_send_x(nj,xj,atypej,indexj,xmax_1,xjf,d_xj);
  switch(tblno){
  case 2:
    vdwforce_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                   nat,(float)xmax,r2min,r2max,d_force);
    break;
  case 3:
    vdwpot_kernel<<< (ni+NTHRE-1)/NTHRE, NTHRE >>>(ni,d_xi,nj,d_xj,
                 nat,(float)xmax,r2min,r2max,d_force);
    break;
  default:
    fprintf(stderr,"** error : not supported tblno = %d **\n",tblno);
    exit(1);
    break;
  }
  CUT_CHECK_ERROR("Kernel execution failed");
  get_result_index(ni,indexi,d_force,forcef,force);
  free_index(indexi);
  free_index(indexj);
}


void MR3calcewald(int *k, int knum_org, double *x, int n, double *q,
                  double alpha, double epsilon, double cell[3][3],
                  double *force, double *tpot, double stress[3][3])
{
  VG_XVEC *d_x,*xf;
  float *d_force,*forcef;
  int knum;
  VG_KVEC *d_k,*kf;
  float *d_bsbc,*bsbcf;
  double cellsize_1[3]={1.0/cell[0][0],1.0/cell[1][1],1.0/cell[2][2]};

  if(knum_org>=0) knum=knum_org;
  else            knum=-knum_org;
  malloc_x(&d_x,&xf,n,NTHRE);
  malloc_f(&d_force,&forcef,n);
  malloc_k(&d_k,&kf,knum,NTHRE);
  malloc_f(&d_bsbc,&bsbcf,knum);
  copy_and_send_xq(n,x,q,cellsize_1,xf,d_x);
  copy_and_send_k(knum,k,alpha,epsilon,cellsize_1,kf,d_k);
  ewald_dft_kernel<<< (knum+NTHRE-1)/NTHRE, NTHRE >>>(d_x,n,d_k,knum,d_bsbc);
  CUDA_SAFE_CALL(cudaMemcpy(bsbcf,d_bsbc,sizeof(float)*knum*3,cudaMemcpyDeviceToHost));
  *tpot=0.0;
  for(int i=0;i<knum;i++){
    *tpot+=(bsbcf[i*3]*bsbcf[i*3]+bsbcf[i*3+1]*bsbcf[i*3+1])*0.5/kf[i].factor1;
  }
  if(knum_org>=0){
    ewaldforce_idft_kernel<<< (n+NTHRE-1)/NTHRE, NTHRE >>>(d_x,n,d_k,knum,d_bsbc,d_force);
    for(int j=0;j<3;j++) cellsize_1[j] *= 2.0 * M_PI;
    get_result_q3(n,d_force,forcef,q,cellsize_1,force);
  }
  else{
    ewaldpot_idft_kernel<<< (n+NTHRE-1)/NTHRE, NTHRE >>>(d_x,n,d_k,knum,d_bsbc,d_force);
    get_result_q(n,d_force,forcef,q,1.0,force);
  }
}

#endif
