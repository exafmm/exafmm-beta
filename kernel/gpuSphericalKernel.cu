#include "kernel.h"
#include "gpu.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;
const int  P4 = P2 * P2;
const real EPS = 1e-6;
double *prefactor, *Anm;
complex *Ynm, *Cnm, I(0.0,1.0);

__device__ void cart2sph(float& r, float& theta, float& phi, float dx, float dy, float dz) {
  r = sqrtf(dx * dx + dy * dy + dz * dz)+EPS;
  theta = acosf(dz / r);
  if( fabs(dx) + fabs(dy) < EPS ) {
    phi = 0;
  } else if( fabs(dx) < EPS ) {
    phi = dy / fabs(dy) * M_PI * 0.5;
  } else if( dx > 0 ) {
    phi = atanf(dy / dx);
  } else {
    phi = atanf(dy / dx) + M_PI;
  }
}

__device__ void M2L_Ynm(float* YnmShrd, float rho, float alpha, float* factShrd) {
  float x = cosf(alpha);
  float s = sqrt(1 - x * x);
  float fact = 1;
  float pn = 1;
  float rhom = 1.0 / rho;
  for( int m=0; m<2*P; ++m ){
    float p = pn;
    int i = m * (m + 1) /2 + m;
    YnmShrd[i] = rhom * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom /= rho;
    float rhon = rhom;
    for( int n=m+1; n<2*P; ++n ){
      i = n * (n + 1) / 2 + m;
      YnmShrd[i] = rhon * p * factShrd[n-m];
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon /= rho;
    }
    pn = -pn * fact * s;
    fact = fact + 2;
  }
}

__device__ void M2L_core(float* target, float  beta, float* factShrd, float* YnmShrd, float* sourceShrd) {
  int j = floor(sqrt(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=-n; m<0; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = anm * ajk * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 - m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] += sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] -= sourceShrd[2*i+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = ODDEVEN((abs(k - m) - k - m) / 2);
      cnm *= anm * ajk * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 + m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] -= sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] += sourceShrd[2*i+1] * CnmReal;
    }
  }
}

__global__ void M2L_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[4*NCOEF];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    int size  = rangeGlob[keys+2*ilist+2];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NCOEF ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho, alpha, beta;
    cart2sph(rho, alpha, beta, d.x, d.y, d.z);
    M2L_Ynm(YnmShrd, rho, alpha, factShrd);
    M2L_core(target, beta, factShrd, YnmShrd, sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__global__ void P2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[4];
  __shared__ float sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  target[0] = targetGlob[4*itarget+0];
  target[1] = targetGlob[4*itarget+1];
  target[2] = targetGlob[4*itarget+2];
  target[3] = targetGlob[4*itarget+3];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+2*ilist+1];
    int size  = rangeGlob[keys+2*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
      __syncthreads();
#pragma unroll 64
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = target[0] - sourceShrd[4*i+0];
        d.y = target[1] - sourceShrd[4*i+1];
        d.z = target[2] - sourceShrd[4*i+2];
        target[3] += sourceShrd[4*i+3] * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
      }
    }
    int iblok = (size-1)/THREADS;
    __syncthreads();
    int isource = begin + iblok * THREADS + threadIdx.x;
    sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
    sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
    sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
    sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = target[0] - sourceShrd[4*i+0];
      d.y = target[1] - sourceShrd[4*i+1];
      d.z = target[2] - sourceShrd[4*i+2];
      target[3] += sourceShrd[4*i+3] * rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
    }
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

void Kernel::initialize() {
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
  cudaThreadSynchronize();                                      // Sync GPU threads
  prefactor = new double  [4*P2];
  Anm       = new double  [4*P2];
  Ynm       = new complex [4*P2];
  Cnm       = new complex [P4];

  for( int n=0; n!=2*P; ++n ) {
    for( int m=-n; m<=n; ++m ) {
      int nm = n*n+n+m;
      int nabsm = abs(m);
      double fnmm = 1.0;
      for( int i=1; i<=n-m; ++i ) fnmm *= i;
      double fnpm = 1.0;
      for( int i=1; i<=n+m; ++i ) fnpm *= i;
      double fnma = 1.0;
      for( int i=1; i<=n-nabsm; ++i ) fnma *= i;
      double fnpa = 1.0;
      for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;
      prefactor[nm] = std::sqrt(fnma/fnpa);
      Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);
    }
  }

  for( int j=0, jk=0, jknm=0; j!=P; ++j ) {
    for( int k=-j; k<=j; ++k, ++jk ){
      for( int n=0, nm=0; n!=P; ++n ) {
        for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {
          const int jnkm = (j+n)*(j+n)+j+n+m-k;
          Cnm[jknm] = std::pow(I,double(abs(k-m)-abs(k)-abs(m)))*(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]);
        }
      }
    }
  }
}

void cart2sph(real& r, real& theta, real& phi, vect dist) {
  r = std::sqrt(norm(dist))+EPS;
  theta = std::acos(dist[2] / r);
  if( std::abs(dist[0]) + std::abs(dist[1]) < EPS ) {
    phi = 0;
  } else if( std::abs(dist[0]) < EPS ) {
    phi = dist[1] / std::abs(dist[1]) * M_PI * 0.5;
  } else if( dist[0] > 0 ) {
    phi = std::atan(dist[1] / dist[0]);
  } else {
    phi = std::atan(dist[1] / dist[0]) + M_PI;
  }
}

void evalMultipole(real rho, real alpha, real beta) {
  double x = std::cos(alpha);
  double s = std::sqrt(1 - x * x);
  double fact = 1;
  double pn = 1;
  double rhom = 1;
  for( int m=0; m!=P; ++m ){
    complex eim = std::exp(I * double(m * beta));
    double p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    Ynm[npn] = rhom * p * prefactor[npn] * eim;
    Ynm[nmn] = std::conj(Ynm[npn]);
    double p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho;
    double rhon = rhom;
    for( int n=m+1; n!=P; ++n ){
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      Ynm[npm] = rhon * p * prefactor[npm] * eim;
      Ynm[nmm] = std::conj(Ynm[npm]);
      double p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= rho;
    }
    pn = -pn * fact * s;
    fact = fact + 2;
  }
}

void evalLocal(real rho, real alpha, real beta) {
  double x = std::cos(alpha);
  double s = std::sqrt(1 - x * x);
  double fact = 1;
  double pn = 1;
  double rhom = 1.0 / rho;
  for( int m=0; m!=2*P; ++m ){
    complex eim = std::exp(I * double(m * beta));
    double p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    Ynm[npn] = rhom * p * prefactor[npn] * eim;
    Ynm[nmn] = std::conj(Ynm[npn]);
    double p1 = p;
    p = x * (2 * m + 1) * p;
    rhom /= rho;
    double rhon = rhom;
    for( int n=m+1; n!=2*P; ++n ){
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      Ynm[npm] = rhon * p * prefactor[npm] * eim;
      Ynm[nmm] = std::conj(Ynm[npm]);
      double p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon /= rho;
    }
    pn = -pn * fact * s;
    fact = fact + 2;
  }
}

void Kernel::P2M() {
  for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {
    vect dist = B->pos - CJ->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        const int nm  = n * n + n + m;
        const int nms = n * (n + 1) / 2 + m;
        CJ->M[nms] += double(B->scal)*Ynm[nm];
      }
    }
  }
}

void Kernel::M2M() {
  vect dist = CI->X - CJ->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,-beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex M = 0;
      for( int n=0; n<=j; ++n ) {
        for( int m=-n; m<=std::min(k-1,n); ++m ) {
          if( j-n >= k-m ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            const int nm    = n * n + n + m;
            M += CJ->M[jnkms]*std::pow(I,double(m-abs(m)))*Ynm[nm]*double(ODDEVEN(n)*Anm[nm]*Anm[jnkm]/Anm[jk]);
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            M += std::conj(CJ->M[jnkms])*Ynm[nm]*double(ODDEVEN(k+n+m)*Anm[nm]*Anm[jnkm]/Anm[jk]);
          }
        }
      }
      CI->M[jks] += M;
    }
  }
}

void Kernel::M2L() {
#if 0                                                           // Emulated GPU kernel
  int numBlocks = keysHost.size();
  for( int blockIdx=0; blockIdx!=numBlocks; ++blockIdx ) {
    int keys = keysHost[blockIdx];
    int numList = rangeHost[keys];
    float sourceShrd[2*THREADS];
    float factShrd[2*P];
    float YnmShrd[4*P*P];
    float fact = 1;
    for( int i=0; i<2*P; ++i ) {
      factShrd[i] = fact;
      fact *= i + 1;
    }
    for( int threadIdx=0; threadIdx!=THREADS; ++threadIdx ) {
      float target[2] = {0.0, 0.0};
      int itarget = blockIdx * THREADS;
      for( int ilist=0; ilist!=numList; ++ilist ) {
        int begin = rangeHost[keys+2*ilist+1];
        int size  = rangeHost[keys+2*ilist+2];
        vect dist;
        if(threadIdx==0) {
          dist[0] = targetHost[2*itarget+0] - sourceHost[begin+0];
          dist[1] = targetHost[2*itarget+1] - sourceHost[begin+1];
          dist[2] = targetHost[2*itarget+2] - sourceHost[begin+2];
        }
        for( int i=0; i<NCOEF; ++i ) {
          sourceShrd[2*i+0] = sourceHost[begin+2*i+3];
          sourceShrd[2*i+1] = sourceHost[begin+2*i+4];
        }
        float rho, alpha, beta;
        cart2sph(rho, alpha, beta, dist);
        float x = cosf(alpha);
        float s = sqrt(1 - x * x);
        float fact = 1;
        float pn = 1;
        float rhom = 1.0 / rho;
        for( int m=0; m<2*P; ++m ){
          float p = pn;
          int i = m * (m + 1) / 2 + m;
          YnmShrd[i] = rhom * p;
          float p1 = p;
          p = x * (2 * m + 1) * p;
          rhom /= rho;
          float rhon = rhom;
          for( int n=m+1; n<2*P; ++n ){
            i = n * (n + 1) / 2 + m;
            YnmShrd[i] = rhon * p * factShrd[n-m];
            float p2 = p1;
            p1 = p;
            p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
            rhon /= rho;
          }
          pn = -pn * fact * s;
          fact = fact + 2;
        }
        int j = floor(sqrt(2*threadIdx+0.25)-0.5);
        int k = 0;
        for( int i=0; i<=j; ++i ) k += i;
        k = threadIdx - k;
        if( threadIdx >= NCOEF ) j = k = 0;
        float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
        for( int n=0; n<P; ++n ) {
          for( int m=-n; m<0; ++m ) {
            int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
            float ere = cosf((m - k) * beta);
            float eim = sinf((m - k) * beta);
            float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
            float cnm = anm * ajk * YnmShrd[jnkm];
            float CnmReal = cnm * ere;
            float CnmImag = cnm * eim;
            int i = n * (n + 1) / 2 - m;
            target[0] += sourceShrd[2*i+0] * CnmReal;
            target[0] += sourceShrd[2*i+1] * CnmImag;
            target[1] += sourceShrd[2*i+0] * CnmImag;
            target[1] -= sourceShrd[2*i+1] * CnmReal;
          }
          for( int m=0; m<=n; ++m ) {
            int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
            float ere = cosf((m - k) * beta);
            float eim = sinf((m - k) * beta);
            float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
            float cnm = ODDEVEN((abs(k - m) - k - m) / 2);
            cnm *= anm * ajk * YnmShrd[jnkm];
            float CnmReal = cnm * ere;
            float CnmImag = cnm * eim;
            int i = n * (n + 1) / 2 + m;
            target[0] += sourceShrd[2*i+0] * CnmReal;
            target[0] -= sourceShrd[2*i+1] * CnmImag;
            target[1] += sourceShrd[2*i+0] * CnmImag;
            target[1] += sourceShrd[2*i+1] * CnmReal;
          }
        }
      }
      itarget = blockIdx * THREADS + threadIdx;
      targetHost[2*itarget+0] = target[0];
      targetHost[2*itarget+1] = target[1];
    }
  }
#else
  int numBlocks = keysHost.size();
  int numRanges = rangeHost.size();
  int numTarget = targetHost.size();
  int numSource = sourceHost.size();
  cudaMalloc( (void**) &keysDevc,   numBlocks*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  numRanges*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, numTarget*sizeof(float) );
  cudaMalloc( (void**) &sourceDevc, numSource*sizeof(float) );
  cudaMemcpy(keysDevc,  &keysHost[0],  numBlocks*sizeof(int),  cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], numRanges*sizeof(int),  cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],numTarget*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],numSource*sizeof(float),cudaMemcpyHostToDevice);
  M2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaMemcpy(&targetHost[0],targetDevc,numTarget*sizeof(float),cudaMemcpyDeviceToHost);
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
#endif
}

void Kernel::M2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CJ->X;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->pot += (CJ->M[nms]*Ynm[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->pot += 2*(CJ->M[nms]*Ynm[nm]).real();
      }
    }
  }
}

void Kernel::P2P() {
#if 0                                                           // Emulated GPU kernel
  int numBlocks = keysHost.size();
  for( int blockIdx=0; blockIdx!=numBlocks; ++blockIdx ) {
    int keys = keysHost[blockIdx];
    int numList = rangeHost[keys];
    vect dist;
    float sourceShrd[4*THREADS];
    float target[4];
    for( int threadIdx=0; threadIdx!=THREADS; ++threadIdx ) {
      int itarget = blockIdx * THREADS + threadIdx;
      target[0] = targetHost[4*itarget+0];
      target[1] = targetHost[4*itarget+1];
      target[2] = targetHost[4*itarget+2];
      target[3] = targetHost[4*itarget+3];
      for( int ilist=0; ilist!=numList; ++ilist ) {
        int begin = rangeHost[keys+2*ilist+1];
        int size  = rangeHost[keys+2*ilist+2];
        for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
          for( int i=0; i<THREADS; i++ ) {
            int isource = begin + iblok * THREADS + i;
            sourceShrd[4*i+0] = sourceHost[4*isource+0];
            sourceShrd[4*i+1] = sourceHost[4*isource+1];
            sourceShrd[4*i+2] = sourceHost[4*isource+2];
            sourceShrd[4*i+3] = sourceHost[4*isource+3];
          }
          for( int i=0; i<THREADS; i++ ) {
            dist[0] = target[0] - sourceShrd[4*i+0];
            dist[1] = target[1] - sourceShrd[4*i+1];
            dist[2] = target[2] - sourceShrd[4*i+2];
            target[3] += sourceShrd[4*i+3] * rsqrtf(norm(dist) + EPS2);
          }
        }
        int iblok = (size-1)/THREADS;
        for( int i=0; i<THREADS; i++ ) {
          int isource = begin + iblok * THREADS + i;
          sourceShrd[4*i+0] = sourceHost[4*isource+0];
          sourceShrd[4*i+1] = sourceHost[4*isource+1];
          sourceShrd[4*i+2] = sourceHost[4*isource+2];
          sourceShrd[4*i+3] = sourceHost[4*isource+3];
        }
        for( int i=0; i<size-iblok*THREADS; i++ ) {
          dist[0] = target[0] - sourceShrd[4*i+0];
          dist[1] = target[1] - sourceShrd[4*i+1];
          dist[2] = target[2] - sourceShrd[4*i+2];
          target[3] += sourceShrd[4*i+3] * rsqrtf(norm(dist) + EPS2);
        }
      }
      targetHost[4*itarget+0] = target[0];
      targetHost[4*itarget+1] = target[1];
      targetHost[4*itarget+2] = target[2];
      targetHost[4*itarget+3] = target[3];
    }
  }
#else
  int numBlocks = keysHost.size();
  int numRanges = rangeHost.size();
  int numTarget = targetHost.size();
  int numSource = sourceHost.size();
  cudaMalloc( (void**) &keysDevc,   numBlocks*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  numRanges*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, numTarget*sizeof(float) );
  cudaMalloc( (void**) &sourceDevc, numSource*sizeof(float) );
  cudaMemcpy(keysDevc,  &keysHost[0],  numBlocks*sizeof(int),  cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], numRanges*sizeof(int),  cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],numTarget*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],numSource*sizeof(float),cudaMemcpyHostToDevice);
  P2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaMemcpy(&targetHost[0],targetDevc,numTarget*sizeof(float),cudaMemcpyDeviceToHost);
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
#endif
}

void Kernel::L2L() {
  vect dist = CI->X - CJ->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex L = 0;
      for( int n=j; n!=P; ++n ) {
        for( int m=j+k-n; m<0; ++m ) {
          const int jnkm = (n - j) * (n - j) + n - j + m - k;
          const int nm   = n * n + n - m;
          const int nms  = n * (n + 1) / 2 - m;
          L += std::conj(CJ->L[nms])*Ynm[jnkm]*double(ODDEVEN(k)*Anm[jnkm]*Anm[jk]/Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
            L += CJ->L[nms]*std::pow(I,double(m-k-abs(m-k)))*Ynm[jnkm]*double(Anm[jnkm]*Anm[jk]/Anm[nm]);
          }
        }
      }
      CI->L[jks] += L;
    }
  }
}

void Kernel::L2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CI->X;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->pot += (CI->L[nms]*Ynm[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->pot += 2*(CI->L[nms]*Ynm[nm]).real();
      }
    }
  }
}

void Kernel::finalize() {
  delete[] prefactor;
  delete[] Anm;
  delete[] Ynm;
  delete[] Cnm;
}
