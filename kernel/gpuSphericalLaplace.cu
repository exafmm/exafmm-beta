#include "kernel.h"
#include "spherical.h"
#include "laplace.h"
#include "gpu.h"

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

__device__ void sph2cart(float r, float theta, float phi, float *spherical, float *cartesian) {
  cartesian[0] += sinf(theta) * cosf(phi) * spherical[0]
                + cosf(theta) * cosf(phi) / r * spherical[1]
                - sinf(phi) / r / sinf(theta) * spherical[2];
  cartesian[1] += sinf(theta) * sinf(phi) * spherical[0]
                + cosf(theta) * sinf(phi) / r * spherical[1]
                + cosf(phi) / r / sinf(theta) * spherical[2];
  cartesian[2] += cosf(theta) * spherical[0]
                - sinf(theta) / r * spherical[1];
}

__device__ void evalMultipole(float *YnmShrd, float rho, float alpha, float *factShrd) {
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
  float fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<P; ++m ){
    float p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    YnmShrd[npn] = rhom * p / factShrd[2*m];
    YnmShrd[nmn] = YnmShrd[npn];
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= -rho;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ){
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      YnmShrd[npm] = rhon * p / factShrd[n+m];
      YnmShrd[nmm] = YnmShrd[npm];
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= -rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
}

__device__ void evalLocal(float *YnmShrd, float rho, float alpha, float *factShrd) {
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
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
    fact += 2;
  }
}

__device__ void P2M_core(float *target, float rho, float alpha, float beta, float source) {
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[NTERM];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int nn = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int mm = 0;
  for( int i=0; i<=nn; ++i ) mm += i;
  mm = threadIdx.x - mm;
  if( threadIdx.x >= NTERM ) nn = mm = 0;
  float x = cosf(alpha);
  float s = sqrtf(1 - x * x);
  fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<=mm; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    YnmShrd[i] = rhom * p * rsqrtf(factShrd[2*m]);
    float p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho;
    float rhon = rhom;
    for( int n=m+1; n<=nn; ++n ) {
      i = n * (n + 1) / 2 + m;
      YnmShrd[i] = rhon * p * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  int i = nn * (nn + 1) / 2 + mm;
  float ere = cosf(-mm * beta);
  float eim = sinf(-mm * beta);
  target[0] += source * YnmShrd[i] * ere;
  target[1] += source * YnmShrd[i] * eim;
}

__global__ void P2M_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float targetShrd[3];
  __shared__ float sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS;
  targetShrd[0] = targetGlob[2*itarget+0];
  targetShrd[1] = targetGlob[2*itarget+1];
  targetShrd[2] = targetGlob[2*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    int size  = rangeGlob[keys+3*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
      __syncthreads();
      for( int i=0; i<THREADS; ++i ) {
        float3 d;
        d.x = sourceShrd[4*i+0] - targetShrd[0];
        d.y = sourceShrd[4*i+1] - targetShrd[1];
        d.z = sourceShrd[4*i+2] - targetShrd[2];
        float rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        P2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
    sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
    sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
    sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      float3 d;
      d.x = sourceShrd[4*i+0] - targetShrd[0];
      d.y = sourceShrd[4*i+1] - targetShrd[1];
      d.z = sourceShrd[4*i+2] - targetShrd[2];
      float rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      P2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void M2M_core(float *target, float beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<=j; ++n ) {
    for( int m=-n; m<=min(k-1,n); ++m ) {
      if( j-n >= k-m ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
        float ere = cosf(-m * beta);
        float eim = sinf(-m * beta);
        float ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        float cnm = ODDEVEN((m-abs(m))/2+j);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] -= sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] += sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
    for( int m=k; m<=n; ++m ) {
      if( j-n >= m-k ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
        float ere = cosf(-m * beta);
        float eim = sinf(-m * beta);
        float ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        float cnm = ODDEVEN(k+j+m);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] += sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
  }
}

__global__ void M2M_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[P*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    M2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

void Kernel::M2M_CPU() {
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
            M += CJ->M[jnkms] * std::pow(I,double(m-abs(m))) * Ynm[nm]
               * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            M += std::conj(CJ->M[jnkms]) * Ynm[nm]
               * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
      }
      CI->M[jks] += M;
    }
  }
}

__device__ void M2L_core(float *target, float  beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
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
  float D0 = -constDevc[0];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[4*NTERM];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
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
            d.x += targetGlob[2*itarget+0] - sourceGlob[begin+0];
            d.y += targetGlob[2*itarget+1] - sourceGlob[begin+1];
            d.z += targetGlob[2*itarget+2] - sourceGlob[begin+2];
            float rho,alpha,beta;
            cart2sph(rho,alpha,beta,d.x,d.y,d.z);
            evalLocal(YnmShrd,rho,alpha,factShrd);
            M2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
          }
        }
      }
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void M2P_core(float *target, float r, float theta, float phi, float *factShrd, float *sourceShrd) {
  float x = cosf(theta);
  float y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  float s = sqrtf(1 - x * x);
  float spherical[3] = {0, 0, 0};
  float cartesian[3] = {0, 0, 0};
  float fact = 1;
  float pn = 1;
  float rhom = 1.0 / r;
  for( int m=0; m<P; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    float ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    float eim = sinf(m * phi);
    float anm = rhom * rsqrtf(factShrd[2*m]);
    float Ynm = anm * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    float YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    float realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    float imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom /= r;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
      target[0] += 2 * Ynm * realj;
      spherical[0] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      rhon /= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,spherical,cartesian);
  target[1] += cartesian[0];
  target[2] += cartesian[1];
  target[3] += cartesian[2];
}

__global__ void M2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float D0 = -constDevc[0];
  float targetPos[3];
  float target[4] = {0, 0, 0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetPos[0] = targetGlob[4*itarget+0];
  targetPos[1] = targetGlob[4*itarget+1];
  targetPos[2] = targetGlob[4*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
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
            d.x += targetPos[0] - sourceGlob[begin+0];
            d.y += targetPos[1] - sourceGlob[begin+1];
            d.z += targetPos[2] - sourceGlob[begin+2];
            float r,theta,phi;
            cart2sph(r,theta,phi,d.x,d.y,d.z);
            M2P_core(target,r,theta,phi,factShrd,sourceShrd);
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

__device__ inline void P2P_core(float *target, float *targetPos, float *sourceShrd, float3 d, int i) {
  d.x += targetPos[0];
  d.x -= sourceShrd[4*i+0];
  d.y += targetPos[1];
  d.y -= sourceShrd[4*i+1];
  d.z += targetPos[2];
  d.z -= sourceShrd[4*i+2];
  float invR = rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
  float invR3 = sourceShrd[4*i+3] * invR * invR * invR;
  target[0] += sourceShrd[4*i+3] * invR;
  target[1] -= d.x * invR3;
  target[2] -= d.y * invR3;
  target[3] -= d.z * invR3;
}

__global__ void P2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float D0 = -constDevc[0];
  float targetPos[3];
  float target[4] = {0, 0, 0, 0};
  __shared__ float sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetPos[0] = targetGlob[4*itarget+0];
  targetPos[1] = targetGlob[4*itarget+1];
  targetPos[2] = targetGlob[4*itarget+2];
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
                P2P_core(target,targetPos,sourceShrd,d,i);
              }
            }
          }
        }
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
    sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
    sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
    sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
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
              P2P_core(target,targetPos,sourceShrd,d,i);
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

__device__ void L2L_core(float *target, float beta, float *factShrd, float *YnmShrd, float *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  float ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=j+k-n; m<0; ++m ) {
      int nms = n * (n + 1) / 2 - m;
      int jnkm = (n - j) * (n - j) + n - j + m - k;
      float ere = cosf((m - k) * beta);
      float eim = sinf((m - k) * beta);
      float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      float cnm = ODDEVEN(k-n) * ajk / anm * YnmShrd[jnkm];
      float CnmReal = cnm * ere;
      float CnmImag = cnm * eim;
      target[0] += sourceShrd[2*nms+0] * CnmReal;
      target[0] += sourceShrd[2*nms+1] * CnmImag;
      target[1] += sourceShrd[2*nms+0] * CnmImag;
      target[1] -= sourceShrd[2*nms+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      if( n-j >= abs(m-k) ) {
        int nms = n * (n + 1) / 2 + m;
        int jnkm = (n - j) * (n - j) + n - j + m - k;
        float ere = cosf((m - k) * beta);
        float eim = sinf((m - k) * beta);
        float anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
        float cnm = ODDEVEN((m-k-abs(m-k)) / 2 - n);
        cnm *= ajk / anm * YnmShrd[jnkm];
        float CnmReal = cnm * ere;
        float CnmImag = cnm * eim;
        target[0] += sourceShrd[2*nms+0] * CnmReal;
        target[0] -= sourceShrd[2*nms+1] * CnmImag;
        target[1] += sourceShrd[2*nms+0] * CnmImag;
        target[1] += sourceShrd[2*nms+1] * CnmReal;
      }
    }
  }
}

__global__ void L2L_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float target[2] = {0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  __shared__ float YnmShrd[P*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    L2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void L2P_core(float *target, float r, float theta, float phi, float *factShrd, float *sourceShrd) {
  float x = cosf(theta);
  float y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  float s = sqrtf(1 - x * x);
  float spherical[3] = {0, 0, 0};
  float cartesian[3] = {0, 0, 0};
  float fact = 1;
  float pn = 1;
  float rhom = 1;
  for( int m=0; m<P; ++m ) {
    float p = pn;
    int i = m * (m + 1) / 2 + m;
    float ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    float eim = sinf(m * phi);
    float anm = rhom * rsqrtf(factShrd[2*m]);
    float Ynm = anm * p;
    float p1 = p;
    p = x * (2 * m + 1) * p;
    float YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    float realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    float imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] += 2 * m / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom *= r;
    float rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      float p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmTheta = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
      realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
      imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
      target[0] += 2 * Ynm * realj;
      spherical[0] += 2 * n / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      rhon *= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,spherical,cartesian);
  target[1] += cartesian[0];
  target[2] += cartesian[1];
  target[3] += cartesian[2];
}

__global__ void L2P_GPU(int *keysGlob, int *rangeGlob, float *targetGlob, float *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  float targetPos[3];
  float target[4] = {0, 0, 0, 0};
  __shared__ float sourceShrd[2*THREADS];
  __shared__ float factShrd[2*P];
  float fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetPos[0] = targetGlob[4*itarget+0];
  targetPos[1] = targetGlob[4*itarget+1];
  targetPos[2] = targetGlob[4*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    float3 d;
    d.x = targetPos[0] - sourceGlob[begin+0];
    d.y = targetPos[1] - sourceGlob[begin+1];
    d.z = targetPos[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    float r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    L2P_core(target,r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

void Kernel::initialize() {
  precalculate();
  startTimer("Init GPU     ");                                  // Start timer
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
  cudaPrintfInit();
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel::allocGPU() {
  cudaThreadSynchronize();
  startTimer("cudaMalloc   ");
  cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );
  cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );
  cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(float) );
  cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(float) );
  cudaThreadSynchronize();
  stopTimer("cudaMalloc   ");
  startTimer("cudaMemcpy   ");
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(float));
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
}

void Kernel::deallocGPU() {
  cudaThreadSynchronize();
  startTimer("cudaMemcpy   ");
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(float),cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  stopTimer("cudaMemcpy   ");
  startTimer("cudaFree     ");
  cudaFree(keysDevc);
  cudaFree(rangeDevc);
  cudaFree(targetDevc);
  cudaFree(sourceDevc);
  cudaThreadSynchronize();
  stopTimer("cudaFree     ");
}

void Kernel::P2M() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("P2M kernel   ");
  int numBlocks = keysHost.size();
  P2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("P2M kernel   ");
  deallocGPU();
}

void Kernel::M2M() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("M2M kernel   ");
  int numBlocks = keysHost.size();
  M2M_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("M2M kernel   ");
  deallocGPU();
}

void Kernel::M2L() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("M2L kernel   ");
  int numBlocks = keysHost.size();
  M2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("M2L kernel   ");
  deallocGPU();
}

void Kernel::M2P() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("M2P kernel   ");
  int numBlocks = keysHost.size();
  M2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("M2P kernel   ");
  deallocGPU();
}

void Kernel::P2P() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("P2P kernel   ");
  int numBlocks = keysHost.size();
  P2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("P2P kernel   ");
  deallocGPU();
}

void Kernel::L2L() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("L2L kernel   ");
  int numBlocks = keysHost.size();
  L2L_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("L2L kernel   ");
  deallocGPU();
}

void Kernel::L2P() {
  allocGPU();
  cudaThreadSynchronize();
  startTimer("L2P kernel   ");
  int numBlocks = keysHost.size();
  L2P_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);
  cudaThreadSynchronize();
  stopTimer("L2P kernel   ");
  deallocGPU();
}

void Kernel::finalize() {
  postcalculate();
  cudaPrintfDisplay(stdout, true);
  cudaPrintfEnd();
}
