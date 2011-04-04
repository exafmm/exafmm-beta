#include "kernel.h"
#include "spherical.h"
#include "laplace.h"
#include "pregpu.h"

void Kernel::LaplaceInit() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
#ifdef CUPRINTF
  cudaPrintfInit();                                             // Initialize cuPrintf
  cudaThreadSynchronize();                                      // Sync GPU threads
#endif
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel::LaplacePre() {
  prefactor = new double  [4*P2];
  Anm       = new double  [4*P2];
  Ynm       = new complex [4*P2];
  YnmTheta  = new complex [4*P2];
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

__device__ void LaplaceP2M_core(double *target, double rho, double alpha, double beta, double source) {
  __shared__ double factShrd[2*P];
  __shared__ double YnmShrd[NTERM];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int nn = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int mm = 0;
  for( int i=0; i<=nn; ++i ) mm += i;
  mm = threadIdx.x - mm;
  if( threadIdx.x >= NTERM ) nn = mm = 0;
  double x = cosf(alpha);
  double s = sqrtf(1 - x * x);
  fact = 1;
  double pn = 1;
  double rhom = 1;
  for( int m=0; m<=mm; ++m ) {
    double p = pn;
    int i = m * (m + 1) / 2 + m;
    YnmShrd[i] = rhom * p * rsqrtf(factShrd[2*m]);
    double p1 = p;
    p = x * (2 * m + 1) * p;
    rhom *= rho;
    double rhon = rhom;
    for( int n=m+1; n<=nn; ++n ) {
      i = n * (n + 1) / 2 + m;
      YnmShrd[i] = rhon * p * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      double p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  int i = nn * (nn + 1) / 2 + mm;
  double ere = cosf(-mm * beta);
  double eim = sinf(-mm * beta);
  target[0] += source * YnmShrd[i] * ere;
  target[1] += source * YnmShrd[i] * eim;
}

__global__ void LaplaceP2M_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[2] = {0, 0};
  __shared__ double targetShrd[3];
  __shared__ double sourceShrd[4*THREADS];
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
        double3 d;
        d.x = sourceShrd[4*i+0] - targetShrd[0];
        d.y = sourceShrd[4*i+1] - targetShrd[1];
        d.z = sourceShrd[4*i+2] - targetShrd[2];
        double rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        LaplaceP2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    if( threadIdx.x < size - iblok * THREADS ) {
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    }
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      double3 d;
      d.x = sourceShrd[4*i+0] - targetShrd[0];
      d.y = sourceShrd[4*i+1] - targetShrd[1];
      d.z = sourceShrd[4*i+2] - targetShrd[2];
      double rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      LaplaceP2M_core(target,rho,alpha,beta,sourceShrd[4*i+3]);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void LaplaceM2M_core(double *target, double beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  double ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<=j; ++n ) {
    for( int m=-n; m<=min(k-1,n); ++m ) {
      if( j-n >= k-m ) {
        int nm = n * n + n + m;
        int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
        double ere = cosf(-m * beta);
        double eim = sinf(-m * beta);
        double ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        double cnm = ODDEVEN((m-abs(m))/2+j);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        double CnmReal = cnm * ere;
        double CnmImag = cnm * eim;
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
        double ere = cosf(-m * beta);
        double eim = sinf(-m * beta);
        double ajnkm = rsqrtf(factShrd[j-n-k+m] * factShrd[j-n+k-m]);
        double cnm = ODDEVEN(k+j+m);
        cnm *= ajnkm / ajk * YnmShrd[nm];
        double CnmReal = cnm * ere;
        double CnmImag = cnm * eim;
        target[0] += sourceShrd[2*jnkms+0] * CnmReal;
        target[0] += sourceShrd[2*jnkms+1] * CnmImag;
        target[1] += sourceShrd[2*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[2*jnkms+1] * CnmReal;
      }
    }
  }
}

__global__ void LaplaceM2M_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[2] = {0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  __shared__ double YnmShrd[P*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    double3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    double rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    LaplaceM2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

void Kernel::LaplaceM2M_CPU() {
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

__device__ void LaplaceM2L_core(double *target, double  beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  double ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=-n; m<0; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 - m + k;
      double ere = cosf((m - k) * beta);
      double eim = sinf((m - k) * beta);
      double anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      double cnm = anm * ajk * YnmShrd[jnkm];
      double CnmReal = cnm * ere;
      double CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 - m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] += sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] -= sourceShrd[2*i+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      int jnkm = (j + n) * (j + n + 1) / 2 + abs(m - k);
      double ere = cosf((m - k) * beta);
      double eim = sinf((m - k) * beta);
      double anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      double cnm = ODDEVEN((abs(k - m) - k - m) / 2);
      cnm *= anm * ajk * YnmShrd[jnkm];
      double CnmReal = cnm * ere;
      double CnmImag = cnm * eim;
      int i = n * (n + 1) / 2 + m;
      target[0] += sourceShrd[2*i+0] * CnmReal;
      target[0] -= sourceShrd[2*i+1] * CnmImag;
      target[1] += sourceShrd[2*i+0] * CnmImag;
      target[1] += sourceShrd[2*i+1] * CnmReal;
    }
  }
}

__global__ void LaplaceM2L_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double target[2] = {0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  __shared__ double YnmShrd[4*NTERM];
  double fact = 1;
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
            double3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            d.x += targetGlob[2*itarget+0] - sourceGlob[begin+0];
            d.y += targetGlob[2*itarget+1] - sourceGlob[begin+1];
            d.z += targetGlob[2*itarget+2] - sourceGlob[begin+2];
            double rho,alpha,beta;
            cart2sph(rho,alpha,beta,d.x,d.y,d.z);
            evalLocal(YnmShrd,rho,alpha,factShrd);
            LaplaceM2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
          }
        }
      }
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void LaplaceM2P_core(double *target, double r, double theta, double phi, double *factShrd, double *sourceShrd) {
  double x = cosf(theta);
  double y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  double s = sqrtf(1 - x * x);
  double spherical[3] = {0, 0, 0};
  double cartesian[3] = {0, 0, 0};
  double fact = 1;
  double pn = 1;
  double rhom = 1.0 / r;
  for( int m=0; m<P; ++m ) {
    double p = pn;
    int i = m * (m + 1) / 2 + m;
    double ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    double eim = sinf(m * phi);
    double anm = rhom * rsqrtf(factShrd[2*m]);
    double Ynm = anm * p;
    double p1 = p;
    p = x * (2 * m + 1) * p;
    double YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    double realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    double imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom /= r;
    double rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      double p2 = p1;
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

__global__ void LaplaceM2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double targetX[3];
  double target[4] = {0, 0, 0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[4*itarget+0];
  targetX[1] = targetGlob[4*itarget+1];
  targetX[2] = targetGlob[4*itarget+2];
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
            double3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            d.x += targetX[0] - sourceGlob[begin+0];
            d.y += targetX[1] - sourceGlob[begin+1];
            d.z += targetX[2] - sourceGlob[begin+2];
            double r,theta,phi;
            cart2sph(r,theta,phi,d.x,d.y,d.z);
            LaplaceM2P_core(target,r,theta,phi,factShrd,sourceShrd);
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

__device__ inline void LaplaceP2P_core(double *target, double *targetX, double *sourceShrd, double3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[4*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[4*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[4*i+2];
  double invR = rsqrtf(d.x * d.x + d.y * d.y + d.z * d.z + EPS2);
  double invR3 = sourceShrd[4*i+3] * invR * invR * invR;
  target[0] += sourceShrd[4*i+3] * invR;
  target[1] -= d.x * invR3;
  target[2] -= d.y * invR3;
  target[3] -= d.z * invR3;
}

__global__ void LaplaceP2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double targetX[3];
  double target[4] = {0, 0, 0, 0};
  __shared__ double sourceShrd[4*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[4*itarget+0];
  targetX[1] = targetGlob[4*itarget+1];
  targetX[2] = targetGlob[4*itarget+2];
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
              double3 d;
              d.x = ix * D0;
              d.y = iy * D0;
              d.z = iz * D0;
#pragma unroll 64
              for( int i=0; i<THREADS; ++i ) {
                LaplaceP2P_core(target,targetX,sourceShrd,d,i);
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
      sourceShrd[4*threadIdx.x+0] = sourceGlob[4*isource+0];
      sourceShrd[4*threadIdx.x+1] = sourceGlob[4*isource+1];
      sourceShrd[4*threadIdx.x+2] = sourceGlob[4*isource+2];
      sourceShrd[4*threadIdx.x+3] = sourceGlob[4*isource+3];
    }
    __syncthreads();
    int I = 0;
    int icounter=0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz, ++I ) {
          if( Iperiodic & (1 << I) ) {
            icounter++;
            double3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            for( int i=0; i<size-iblok*THREADS; ++i ) {
              LaplaceP2P_core(target,targetX,sourceShrd,d,i);
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

__device__ void LaplaceL2L_core(double *target, double beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
  int j = floor(sqrtf(2*threadIdx.x+0.25)-0.5);
  int k = 0;
  for( int i=0; i<=j; ++i ) k += i;
  k = threadIdx.x - k;
  if( threadIdx.x >= NTERM ) j = k = 0;
  double ajk = ODDEVEN(j) * rsqrtf(factShrd[j-k] * factShrd[j+k]);
  for( int n=0; n<P; ++n ) {
    for( int m=j+k-n; m<0; ++m ) {
      int nms = n * (n + 1) / 2 - m;
      int jnkm = (n - j) * (n - j) + n - j + m - k;
      double ere = cosf((m - k) * beta);
      double eim = sinf((m - k) * beta);
      double anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
      double cnm = ODDEVEN(k-n) * ajk / anm * YnmShrd[jnkm];
      double CnmReal = cnm * ere;
      double CnmImag = cnm * eim;
      target[0] += sourceShrd[2*nms+0] * CnmReal;
      target[0] += sourceShrd[2*nms+1] * CnmImag;
      target[1] += sourceShrd[2*nms+0] * CnmImag;
      target[1] -= sourceShrd[2*nms+1] * CnmReal;
    }
    for( int m=0; m<=n; ++m ) {
      if( n-j >= abs(m-k) ) {
        int nms = n * (n + 1) / 2 + m;
        int jnkm = (n - j) * (n - j) + n - j + m - k;
        double ere = cosf((m - k) * beta);
        double eim = sinf((m - k) * beta);
        double anm = rsqrtf(factShrd[n-m] * factShrd[n+m]);
        double cnm = ODDEVEN((m-k-abs(m-k)) / 2 - n);
        cnm *= ajk / anm * YnmShrd[jnkm];
        double CnmReal = cnm * ere;
        double CnmImag = cnm * eim;
        target[0] += sourceShrd[2*nms+0] * CnmReal;
        target[0] -= sourceShrd[2*nms+1] * CnmImag;
        target[1] += sourceShrd[2*nms+0] * CnmImag;
        target[1] += sourceShrd[2*nms+1] * CnmReal;
      }
    }
  }
}

__global__ void LaplaceL2L_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[2] = {0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  __shared__ double YnmShrd[P*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS;
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    double3 d;
    d.x = targetGlob[2*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[2*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[2*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    double rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    LaplaceL2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[2*itarget+0] = target[0];
  targetGlob[2*itarget+1] = target[1];
}

__device__ void LaplaceL2P_core(double *target, double r, double theta, double phi, double *factShrd, double *sourceShrd) {
  double x = cosf(theta);
  double y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  double s = sqrtf(1 - x * x);
  double spherical[3] = {0, 0, 0};
  double cartesian[3] = {0, 0, 0};
  double fact = 1;
  double pn = 1;
  double rhom = 1;
  for( int m=0; m<P; ++m ) {
    double p = pn;
    int i = m * (m + 1) / 2 + m;
    double ere = cosf(m * phi);
    if( m == 0 ) ere = 0.5;
    double eim = sinf(m * phi);
    double anm = rhom * rsqrtf(factShrd[2*m]);
    double Ynm = anm * p;
    double p1 = p;
    p = x * (2 * m + 1) * p;
    double YnmTheta = anm * (p - (m + 1) * x * p1) / y;
    double realj = ere * sourceShrd[2*i+0] - eim * sourceShrd[2*i+1];
    double imagj = eim * sourceShrd[2*i+0] + ere * sourceShrd[2*i+1];
    target[0] += 2 * Ynm * realj;
    spherical[0] += 2 * m / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    rhom *= r;
    double rhon = rhom;
    for( int n=m+1; n<P; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      Ynm = anm * p;
      double p2 = p1;
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

__global__ void LaplaceL2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double targetX[3];
  double target[4] = {0, 0, 0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[4*itarget+0];
  targetX[1] = targetGlob[4*itarget+1];
  targetX[2] = targetGlob[4*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    double3 d;
    d.x = targetX[0] - sourceGlob[begin+0];
    d.y = targetX[1] - sourceGlob[begin+1];
    d.z = targetX[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[2*threadIdx.x+0] = sourceGlob[begin+2*threadIdx.x+3];
      sourceShrd[2*threadIdx.x+1] = sourceGlob[begin+2*threadIdx.x+4];
    }
    __syncthreads();
    double r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    LaplaceL2P_core(target,r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[4*itarget+0] = target[0];
  targetGlob[4*itarget+1] = target[1];
  targetGlob[4*itarget+2] = target[2];
  targetGlob[4*itarget+3] = target[3];
}

void Kernel::LaplacePost() {
  delete[] prefactor;
  delete[] Anm;
  delete[] Ynm;
  delete[] YnmTheta;
  delete[] Cnm;
}

void Kernel::LaplaceFinal() {
#ifdef CUPRINTF
  cudaPrintfDisplay(stdout, true);                              // Print cuPrintf buffer to display
  cudaPrintfEnd();                                              // Finalize cuPrintf
#endif
}

#include "gpu.h"
