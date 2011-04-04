#include "kernel.h"
#include "spherical.h"
#include "stretching.h"
#include "pregpu.h"

void Kernel::StretchingInit() {
  startTimer("Init GPU     ");                                  // Start timer
  cudaThreadExit();                                             // Exit GPU thread
  cudaSetDevice(MPIRANK % GPUS);                                // Set GPU device
#ifdef CUPRINTF
  cudaPrintfInit();                                             // Initialize cuPrintf
#endif
  cudaThreadSynchronize();                                      // Sync GPU threads
  stopTimer("Init GPU     ",MPIRANK==0);                        // Stop timer & print
  eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel::StretchingPre() {
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

__device__ void StretchingP2M_core(double *target, double rho, double alpha, double beta, double *sourceShrd, int ithread) {
  __shared__ double factShrd[2*P];
  __shared__ double YnmShrd[NTERM];
  __shared__ double YnmAlphaShrd[NTERM];
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
  double y = sinf(alpha);
  if( fabs(y) < EPS ) y = 1 / EPS;
  double s = sqrtf(1 - x * x);
  fact = 1;
  double pn = 1;
  double rhom = 1;
  for( int m=0; m<=mm; ++m ) {
    double p = pn;
    int i = m * (m + 1) / 2 + m;
    double anm = rhom * rsqrtf(factShrd[2*m]);
    YnmShrd[i] = anm * p;
    double p1 = p;
    p = x * (2 * m + 1) * p;
    YnmAlphaShrd[i] = anm * (p - (m + 1) * x * p1) / y;
    rhom *= rho;
    double rhon = rhom;
    for( int n=m+1; n<=nn; ++n ) {
      i = n * (n + 1) / 2 + m;
      anm = rhon * rsqrtf(factShrd[n+m] / factShrd[n-m]);
      YnmShrd[i] = anm * p;
      double p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmAlphaShrd[i] = anm * ((n - m + 1) * p - (n + 1) * x * p1) / y;
      rhon *= rho;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  int i = nn * (nn + 1) / 2 + mm;
  double ere = cosf(-mm * beta);
  double eim = sinf(-mm * beta);
  double spherical[6];
  double cartesian[6];
  spherical[0] = YnmShrd[i] * nn / rho * ere;
  spherical[1] = YnmAlphaShrd[i] * ere;
  spherical[2] = YnmShrd[i] * mm * eim;
  spherical[3] = YnmShrd[i] * nn / rho * eim;
  spherical[4] = YnmAlphaShrd[i] * eim;
  spherical[5] = -YnmShrd[i] * mm * ere;
  sph2cart(rho,alpha,beta,&spherical[0],&cartesian[0]);
  sph2cart(rho,alpha,beta,&spherical[3],&cartesian[3]);
  target[0] += sourceShrd[6*ithread+4] * cartesian[2] - sourceShrd[6*ithread+5] * cartesian[1];
  target[1] += sourceShrd[6*ithread+4] * cartesian[5] - sourceShrd[6*ithread+5] * cartesian[4];
  target[2] += sourceShrd[6*ithread+5] * cartesian[0] - sourceShrd[6*ithread+3] * cartesian[2];
  target[3] += sourceShrd[6*ithread+5] * cartesian[3] - sourceShrd[6*ithread+3] * cartesian[5];
  target[4] += sourceShrd[6*ithread+3] * cartesian[1] - sourceShrd[6*ithread+4] * cartesian[0];
  target[5] += sourceShrd[6*ithread+3] * cartesian[4] - sourceShrd[6*ithread+4] * cartesian[3];
}

__global__ void StretchingP2M_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ double targetShrd[3];
  __shared__ double sourceShrd[6*THREADS];
  int itarget = blockIdx.x * THREADS;
  targetShrd[0] = targetGlob[6*itarget+0];
  targetShrd[1] = targetGlob[6*itarget+1];
  targetShrd[2] = targetGlob[6*itarget+2];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    int size  = rangeGlob[keys+3*ilist+2];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[6*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[7*isource+5];
      __syncthreads();
      for( int i=0; i<THREADS; ++i ) {
        double3 d;
        d.x = sourceShrd[6*i+0] - targetShrd[0];
        d.y = sourceShrd[6*i+1] - targetShrd[1];
        d.z = sourceShrd[6*i+2] - targetShrd[2];
        double rho,alpha,beta;
        cart2sph(rho,alpha,beta,d.x,d.y,d.z);
        StretchingP2M_core(target,rho,alpha,beta,sourceShrd,i);
      }
    }
    int iblok = (size-1)/THREADS;
    int isource = begin + iblok * THREADS + threadIdx.x;
    __syncthreads();
    if( threadIdx.x < size - iblok * THREADS ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[7*isource+5];
    }
    __syncthreads();
    for( int i=0; i<size-iblok*THREADS; ++i ) {
      double3 d;
      d.x = sourceShrd[6*i+0] - targetShrd[0];
      d.y = sourceShrd[6*i+1] - targetShrd[1];
      d.z = sourceShrd[6*i+2] - targetShrd[2];
      double rho,alpha,beta;
      cart2sph(rho,alpha,beta,d.x,d.y,d.z);
      StretchingP2M_core(target,rho,alpha,beta,sourceShrd,i);
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

__device__ void StretchingM2M_core(double *target, double beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
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
        target[0] += sourceShrd[6*jnkms+0] * CnmReal;
        target[0] -= sourceShrd[6*jnkms+1] * CnmImag;
        target[1] += sourceShrd[6*jnkms+0] * CnmImag;
        target[1] += sourceShrd[6*jnkms+1] * CnmReal;
        target[2] += sourceShrd[6*jnkms+2] * CnmReal;
        target[2] -= sourceShrd[6*jnkms+3] * CnmImag;
        target[3] += sourceShrd[6*jnkms+2] * CnmImag;
        target[3] += sourceShrd[6*jnkms+3] * CnmReal;
        target[4] += sourceShrd[6*jnkms+4] * CnmReal;
        target[4] -= sourceShrd[6*jnkms+5] * CnmImag;
        target[5] += sourceShrd[6*jnkms+4] * CnmImag;
        target[5] += sourceShrd[6*jnkms+5] * CnmReal;
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
        target[0] += sourceShrd[6*jnkms+0] * CnmReal;
        target[0] += sourceShrd[6*jnkms+1] * CnmImag;
        target[1] += sourceShrd[6*jnkms+0] * CnmImag;
        target[1] -= sourceShrd[6*jnkms+1] * CnmReal;
        target[2] += sourceShrd[6*jnkms+2] * CnmReal;
        target[2] += sourceShrd[6*jnkms+3] * CnmImag;
        target[3] += sourceShrd[6*jnkms+2] * CnmImag;
        target[3] -= sourceShrd[6*jnkms+3] * CnmReal;
        target[4] += sourceShrd[6*jnkms+4] * CnmReal;
        target[4] += sourceShrd[6*jnkms+5] * CnmImag;
        target[5] += sourceShrd[6*jnkms+4] * CnmImag;
        target[5] -= sourceShrd[6*jnkms+5] * CnmReal;
      }
    }
  }
}

__global__ void StretchingM2M_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ double sourceShrd[6*THREADS];
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
    d.x = targetGlob[6*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[6*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[6*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    double rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    StretchingM2M_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

void Kernel::StretchingM2M_CPU() {
  vect dist = CI->X - CJ->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,-beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex M[3] = {0, 0, 0};
      for( int n=0; n<=j; ++n ) {
        for( int m=-n; m<=std::min(k-1,n); ++m ) {
          if( j-n >= k-m ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            const int nm    = n * n + n + m;
            for( int d=0; d!=3; ++d ) {
              M[d] += CJ->M[3*jnkms+d] * std::pow(I,double(m-abs(m))) * Ynm[nm]
                    * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            for( int d=0; d!=3; ++d ) {
              M[d] += std::conj(CJ->M[3*jnkms+d]) * Ynm[nm]
                    * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
      }
      for( int d=0; d!=3; ++d ) {
        CI->M[3*jks+d] += M[d];
      }
    }
  }
}

__device__ void StretchingM2L_core(double *target, double  beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
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
      target[0] += sourceShrd[6*i+0] * CnmReal;
      target[0] += sourceShrd[6*i+1] * CnmImag;
      target[1] += sourceShrd[6*i+0] * CnmImag;
      target[1] -= sourceShrd[6*i+1] * CnmReal;
      target[2] += sourceShrd[6*i+2] * CnmReal;
      target[2] += sourceShrd[6*i+3] * CnmImag;
      target[3] += sourceShrd[6*i+2] * CnmImag;
      target[3] -= sourceShrd[6*i+3] * CnmReal;
      target[4] += sourceShrd[6*i+4] * CnmReal;
      target[4] += sourceShrd[6*i+5] * CnmImag;
      target[5] += sourceShrd[6*i+4] * CnmImag;
      target[5] -= sourceShrd[6*i+5] * CnmReal;
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
      target[0] += sourceShrd[6*i+0] * CnmReal;
      target[0] -= sourceShrd[6*i+1] * CnmImag;
      target[1] += sourceShrd[6*i+0] * CnmImag;
      target[1] += sourceShrd[6*i+1] * CnmReal;
      target[2] += sourceShrd[6*i+2] * CnmReal;
      target[2] -= sourceShrd[6*i+3] * CnmImag;
      target[3] += sourceShrd[6*i+2] * CnmImag;
      target[3] += sourceShrd[6*i+3] * CnmReal;
      target[4] += sourceShrd[6*i+4] * CnmReal;
      target[4] -= sourceShrd[6*i+5] * CnmImag;
      target[5] += sourceShrd[6*i+4] * CnmImag;
      target[5] += sourceShrd[6*i+5] * CnmReal;
    }
  }
}

__global__ void StretchingM2L_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ double sourceShrd[6*THREADS];
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
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
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
            d.x += targetGlob[6*itarget+0] - sourceGlob[begin+0];
            d.y += targetGlob[6*itarget+1] - sourceGlob[begin+1];
            d.z += targetGlob[6*itarget+2] - sourceGlob[begin+2];
            double rho,alpha,beta;
            cart2sph(rho,alpha,beta,d.x,d.y,d.z);
            evalLocal(YnmShrd,rho,alpha,factShrd);
            StretchingM2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
          }
        }
      }
    }
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

__device__ void StretchingM2P_core(double *target, double *targetQ, double r, double theta, double phi,
                                   double *factShrd, double *sourceShrd) {
  double x = cosf(theta);
  double y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  double s = sqrtf(1 - x * x);
  double spherical[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double cartesian[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    double realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
    double imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
    spherical[0] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
    imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
    spherical[3] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[4] += 2 * YnmTheta * realj;
    spherical[5] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
    imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
    spherical[6] -= 2 * (m + 1) / r * Ynm * realj;
    spherical[7] += 2 * YnmTheta * realj;
    spherical[8] -= 2 * m * Ynm * imagj;
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
      realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
      imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
      spherical[0] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
      imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
      spherical[3] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[4] += 2 * YnmTheta * realj;
      spherical[5] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
      imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
      spherical[6] -= 2 * (n + 1) / r * Ynm * realj;
      spherical[7] += 2 * YnmTheta * realj;
      spherical[8] -= 2 * m * Ynm * imagj;
      rhon /= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,&spherical[0],&cartesian[0]);
  sph2cart(r,theta,phi,&spherical[3],&cartesian[3]);
  sph2cart(r,theta,phi,&spherical[6],&cartesian[6]);
  target[0] -= 0.25 / M_PI * (targetQ[0] * cartesian[0] + targetQ[1] * cartesian[1] + targetQ[2] * cartesian[2]);
  target[1] -= 0.25 / M_PI * (targetQ[0] * cartesian[3] + targetQ[1] * cartesian[4] + targetQ[2] * cartesian[5]);
  target[2] -= 0.25 / M_PI * (targetQ[0] * cartesian[6] + targetQ[1] * cartesian[7] + targetQ[2] * cartesian[8]);
}

__global__ void StretchingM2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double targetX[3], targetQ[3];
  double target[3] = {0, 0, 0};
  __shared__ double sourceShrd[2*THREADS];
  __shared__ double factShrd[2*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  targetQ[0] = targetGlob[6*itarget+3];
  targetQ[1] = targetGlob[6*itarget+4];
  targetQ[2] = targetGlob[6*itarget+5];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
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
            StretchingM2P_core(target,targetQ,r,theta,phi,factShrd,sourceShrd);
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

__device__ inline void StretchingP2P_core(double *target, double *targetX, double *targetQ, double *sourceShrd, double3 d, int i) {
  d.x += targetX[0];
  d.x -= sourceShrd[7*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[7*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[7*i+2];
  double R2 = d.x * d.x + d.y * d.y + d.z * d.z + EPS2;
#if 0
  double S2 = 2 * sourceShrd[7*i+6] * sourceShrd[7*i+6];
  double RS = R2 / S2;
  double cutoff = 0.25 / M_PI / R2 / sqrtf(R2) * (erff( sqrtf(RS) )
               - sqrtf(4 / M_PI * RS) * expf(-RS));
  target[0] += (targetQ[1] * sourceShrd[7*i+5] - targetQ[2] * sourceShrd[7*i+4]) * cutoff;
  target[1] += (targetQ[2] * sourceShrd[7*i+3] - targetQ[0] * sourceShrd[7*i+5]) * cutoff;
  target[2] += (targetQ[0] * sourceShrd[7*i+4] - targetQ[1] * sourceShrd[7*i+3]) * cutoff;
  cutoff = 0.25 / M_PI / R2 / R2 / sqrtf(R2) * (3 * erff( sqrtf(RS) )
         - (2 * RS + 3) * sqrtf(4 / M_PI * RS) * expf(-RS))
         * (targetQ[0] * d.x + targetQ[1] * d.y + targetQ[2] * d.z);
  target[0] += (sourceShrd[7*i+4] * d.z - sourceShrd[7*i+5] * d.y) * cutoff;
  target[1] += (sourceShrd[7*i+5] * d.x - sourceShrd[7*i+3] * d.z) * cutoff;
  target[2] += (sourceShrd[7*i+3] * d.y - sourceShrd[7*i+4] * d.x) * cutoff;
#else
  const double SQRT4PI = M_2_SQRTPI;
  const double FOURPI = 0.25 * M_1_PI;
  double SQRT_R2_1 = rsqrtf(R2);
  double RS = R2 * sourceShrd[7*i+6];
  double SQRT_RS = sqrtf(RS);
  double z = SQRT_RS,t,ERF_SQRT_RS;
  (t)=1.0f/(1.0f+0.5f*(z));
  ERF_SQRT_RS=1.0f - (t)*expf(-(z)*(z)-1.26551223f+(t)*(1.00002368f+(t)*(0.37409196f+(t)*(0.09678418f+
      (t)*(-0.18628806f+(t)*(0.27886807f+(t)*(-1.13520398f+(t)*(1.48851587f+
      (t)*(-0.82215223f+(t)*0.17087277f)))))))));
  double EXP_RS = expf(-RS);
  double cutoff = FOURPI * SQRT_R2_1 * SQRT_R2_1 * SQRT_R2_1 * (ERF_SQRT_RS
               - SQRT4PI * SQRT_RS * EXP_RS);
  target[0] += (targetQ[1] * sourceShrd[7*i+5] - targetQ[2] * sourceShrd[7*i+4]) * cutoff;
  target[1] += (targetQ[2] * sourceShrd[7*i+3] - targetQ[0] * sourceShrd[7*i+5]) * cutoff;
  target[2] += (targetQ[0] * sourceShrd[7*i+4] - targetQ[1] * sourceShrd[7*i+3]) * cutoff;
  double cutoff2 = FOURPI * SQRT_R2_1 * SQRT_R2_1 * SQRT_R2_1 * SQRT_R2_1 * SQRT_R2_1 * (3.0f * ERF_SQRT_RS
         - (2.0f * RS + 3.0f) * SQRT4PI * SQRT_RS * EXP_RS)
         * (targetQ[0] * d.x + targetQ[1] * d.y + targetQ[2] * d.z);
  target[0] += (sourceShrd[7*i+4] * d.z - sourceShrd[7*i+5] * d.y) * cutoff2;
  target[1] += (sourceShrd[7*i+5] * d.x - sourceShrd[7*i+3] * d.z) * cutoff2;
  target[2] += (sourceShrd[7*i+3] * d.y - sourceShrd[7*i+4] * d.x) * cutoff2;
#endif
}

__global__ void StretchingP2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double D0 = -constDevc[0];
  double targetX[3], targetQ[3];
  double target[3] = {0, 0, 0};
  __shared__ double sourceShrd[7*THREADS];
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  targetQ[0] = targetGlob[6*itarget+3];
  targetQ[1] = targetGlob[6*itarget+4];
  targetQ[2] = targetGlob[6*itarget+5];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin     = rangeGlob[keys+3*ilist+1];
    int size      = rangeGlob[keys+3*ilist+2];
    int Iperiodic = rangeGlob[keys+3*ilist+3];
    for( int iblok=0; iblok<(size-1)/THREADS; ++iblok ) {
      int isource = begin + iblok * THREADS + threadIdx.x;
      __syncthreads();
      sourceShrd[7*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[7*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[7*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[7*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[7*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[7*threadIdx.x+5] = sourceGlob[7*isource+5];
//      sourceShrd[7*threadIdx.x+6] = sourceGlob[7*isource+6];
      sourceShrd[7*threadIdx.x+6] = 0.5f / (sourceGlob[7*isource+6] * sourceGlob[7*isource+6]);
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
                StretchingP2P_core(target,targetX,targetQ,sourceShrd,d,i);
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
      sourceShrd[7*threadIdx.x+0] = sourceGlob[7*isource+0];
      sourceShrd[7*threadIdx.x+1] = sourceGlob[7*isource+1];
      sourceShrd[7*threadIdx.x+2] = sourceGlob[7*isource+2];
      sourceShrd[7*threadIdx.x+3] = sourceGlob[7*isource+3];
      sourceShrd[7*threadIdx.x+4] = sourceGlob[7*isource+4];
      sourceShrd[7*threadIdx.x+5] = sourceGlob[7*isource+5];
//      sourceShrd[7*threadIdx.x+6] = sourceGlob[7*isource+6];
      sourceShrd[7*threadIdx.x+6] = 0.5f / (sourceGlob[7*isource+6] * sourceGlob[7*isource+6]);
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
              StretchingP2P_core(target,targetX,targetQ,sourceShrd,d,i);
            }
          }
        }
      }
    }
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

__device__ void StretchingL2L_core(double *target, double beta, double *factShrd, double *YnmShrd, double *sourceShrd) {
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
      target[0] += sourceShrd[6*nms+0] * CnmReal;
      target[0] += sourceShrd[6*nms+1] * CnmImag;
      target[1] += sourceShrd[6*nms+0] * CnmImag;
      target[1] -= sourceShrd[6*nms+1] * CnmReal;
      target[2] += sourceShrd[6*nms+2] * CnmReal;
      target[2] += sourceShrd[6*nms+3] * CnmImag;
      target[3] += sourceShrd[6*nms+2] * CnmImag;
      target[3] -= sourceShrd[6*nms+3] * CnmReal;
      target[4] += sourceShrd[6*nms+4] * CnmReal;
      target[4] += sourceShrd[6*nms+5] * CnmImag;
      target[5] += sourceShrd[6*nms+4] * CnmImag;
      target[5] -= sourceShrd[6*nms+5] * CnmReal;
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
        target[0] += sourceShrd[6*nms+0] * CnmReal;
        target[0] -= sourceShrd[6*nms+1] * CnmImag;
        target[1] += sourceShrd[6*nms+0] * CnmImag;
        target[1] += sourceShrd[6*nms+1] * CnmReal;
        target[2] += sourceShrd[6*nms+2] * CnmReal;
        target[2] -= sourceShrd[6*nms+3] * CnmImag;
        target[3] += sourceShrd[6*nms+2] * CnmImag;
        target[3] += sourceShrd[6*nms+3] * CnmReal;
        target[4] += sourceShrd[6*nms+4] * CnmReal;
        target[4] -= sourceShrd[6*nms+5] * CnmImag;
        target[5] += sourceShrd[6*nms+4] * CnmImag;
        target[5] += sourceShrd[6*nms+5] * CnmReal;
      }
    }
  }
}

__global__ void StretchingL2L_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double target[6] = {0, 0, 0, 0, 0, 0};
  __shared__ double sourceShrd[6*THREADS];
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
    d.x = targetGlob[6*itarget+0] - sourceGlob[begin+0];
    d.y = targetGlob[6*itarget+1] - sourceGlob[begin+1];
    d.z = targetGlob[6*itarget+2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    double rho,alpha,beta;
    cart2sph(rho,alpha,beta,d.x,d.y,d.z);
    evalMultipole(YnmShrd,rho,alpha,factShrd);
    StretchingL2L_core(target,beta,factShrd,YnmShrd,sourceShrd);
  }
  itarget = blockIdx.x * THREADS + threadIdx.x;
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
  targetGlob[6*itarget+3] = target[3];
  targetGlob[6*itarget+4] = target[4];
  targetGlob[6*itarget+5] = target[5];
}

__device__ void StretchingL2P_core(double *target, double *targetQ, double r, double theta, double phi,
                                   double *factShrd, double *sourceShrd) {
  double x = cosf(theta);
  double y = sinf(theta);
  if( fabs(y) < EPS ) y = 1 / EPS;
  double s = sqrtf(1 - x * x);
  double spherical[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double cartesian[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    double realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
    double imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
    spherical[0] += 2 * m / r * Ynm * realj;
    spherical[1] += 2 * YnmTheta * realj;
    spherical[2] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
    imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
    spherical[3] += 2 * m / r * Ynm * realj;
    spherical[4] += 2 * YnmTheta * realj;
    spherical[5] -= 2 * m * Ynm * imagj;
    realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
    imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
    spherical[6] += 2 * m / r * Ynm * realj;
    spherical[7] += 2 * YnmTheta * realj;
    spherical[8] -= 2 * m * Ynm * imagj;
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
      realj = ere * sourceShrd[6*i+0] - eim * sourceShrd[6*i+1];
      imagj = eim * sourceShrd[6*i+0] + ere * sourceShrd[6*i+1];
      spherical[0] += 2 * n / r * Ynm * realj;
      spherical[1] += 2 * YnmTheta * realj;
      spherical[2] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+2] - eim * sourceShrd[6*i+3];
      imagj = eim * sourceShrd[6*i+2] + ere * sourceShrd[6*i+3];
      spherical[3] += 2 * n / r * Ynm * realj;
      spherical[4] += 2 * YnmTheta * realj;
      spherical[5] -= 2 * m * Ynm * imagj;
      realj = ere * sourceShrd[6*i+4] - eim * sourceShrd[6*i+5];
      imagj = eim * sourceShrd[6*i+4] + ere * sourceShrd[6*i+5];
      spherical[6] += 2 * n / r * Ynm * realj;
      spherical[7] += 2 * YnmTheta * realj;
      spherical[8] -= 2 * m * Ynm * imagj;
      rhon *= r;
    }
    pn = -pn * fact * s;
    fact += 2;
  }
  sph2cart(r,theta,phi,&spherical[0],&cartesian[0]);
  sph2cart(r,theta,phi,&spherical[3],&cartesian[3]);
  sph2cart(r,theta,phi,&spherical[6],&cartesian[6]);
  target[0] -= 0.25 / M_PI * (targetQ[0] * cartesian[0] + targetQ[1] * cartesian[1] + targetQ[2] * cartesian[2]);
  target[1] -= 0.25 / M_PI * (targetQ[0] * cartesian[3] + targetQ[1] * cartesian[4] + targetQ[2] * cartesian[5]);
  target[2] -= 0.25 / M_PI * (targetQ[0] * cartesian[6] + targetQ[1] * cartesian[7] + targetQ[2] * cartesian[8]);
}

__global__ void StretchingL2P_GPU(int *keysGlob, int *rangeGlob, double *targetGlob, double *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  double targetX[3], targetQ[3];
  double target[3] = {0, 0, 0};
  __shared__ double sourceShrd[6*THREADS];
  __shared__ double factShrd[2*P];
  double fact = 1;
  for( int i=0; i<2*P; ++i ) {
    factShrd[i] = fact;
    fact *= i + 1;
  }
  int itarget = blockIdx.x * THREADS + threadIdx.x;
  targetX[0] = targetGlob[6*itarget+0];
  targetX[1] = targetGlob[6*itarget+1];
  targetX[2] = targetGlob[6*itarget+2];
  targetQ[0] = targetGlob[6*itarget+3];
  targetQ[1] = targetGlob[6*itarget+4];
  targetQ[2] = targetGlob[6*itarget+5];
  for( int ilist=0; ilist<numList; ++ilist ) {
    int begin = rangeGlob[keys+3*ilist+1];
    double3 d;
    d.x = targetX[0] - sourceGlob[begin+0];
    d.y = targetX[1] - sourceGlob[begin+1];
    d.z = targetX[2] - sourceGlob[begin+2];
    __syncthreads();
    if( threadIdx.x < NTERM ) {
      sourceShrd[6*threadIdx.x+0] = sourceGlob[begin+6*threadIdx.x+3];
      sourceShrd[6*threadIdx.x+1] = sourceGlob[begin+6*threadIdx.x+4];
      sourceShrd[6*threadIdx.x+2] = sourceGlob[begin+6*threadIdx.x+5];
      sourceShrd[6*threadIdx.x+3] = sourceGlob[begin+6*threadIdx.x+6];
      sourceShrd[6*threadIdx.x+4] = sourceGlob[begin+6*threadIdx.x+7];
      sourceShrd[6*threadIdx.x+5] = sourceGlob[begin+6*threadIdx.x+8];
    }
    __syncthreads();
    double r,theta,phi;
    cart2sph(r,theta,phi,d.x,d.y,d.z);
    StretchingL2P_core(target,targetQ,r,theta,phi,factShrd,sourceShrd);
  }
  targetGlob[6*itarget+0] = target[0];
  targetGlob[6*itarget+1] = target[1];
  targetGlob[6*itarget+2] = target[2];
}

void Kernel::StretchingPost() {
  delete[] prefactor;
  delete[] Anm;
  delete[] Ynm;
  delete[] YnmTheta;
  delete[] Cnm;
}

void Kernel::StretchingFinal() {
#ifdef CUPRINTF
  cudaPrintfDisplay(stdout, true);                              // Print cuPrintf buffer to display
  cudaPrintfEnd();                                              // Finalize cuPrintf
#endif
}

#include "gpu.h"
