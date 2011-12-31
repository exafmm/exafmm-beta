/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

__device__ inline void LaplaceEwaldReal_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd,
                                             float3 d, int i, gpureal alpha) {
  d.x += targetX[0];
  d.x -= sourceShrd[4*i+0];
  d.y += targetX[1];
  d.y -= sourceShrd[4*i+1];
  d.z += targetX[2];
  d.z -= sourceShrd[4*i+2];
  gpureal R2 = d.x * d.x + d.y * d.y + d.z * d.z;
  if( R2 == 0 ) return;
  gpureal R2s = R2 * alpha * alpha;
  gpureal Rs = sqrtf(R2s);
  gpureal invRs = 1 / Rs;
  gpureal invR2s = invRs * invRs;
  gpureal invR3s = invR2s * invRs;
  gpureal ptmp = sourceShrd[4*i+3] * erfcf(Rs) * invRs * alpha;
  gpureal ftmp = sourceShrd[4*i+3] * (M_2_SQRTPI * expf(-R2s) * invR2s + erfcf(Rs) * invR3s);
  ftmp *= alpha * alpha * alpha;
  target[0] += ptmp;
  target[1] -= d.x * ftmp;
  target[2] -= d.y * ftmp;
  target[3] -= d.z * ftmp;
}

__global__ void LaplaceEwaldReal_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob) {
  int keys = keysGlob[blockIdx.x];
  int numList = rangeGlob[keys];
  gpureal D0 = -constDevc[0];
  gpureal alpha = constDevc[1];
  gpureal targetX[3];
  gpureal target[4] = {0, 0, 0, 0};
  __shared__ gpureal sourceShrd[4*THREADS];
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
              float3 d;
              d.x = ix * D0;
              d.y = iy * D0;
              d.z = iz * D0;
#pragma unroll 64
              for( int i=0; i<THREADS; ++i ) {
                LaplaceEwaldReal_core(target,targetX,sourceShrd,d,i,alpha);
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
            float3 d;
            d.x = ix * D0;
            d.y = iy * D0;
            d.z = iz * D0;
            for( int i=0; i<size-iblok*THREADS; ++i ) {
              LaplaceEwaldReal_core(target,targetX,sourceShrd,d,i,alpha);
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

template<>
void Kernel<Laplace>::EwaldReal() {
  cudaThreadSynchronize();
  startTimer("EwaldReal GPU");
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    LaplaceEwaldReal_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer("EwaldReal GPU");\
}

namespace {
void dft(Ewalds &ewalds, Bodies &bodies, real R0) {
  real scale = M_PI / R0;
  for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
    E->REAL = E->IMAG = 0;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      real th = 0;
      for( int d=0; d<3; d++ ) th += E->K[d] * B->X[d] * scale;
      E->REAL += B->SRC * cos(th);
      E->IMAG += B->SRC * sin(th);
    }
  }
}

void idft(Ewalds &ewalds, Bodies &bodies, real R0) {
  real scale = M_PI / R0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    vec<4,real> TRG = 0;
    for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
      real th = 0;
      for( int d=0; d<3; d++ ) th += E->K[d] * B->X[d] * scale;
      real dtmp = E->REAL * sin(th) - E->IMAG * cos(th);
      TRG[0]   += E->REAL * cos(th) + E->IMAG * sin(th);
      for( int d=0; d<3; d++ ) TRG[d+1] -= dtmp * E->K[d];
    }
    B->TRG = TRG;
  }
}
}

template<>
void Kernel<Laplace>::EwaldWave(Bodies &bodies) const {         // Ewald wave part on CPU
  real scale = M_PI / R0;
  real coef = .25 / M_PI / M_PI / SIGMA / R0;
  real coef2 = scale * scale / (4 * ALPHA * ALPHA);

  Ewalds ewalds;
  real kmaxsq = KSIZE * KSIZE;
  int kmax = KSIZE;
  for( int l=0; l<=kmax; l++ ) {
    int mmin = -kmax;
    if( l==0 ) mmin = 0;
    for( int m=mmin; m<=kmax; m++ ) {
      int nmin = -kmax;
      if( l==0 && m==0 ) nmin=1;
      for( int n=nmin; n<=kmax; n++ ) {
        real ksq = l * l + m * m + n * n;
        if( ksq <= kmaxsq ) {
          Ewald ewald;
          ewald.K[0] = l;
          ewald.K[1] = m;
          ewald.K[2] = n;
          ewalds.push_back(ewald);
        }
      }
    }
  }

  dft(ewalds,bodies,R0);
  for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
    real R2 = norm(E->K);
    real factor = coef * exp(-R2 * coef2) / R2;
    E->REAL *= factor;
    E->IMAG *= factor;
  }
  idft(ewalds,bodies,R0);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int d=0; d<3; d++ ) B->TRG[d+1] *= scale;
  }

  vect dipole = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    dipole += (B->X - R0) * B->SRC;
  }
  coef = M_PI / (6 * R0 * R0 * R0);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->TRG[0] += coef * norm(dipole) / bodies.size() / B->SRC;
    for( int d=0; d!=3; ++d ) {
      B->TRG[d+1] += coef * dipole[d];
    }
  }
}
