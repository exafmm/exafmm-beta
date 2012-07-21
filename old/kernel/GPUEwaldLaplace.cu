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

__global__ void dft(gpureal *ewaldsGlob, gpureal *bodiesGlob, const int numEwalds, const int numBodies, const real R0) {
  int i = blockIdx.x * THREADS + threadIdx.x;
  gpureal scale = M_PI / R0;
  gpureal REAL = 0, IMAG = 0;
  __shared__ gpureal bodiesShrd[4*THREADS];
  for( int iblok=0; iblok<(numBodies-1)/THREADS; ++iblok ) {
    int ibodies = iblok * THREADS + threadIdx.x;
    __syncthreads();
    bodiesShrd[4*threadIdx.x+0] = bodiesGlob[8*ibodies+0];
    bodiesShrd[4*threadIdx.x+1] = bodiesGlob[8*ibodies+1];
    bodiesShrd[4*threadIdx.x+2] = bodiesGlob[8*ibodies+2];
    bodiesShrd[4*threadIdx.x+3] = bodiesGlob[8*ibodies+3];
    __syncthreads();
    for( int j=0; j<THREADS; ++j ) {
      gpureal th = 0;
      for( int d=0; d<3; d++ ) th += ewaldsGlob[5*i+d] * bodiesShrd[4*j+d] * scale;
      REAL += bodiesShrd[4*j+3] * cosf(th);
      IMAG += bodiesShrd[4*j+3] * sinf(th);
    }
  }
  int iblok = (numBodies-1)/THREADS;
  int ibodies = iblok * THREADS + threadIdx.x;
  __syncthreads();
  if( threadIdx.x < numBodies - iblok * THREADS ) {
    bodiesShrd[4*threadIdx.x+0] = bodiesGlob[8*ibodies+0];
    bodiesShrd[4*threadIdx.x+1] = bodiesGlob[8*ibodies+1];
    bodiesShrd[4*threadIdx.x+2] = bodiesGlob[8*ibodies+2];
    bodiesShrd[4*threadIdx.x+3] = bodiesGlob[8*ibodies+3];
  }
  __syncthreads();
  for( int j=0; j<numBodies-iblok*THREADS; ++j ) {
    gpureal th = 0;
    for( int d=0; d<3; d++ ) th += ewaldsGlob[5*i+d] * bodiesShrd[4*j+d] * scale;
    REAL += bodiesShrd[4*j+3] * cosf(th);
    IMAG += bodiesShrd[4*j+3] * sinf(th);
  }
  ewaldsGlob[5*i+3] = REAL;
  ewaldsGlob[5*i+4] = IMAG;
}

__global__ void idft(gpureal *ewaldsGlob, gpureal *bodiesGlob, const int numEwalds, const int numBodies, const real R0) {
  int i = blockIdx.x * THREADS + threadIdx.x;
  gpureal scale = M_PI / R0;
  gpureal TRG[4] = {0,0,0,0};
  __shared__ gpureal ewaldsShrd[5*THREADS];
  for( int iblok=0; iblok<(numEwalds-1)/THREADS; ++iblok ) {
    int iewalds = iblok * THREADS + threadIdx.x;
    __syncthreads();
    ewaldsShrd[5*threadIdx.x+0] = ewaldsGlob[5*iewalds+0];
    ewaldsShrd[5*threadIdx.x+1] = ewaldsGlob[5*iewalds+1];
    ewaldsShrd[5*threadIdx.x+2] = ewaldsGlob[5*iewalds+2];
    ewaldsShrd[5*threadIdx.x+3] = ewaldsGlob[5*iewalds+3];
    ewaldsShrd[5*threadIdx.x+4] = ewaldsGlob[5*iewalds+4];
    __syncthreads();
    for( int j=0; j<THREADS; ++j ) {
      gpureal th = 0;
      for( int d=0; d<3; d++ ) th += ewaldsShrd[5*j+d] * bodiesGlob[8*i+d] * scale;
      gpureal ftmp = ewaldsShrd[5*j+3] * sinf(th) - ewaldsShrd[5*j+4] * cosf(th);
      TRG[0]      += ewaldsShrd[5*j+3] * cosf(th) + ewaldsShrd[5*j+4] * sinf(th);
      for( int d=0; d<3; d++ ) TRG[d+1] -= ftmp * ewaldsShrd[5*j+d];
    }
  }
  int iblok = (numEwalds-1)/THREADS;
  int iewalds = iblok * THREADS + threadIdx.x;
  __syncthreads();
  if( threadIdx.x < numEwalds - iblok * THREADS ) {
    ewaldsShrd[5*threadIdx.x+0] = ewaldsGlob[5*iewalds+0];
    ewaldsShrd[5*threadIdx.x+1] = ewaldsGlob[5*iewalds+1];
    ewaldsShrd[5*threadIdx.x+2] = ewaldsGlob[5*iewalds+2];
    ewaldsShrd[5*threadIdx.x+3] = ewaldsGlob[5*iewalds+3];
    ewaldsShrd[5*threadIdx.x+4] = ewaldsGlob[5*iewalds+4];
  }
  __syncthreads();
  for( int j=0; j<numEwalds-iblok*THREADS; ++j ) {
    gpureal th = 0;
    for( int d=0; d<3; d++ ) th += ewaldsShrd[5*j+d] * bodiesGlob[8*i+d] * scale;
    gpureal ftmp = ewaldsShrd[5*j+3] * sinf(th) - ewaldsShrd[5*j+4] * cosf(th);
    TRG[0]      += ewaldsShrd[5*j+3] * cosf(th) + ewaldsShrd[5*j+4] * sinf(th);
    for( int d=0; d<3; d++ ) TRG[d+1] -= ftmp * ewaldsShrd[5*j+d];
  }
  for( int d=0; d<4; d++ ) bodiesGlob[8*i+d+4] = TRG[d];
}

__global__ void factor(gpureal *ewaldsGlob, const gpureal coef, const gpureal coef2) {
  int i = blockIdx.x * THREADS + threadIdx.x;
  real R2 = ewaldsGlob[5*i+0] * ewaldsGlob[5*i+0]
          + ewaldsGlob[5*i+1] * ewaldsGlob[5*i+1]
          + ewaldsGlob[5*i+2] * ewaldsGlob[5*i+2];
  real factor = coef * expf(-R2 * coef2) / R2;
  ewaldsGlob[5*i+3] *= factor;
  ewaldsGlob[5*i+4] *= factor;
}

template<>
void Kernel<Laplace>::EwaldWave(Bodies &bodies) const {     // Ewald wave part on CPU
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

  gpureal *bodiesDevc;
  gpureal *ewaldsDevc;
  int numBodies = bodies.size();
  int numEwalds = ewalds.size();
  int paddedBodies = ((numBodies-1) / THREADS + 1) * THREADS;
  int paddedEwalds = ((numEwalds-1) / THREADS + 1) * THREADS;
  gpureal *bodiesHost = new gpureal [8*paddedBodies];
  gpureal *ewaldsHost = new gpureal [5*paddedEwalds];
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    bodiesHost[8*i+0] = B->X[0];
    bodiesHost[8*i+1] = B->X[1];
    bodiesHost[8*i+2] = B->X[2];
    bodiesHost[8*i+3] = B->SRC;
    bodiesHost[8*i+4] = B->TRG[0];
    bodiesHost[8*i+5] = B->TRG[1];
    bodiesHost[8*i+6] = B->TRG[2];
    bodiesHost[8*i+7] = B->TRG[3];
  }
  for( int i=numBodies; i!=paddedBodies; ++i ) {
    for( int d=0; d!=8; ++d ) bodiesHost[8*i+d] = 0;
  }
  for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
    int i = E-ewalds.begin();
    ewaldsHost[5*i+0] = E->K[0];
    ewaldsHost[5*i+1] = E->K[1];
    ewaldsHost[5*i+2] = E->K[2];
    ewaldsHost[5*i+3] = E->REAL;
    ewaldsHost[5*i+4] = E->IMAG;
  }
  for( int i=numEwalds; i!=paddedEwalds; ++i ) {
    for( int d=0; d!=5; ++d ) ewaldsHost[5*i+d] = 0;
  }
  cudaMalloc((void**)&bodiesDevc,sizeof(gpureal)*8*paddedBodies);
  cudaMalloc((void**)&ewaldsDevc,sizeof(gpureal)*5*paddedEwalds);
  cudaMemcpy(bodiesDevc,bodiesHost,sizeof(gpureal)*8*paddedBodies,cudaMemcpyHostToDevice);
  cudaMemcpy(ewaldsDevc,ewaldsHost,sizeof(gpureal)*5*paddedEwalds,cudaMemcpyHostToDevice);
  int BLOCKS = (numEwalds - 1) / THREADS + 1;
  dft<<<BLOCKS,THREADS>>>(ewaldsDevc,bodiesDevc,numEwalds,numBodies,R0);
  factor<<<BLOCKS,THREADS>>>(ewaldsDevc,coef,coef2);
  BLOCKS = (numBodies - 1) / THREADS + 1;
  idft<<<BLOCKS,THREADS>>>(ewaldsDevc,bodiesDevc,numEwalds,paddedBodies,R0);
  cudaMemcpy(bodiesHost,bodiesDevc,sizeof(gpureal)*8*paddedBodies,cudaMemcpyDeviceToHost);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->TRG[0] = bodiesHost[8*i+4];
    B->TRG[1] = bodiesHost[8*i+5];
    B->TRG[2] = bodiesHost[8*i+6];
    B->TRG[3] = bodiesHost[8*i+7];
  }
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
