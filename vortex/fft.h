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
#ifndef fft_h
#define fft_h
#include <fftw3.h>
#include "parallelfmm.h"

class FastFourierTransform : public ParallelFMM {
private:
  int nxLocal;
  int numSend;
  int   *Kk;
  int   *NkSend;
  int   *NkRecv;
  float *EkSend;
  float *EkRecv;
  fftw_complex *vec1d;
  fftw_complex *vec2d;
  fftw_plan forward1d;
  fftw_plan forward2d;
  fftw_plan backward1d;
  fftw_plan backward2d;

protected:
  float *realSend;
  float *realRecv;
  float *imagSend;
  float *imagRecv;

public:
  const int nx;
  int numBodies;
  int numGlobal;

public:
  FastFourierTransform(int N) : nx(N) {
    nxLocal    = nx / MPISIZE;
    numGlobal  = nx * nx * nx;
    numBodies  = nx * nx * nxLocal;
    numSend    = nx * nxLocal * nxLocal;
    Kk         = new int   [nx];
    NkSend     = new int   [nx];
    NkRecv     = new int   [nx];
    EkSend     = new float [nx];
    EkRecv     = new float [nx];
    realSend   = new float [numBodies];
    realRecv   = new float [numBodies];
    imagSend   = new float [numBodies];
    imagRecv   = new float [numBodies];
    vec1d      = (fftw_complex*) fftw_malloc(nx *      sizeof(fftw_complex));
    vec2d      = (fftw_complex*) fftw_malloc(nx * nx * sizeof(fftw_complex));
    forward1d  = fftw_plan_dft_1d(nx,     vec1d, vec1d, FFTW_FORWARD,  FFTW_ESTIMATE);
    forward2d  = fftw_plan_dft_2d(nx, nx, vec2d, vec2d, FFTW_FORWARD,  FFTW_ESTIMATE);
    backward1d = fftw_plan_dft_1d(nx,     vec1d, vec1d, FFTW_BACKWARD, FFTW_ESTIMATE);
    backward2d = fftw_plan_dft_2d(nx, nx, vec2d, vec2d, FFTW_BACKWARD, FFTW_ESTIMATE);
    for( int k=0; k!=nx/2; ++k ) {
      Kk[k] = k;
      Kk[k+nx/2] = k - nx/2;
    }
  }

  ~FastFourierTransform() {
    fftw_destroy_plan(forward1d);
    fftw_destroy_plan(forward2d);
    fftw_destroy_plan(backward1d);
    fftw_destroy_plan(backward2d);
    fftw_free(vec1d);
    fftw_free(vec2d);
    delete[] Kk;
    delete[] NkSend;
    delete[] NkRecv;
    delete[] EkSend;
    delete[] EkRecv;
    delete[] realSend;
    delete[] realRecv;
    delete[] imagSend;
    delete[] imagRecv;
  }

  void forwardFFT() {
    for( int iz=0; iz<nxLocal; ++iz ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int i = iz * nx * nx + iy * nx + ix;
          vec2d[ix+iy*nx][0] = realRecv[i] / nx / nx;
          vec2d[ix+iy*nx][1] = 0;
        }
      }
      fftw_execute(forward2d);
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int i = ix * nx * nxLocal + iy * nxLocal + iz;
          realSend[i] = vec2d[ix+iy*nx][0];
          imagSend[i] = vec2d[ix+iy*nx][1];
        }
      }
    }
    MPI_Alltoall(realSend,numSend,MPI_FLOAT,realRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Alltoall(imagSend,numSend,MPI_FLOAT,imagRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    for( int iz=0; iz<nxLocal; ++iz ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int iix = ix % nxLocal;
          int iiz = iz + (ix / nxLocal) * nxLocal;
          int i = iiz * nx * nxLocal + iy * nxLocal + iix;
          vec1d[ix][0] = realRecv[i] / nx;
          vec1d[ix][1] = imagRecv[i] / nx;
        }
        fftw_execute(forward1d);
        for( int ix=0; ix<nx; ++ix ) {
          int i = iz * nx * nx + iy * nx + ix;
          realSend[i] = vec1d[ix][0];
          imagSend[i] = vec1d[ix][1];
        }
      }
    }
  }

  void backwardFFT() {
    for( int iz=0; iz<nxLocal; ++iz ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int i = iz * nx * nx + iy * nx + ix;
          vec1d[ix][0] = realRecv[i];
          vec1d[ix][1] = imagRecv[i];
        }
        fftw_execute(backward1d);
        for( int ix=0; ix<nx; ++ix ) {
          int i = ix * nx * nxLocal + iy * nxLocal + iz;
          realSend[i] = vec1d[ix][0];
          imagSend[i] = vec1d[ix][1];
        }
      }
    }
    MPI_Alltoall(realSend,numSend,MPI_FLOAT,realRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Alltoall(imagSend,numSend,MPI_FLOAT,imagRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    for( int iz=0; iz<nxLocal; ++iz ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int iix = ix % nxLocal;
          int iiz = iz + (ix / nxLocal) * nxLocal;
          int i = iiz * nx * nxLocal + iy * nxLocal + iix;
          vec2d[ix+iy*nx][0] = realRecv[i];
          vec2d[ix+iy*nx][1] = imagRecv[i];
        }
      }
      fftw_execute(backward2d);
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int i = iz * nx * nx + iy * nx + ix;
          realSend[i] = vec2d[ix+iy*nx][0];
          imagSend[i] = vec2d[ix+iy*nx][1];
        }
      }
    }
  }

  void xDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      realRecv[i] = -Kk[ix] * imagSend[i];
      imagRecv[i] =  Kk[ix] * realSend[i];
    }
    backwardFFT();
  }

  void yDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      realRecv[i] = -Kk[iy] * imagSend[i];
      imagRecv[i] =  Kk[iy] * realSend[i];
    }
    backwardFFT();
  }

  void zDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int iz = (i + numBodies * MPIRANK) % nx;
      realRecv[i] = -Kk[iz] * imagSend[i];
      imagRecv[i] =  Kk[iz] * realSend[i];
    }
    backwardFFT();
  }

  void initSpectrum() {
    for( int k=0; k<nx; ++k ) {
      EkSend[k] = 0;
      NkSend[k] = 0;
    }
  }

#if 0
  void addSpectrum() {
    for( int iz=0; iz<nxLocal; ++iz ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int ix=0; ix<nx; ++ix ) {
          int iiz = iz + nxLocal * MPIRANK;
          int i = iz * nx * nx + iy * nx + ix;
          int k = floor(sqrtf(Kk[iiz] * Kk[iiz] + Kk[iy] * Kk[iy] + Kk[ix] * Kk[ix]));
          EkSend[k] += (realSend[i] * realSend[i] + imagSend[i] * imagSend[i]);
        }
      }
    }
  }
#else
  void addSpectrum() {
    for( int k=0; k<nx; ++k ) NkSend[k] = 0;
    for( int ix=1; ix<nxLocal; ++ix ) {
      for( int iy=1; iy<=nx/2; ++iy ) {
        for( int iz=1; iz<=nx/2; ++iz ) {
          int iix = ix + nxLocal * MPIRANK;
          if( iix <= nx/2 ) {
            int i = ix * nx * nx + iy * nx + iz;
            float kf = sqrtf(iix * iix + iy * iy + iz * iz);
            int k = floor(kf);
            EkSend[k] += (realSend[i] * realSend[i] + imagSend[i] * imagSend[i]) * 4 * M_PI * kf * kf;
            NkSend[k]++;
          }
        }
      }
    }
  }
#endif

  void writeSpectrum() {
    MPI_Reduce(NkSend,NkRecv,nx,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(EkSend,EkRecv,nx,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    if( MPIRANK == 0 ) {
      std::ofstream fid("statistics.dat",std::ios::out | std::ios::app);
      for( int k=0; k<nx; ++k ) {
        if( NkRecv[k] == 0 ) NkRecv[k] = 1;
        fid << EkRecv[k]/NkRecv[k] << std::endl;
      }
      fid.close();
    }
  }

};

#endif
