#ifndef fft_h
#define fft_h
#include <fftw3.h>
#include "let.h"

class FastFourierTransform : public LocalEssentialTree {
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

public:
  FastFourierTransform(int N) : nx(N) {
    nxLocal    = nx / SIZE;
    numBodies  = nx * nx * nxLocal;
    numSend    = nx * nxLocal * nxLocal;
    Kk         = new int   [nx];
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
    delete[] EkSend;
    delete[] EkRecv;
    delete[] realSend;
    delete[] realRecv;
    delete[] imagSend;
    delete[] imagRecv;
  }

  void forwardFFT() {
    for( int ix=0; ix<nxLocal; ++ix ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int i = ix * nx * nx + iy * nx + iz;
          vec2d[iz+iy*nx][0] = realRecv[i] / nx / nx;
          vec2d[iz+iy*nx][1] = 0;
        }
      }
      fftw_execute(forward2d);
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int i = iz * nx * nxLocal + iy * nxLocal + ix;
          realSend[i] = vec2d[iz+iy*nx][0];
          imagSend[i] = vec2d[iz+iy*nx][1];
        }
      }
    }
    MPI_Alltoall(realSend,numSend,MPI_FLOAT,realRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Alltoall(imagSend,numSend,MPI_FLOAT,imagRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    for( int ix=0; ix<nxLocal; ++ix ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int iiz = iz % nxLocal;
          int iix = ix + (iz / nxLocal) * nxLocal;
          int i = iix * nx * nxLocal + iy * nxLocal + iiz;
          vec1d[iz][0] = realRecv[i] / nx;
          vec1d[iz][1] = imagRecv[i] / nx;
        }
        fftw_execute(forward1d);
        for( int iz=0; iz<nx; ++iz ) {
          int i = ix * nx * nx + iy * nx + iz;
          realSend[i] = vec1d[iz][0];
          imagSend[i] = vec1d[iz][1];
        }
      }
    }
  }

  void backwardFFT() {
    for( int ix=0; ix<nxLocal; ++ix ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int i = ix * nx * nx + iy * nx + iz;
          vec1d[iz][0] = realRecv[i];
          vec1d[iz][1] = imagRecv[i];
        }
        fftw_execute(backward1d);
        for( int iz=0; iz<nx; ++iz ) {
          int i = iz * nx * nxLocal + iy * nxLocal + ix;
          realSend[i] = vec1d[iz][0];
          imagSend[i] = vec1d[iz][1];
        }
      }
    }
    MPI_Alltoall(realSend,numSend,MPI_FLOAT,realRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Alltoall(imagSend,numSend,MPI_FLOAT,imagRecv,numSend,MPI_FLOAT,MPI_COMM_WORLD);
    for( int ix=0; ix<nxLocal; ++ix ) {
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int iiz = iz % nxLocal;
          int iix = ix + (iz / nxLocal) * nxLocal;
          int i = iix * nx * nxLocal + iy * nxLocal + iiz;
          vec2d[iz+iy*nx][0] = realRecv[i];
          vec2d[iz+iy*nx][1] = imagRecv[i];
        }
      }
      fftw_execute(backward2d);
      for( int iy=0; iy<nx; ++iy ) {
        for( int iz=0; iz<nx; ++iz ) {
          int i = ix * nx * nx + iy * nx + iz;
          realSend[i] = vec2d[iz+iy*nx][0];
          imagSend[i] = vec2d[iz+iy*nx][1];
        }
      }
    }
  }

  void xDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int ix = (i + numBodies * RANK) % nx;
      realRecv[i] = -Kk[ix] * imagSend[i];
      imagRecv[i] =  Kk[ix] * realSend[i];
    }
    backwardFFT();
  }

  void yDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int iy = (i + numBodies * RANK) / nx % nx;
      realRecv[i] = -Kk[iy] * imagSend[i];
      imagRecv[i] =  Kk[iy] * realSend[i];
    }
    backwardFFT();
  }

  void zDerivative() {
    forwardFFT();
    for( int i=0; i!=numBodies; ++i ) {
      int iz = (i + numBodies * RANK) / nx / nx;
      realRecv[i] = -Kk[iz] * imagSend[i];
      imagRecv[i] =  Kk[iz] * realSend[i];
    }
    backwardFFT();
  }

  void initSpectrum() {
    for( int k=0; k<nx; ++k ) {
      EkSend[k] = 0;
    }
  }

  void addSpectrum() {
    for( int ix=0; ix<nxLocal; ++ix ) {
      for( int iy=0; iy<nx/2; ++iy ) {
        for( int iz=0; iz<nx/2; ++iz ) {
          int iix = ix + nxLocal * RANK;
          int i = ix * nx * nx + iy * nx + iz;
          int k = floor(sqrtf(Kk[iix] * Kk[iix] + Kk[iy] * Kk[iy] + Kk[iz] * Kk[iz]));
          EkSend[k] += (realSend[i] * realSend[i] + imagSend[i] * imagSend[i]);
        }
      }
    }
  }

  void writeSpectrum() {
    MPI_Reduce(EkSend,EkRecv,nx,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    if( RANK == 0 ) {
      std::ofstream fid("statistics.dat",std::ios::out | std::ios::app);
      for( int k=0; k<nx; ++k ) {
        fid << EkRecv[k] << std::endl;
      }
      fid.close();
    }
  }

};

#endif
