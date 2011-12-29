#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include "serialfmm.h"

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void ewaldreal(Bodies &bodies, real rscale, real xmax) {
  for( B_iter Bi=bodies.begin(); Bi!=bodies.end(); ++Bi ) {
    for( B_iter Bj=bodies.begin(); Bj!=bodies.end(); ++Bj ) {
      vect dX = Bi->X - Bj->X;
      for( int d=0; d<3; d++ ) {
        if( dX[d] < -xmax/2 ) {
          dX[d] += xmax;
        }
        if( dX[d] >= xmax/2 ) {
          dX[d] -= xmax;
        }
      }
      real R2 = norm(dX);
      if( R2 != 0 ) {
        real R2s = R2 * rscale * rscale;
        real Rs = std::sqrt(R2s);
        real invRs = 1 / Rs;
        real invR2s = invRs * invRs;
        real invR3s = invR2s * invRs;
        real dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
        dtmp *= rscale * rscale * rscale;
        Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * rscale;
        for( int d=0; d<3; d++ ) {
          Bi->TRG[d+1] -= dX[d] * dtmp;
        }
      }
    }
  }
}

void MR3calcewald_dft_host(Ewalds &ewalds, Bodies &bodies, real xmax) {
  real scale = 2 * M_PI / xmax;
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

void MR3calcewald_idft_host(Ewalds &ewalds, Bodies &bodies, real xmax) {
  real scale = 2 * M_PI / xmax;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    vec<4,real> TRG = 0;
    for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
      real th = 0;
      for( int d=0; d<3; d++ ) th += E->K[d] * B->X[d] * scale;
      real dtmp = E->REAL * sin(th) - E->IMAG * cos(th);
      TRG[0]    += E->REAL * cos(th) + E->IMAG * sin(th);
      for( int d=0; d<3; d++ ) {
        TRG[d+1] -= dtmp * E->K[d];
      }
    }
    B->TRG = TRG;
  }
}

void MR3calcewald_host(Bodies &bodies, real ksize, real alpha, real epsilon, real xmax) {
  real scale = 2 * M_PI / xmax;
  real coef = .5 / M_PI / M_PI / epsilon / xmax;
  real coef2 = scale * scale / (4 * alpha * alpha);

  Ewalds ewalds;
  real kmaxsq = ksize * ksize;
  int kmax = ksize;
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

  MR3calcewald_dft_host(ewalds,bodies,xmax);
  for( E_iter E=ewalds.begin(); E!=ewalds.end(); ++E ) {
    real R2 = norm(E->K);
    real factor = coef * exp(-R2 * coef2) / R2;
    E->REAL *= factor;
    E->IMAG *= factor;
  }

  MR3calcewald_idft_host(ewalds,bodies,xmax);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int d=0; d<3; d++ ) B->TRG[d+1] *= scale;
  }
}

int main() {
  const int numBodies = 1000;
  const real xmax = 100.0;
  const real ksize = 11.0;
  const real alpha = 0.1;
  IMAGES = 2;
  Bodies bodies(numBodies);
  Bodies jbodies;
  SerialFMM<Laplace> FMM;

  srand48(2);
  // set positions and types
  real average = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int d=0; d!=3; ++d ) {
      B->X[d] = drand48() * xmax;
    }
    B->SRC = drand48();
    average += B->SRC;
    B->TRG = 0;
  }
  average /= numBodies;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->SRC -= average;
  }

  double tic = get_time();
  // Ewald wave part
  MR3calcewald_host(bodies,ksize,alpha,1.0/(M_PI*4.0),xmax);
  // Ewald real part
  ewaldreal(bodies,alpha,xmax);
  // Potential energy self term
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->TRG[0] -= M_2_SQRTPI * B->SRC * alpha;
  }
  // Dipole correction
  vect dipole = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    dipole += (B->X - 0.5 * xmax) * B->SRC;
  }
  real coef = 4 * M_PI / (3 * xmax * xmax * xmax);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->TRG[0] += coef * norm(dipole) / numBodies / B->SRC;
    for( int d=0; d!=3; ++d ) {
      B->TRG[d+1] += coef * dipole[d];
    }
  }
  double toc = get_time();
  std::cout << "ewald  : " << toc-tic << std::endl;

  tic = get_time();
  // Direct with images
  FMM.setDomain(bodies,xmax/2,xmax/2);
  jbodies = FMM.periodicBodies(bodies);
  FMM.buffer = bodies;
  FMM.initTarget(FMM.buffer);
  FMM.evalP2P(FMM.buffer,jbodies);
  toc = get_time();
  std::cout << "direct : " << toc-tic << std::endl;

  // Calculate error
  real p = 0, p2 = 0, diff = 0, norm = 0;
  B_iter B2 = FMM.buffer.begin();
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {
    p += B->TRG[0] * B->SRC;
    p2 += B2->TRG[0] * B2->SRC;
    for( int d=0; d!=3; ++d ) {
      diff += (B->TRG[d+1] - B2->TRG[d+1]) * (B->TRG[d+1] - B2->TRG[d+1]);
      norm += B2->TRG[d+1] * B2->TRG[d+1];
    }
  }
  std::cout << "p err  : " << std::abs((p-p2)/p2) << std::endl;
  std::cout << "f err  : " << std::sqrt(diff/norm) << std::endl;
}
