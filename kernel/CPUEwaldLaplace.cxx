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
#define KERNEL
#include "kernel.h"
#undef KERNEL

namespace {
void dft(Ewalds &ewalds, Bodies &bodies, real R0) {
  real scale = M_PI / R0;
//#pragma omp parallel for
  for( int i=0; i<int(ewalds.size()); ++i ) {
    E_iter E = ewalds.begin() + i;
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
//#pragma omp parallel for
  for( int i=0; i<int(bodies.size()); ++i ) {
    B_iter B = bodies.begin() + i;
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
void Kernel<Laplace>::EwaldReal(C_iter Ci, C_iter Cj) const {   // Ewald real part on CPU
  for( int i=0; i<Ci->NDLEAF; ++i ) {                           // Loop over target bodies
    B_iter Bi = Ci->LEAF + i;                                   //  Target body iterator
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      for( int d=0; d<3; d++ ) {                                //   Loop over dimensions
        if( dist[d] < -R0 ) dist[d] += 2 * R0;                  //    Wrap domain so that target is always at
        if( dist[d] >= R0 ) dist[d] -= 2 * R0;                  //    the center of a [-R0,R0]^3 source cube
      }                                                         //   End loop over dimensions
      real R2 = norm(dist);                                     //   R^2
      if( R2 != 0 ) {                                           //   Exclude self interaction
        real R2s = R2 * ALPHA * ALPHA;                          //    (R * alpha)^2
        real Rs = std::sqrt(R2s);                               //    R * alpha
        real invRs = 1 / Rs;                                    //    1 / (R * alpha)
        real invR2s = invRs * invRs;                            //    1 / (R * alpha)^2
        real invR3s = invR2s * invRs;                           //    1 / (R * alpha)^3
        real dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
        dtmp *= ALPHA * ALPHA * ALPHA;                          //    Scale temporary value
        Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * ALPHA;       //    Ewald real potential
        Bi->TRG[1] -= dist[0] * dtmp;                           //    x component of Ewald real force
        Bi->TRG[2] -= dist[1] * dtmp;                           //    y component of Ewald real force
        Bi->TRG[3] -= dist[2] * dtmp;                           //    z component of Ewald real force
      }                                                         //   End if for self interaction
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

template<>
void Kernel<Laplace>::EwaldWave(Bodies &bodies) const {         // Ewald wave part on CPU
  real scale = M_PI / R0;
  real coef = .25 / M_PI / M_PI / SIGMA / R0;
  real coef2 = scale * scale / (4 * ALPHA * ALPHA);

  Ewalds ewalds;
  real kmaxsq = KSIZE * KSIZE;
  int kmax = int(KSIZE);
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
          ewald.REAL = ewald.IMAG = 0;
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

#if 0
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
#endif
}
