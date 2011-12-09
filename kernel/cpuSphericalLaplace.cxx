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
#include "laplace.h"

template<>
void Kernel<Laplace>::initialize() {}

template<>
void Kernel<Laplace>::P2M(C_iter Ci) const {
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Ci->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        const int nm  = n * n + n + m;
        const int nms = n * (n + 1) / 2 + m;
        Ci->M[3*nms] += double(B->SRC[0]) * Ynm[nm];
      }
    }
  }
}

template<>
void Kernel<Laplace>::M2M(C_iter Ci, C_iter Cj) const {
  const complex I(0.,1.);                                       // Imaginary unit
  vect dist = Ci->X - Cj->X;
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
            M += Cj->M[3*jnkms] * std::pow(I,double(m-abs(m))) * Ynm[nm]
               * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            M += std::conj(Cj->M[3*jnkms]) * Ynm[nm]
               * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
          }
        }
      }
      Ci->M[3*jks] += M;
    }
  }
}

template<>
void Kernel<Laplace>::M2L(C_iter Ci, C_iter Cj) const {
  vect dist = Ci->X - Cj->X - Xperiodic;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalLocal(rho,alpha,beta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      const int jk = j * j + j + k;
      const int jks = j * (j + 1) / 2 + k;
      complex L = 0;
      for( int n=0; n!=P; ++n ) {
        for( int m=-n; m<0; ++m ) {
          const int nm   = n * n + n + m;
          const int nms  = n * (n + 1) / 2 - m;
          const int jknm = jk * P2 + nm;
          const int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += std::conj(Cj->M[3*nms]) * Cnm[jknm] * Ynm[jnkm];
        }
        for( int m=0; m<=n; ++m ) {
          const int nm   = n * n + n + m;
          const int nms  = n * (n + 1) / 2 + m;
          const int jknm = jk * P2 + nm;
          const int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += Cj->M[3*nms] * Cnm[jknm] * Ynm[jnkm];
        }
      }
      Ci->L[3*jks] += L;
    }
  }
}

template<>
void Kernel<Laplace>::M2P(C_iter Ci, C_iter Cj) const {
  const complex I(0.,1.);                                       // Imaginary unit
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Cj->X - Xperiodic;
    vect spherical = 0;
    vect cartesian = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->TRG[0] += (Cj->M[3*nms] * Ynm[nm]).real();
      spherical[0] -= (Cj->M[3*nms] * Ynm[nm]).real() / r * (n+1);
      spherical[1] += (Cj->M[3*nms] * YnmTheta[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->TRG[0] += 2 * (Cj->M[3*nms] * Ynm[nm]).real();
        spherical[0] -= 2 * (Cj->M[3*nms] *Ynm[nm]).real() / r * (n+1);
        spherical[1] += 2 * (Cj->M[3*nms] *YnmTheta[nm]).real();
        spherical[2] += 2 * (Cj->M[3*nms] *Ynm[nm] * I).real() * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    B->TRG[1] += cartesian[0];
    B->TRG[2] += cartesian[1];
    B->TRG[3] += cartesian[2];
  }
}

template<>
void Kernel<Laplace>::L2L(C_iter Ci, C_iter Cj) const {
  const complex I(0.,1.);                                       // Imaginary unit
  vect dist = Ci->X - Cj->X;
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
          L += std::conj(Cj->L[3*nms]) * Ynm[jnkm]
             * double(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[3*nms] * std::pow(I,double(m-k-abs(m-k)))
               * Ynm[jnkm] * double(Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
        }
      }
      Ci->L[3*jks] += L;
    }
  }
}

template<>
void Kernel<Laplace>::L2P(C_iter Ci) const {
  const complex I(0.,1.);                                       // Imaginary unit
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Ci->X;
    vect spherical = 0;
    vect cartesian = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->TRG[0] += (Ci->L[3*nms] * Ynm[nm]).real();
      spherical[0] += (Ci->L[3*nms] * Ynm[nm]).real() / r * n;
      spherical[1] += (Ci->L[3*nms] * YnmTheta[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->TRG[0] += 2 * (Ci->L[3*nms] * Ynm[nm]).real();
        spherical[0] += 2 * (Ci->L[3*nms] * Ynm[nm]).real() / r * n;
        spherical[1] += 2 * (Ci->L[3*nms] * YnmTheta[nm]).real();
        spherical[2] += 2 * (Ci->L[3*nms] * Ynm[nm] * I).real() * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    B->TRG[1] += cartesian[0];
    B->TRG[2] += cartesian[1];
    B->TRG[3] += cartesian[2];
  }
}

template<>
void Kernel<Laplace>::finalize() {}
