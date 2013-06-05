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

namespace{
//! Get r,theta,phi from x,y,z
void cart2sph(real& r, real& theta, real& phi, vect dist) {
  r = sqrt(norm(dist)) * (1 + EPS);                           // r = sqrt(x^2 + y^2 + z^2)
  if( r < EPS ) {                                             // If r == 0
    theta = 0;                                                //  theta can be anything so we set it to 0
  } else {                                                    // If r != 0
    theta = acos(dist[2] / r);                                //  theta = acos(z / r)
  }                                                           // End if for r == 0
  if( fabs(dist[0]) + fabs(dist[1]) < EPS ) {                 // If |x| < eps & |y| < eps
    phi = 0;                                                  //  phi can be anything so we set it to 0
  } else if( fabs(dist[0]) < EPS ) {                          // If |x| < eps
    phi = dist[1] / fabs(dist[1]) * M_PI * 0.5;               //  phi = sign(y) * pi / 2
  } else if( dist[0] > 0 ) {                                  // If x > 0
    phi = atan(dist[1] / dist[0]);                            //  phi = atan(y / x)
  } else {                                                    // If x < 0
    phi = atan(dist[1] / dist[0]) + M_PI;                     //  phi = atan(y / x) + pi
  }                                                           // End if for x,y cases
}

//! Spherical to cartesian coordinates
template<typename T>
void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) {
  cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
               + cos(theta) * cos(phi) / r * spherical[1]
               - sin(phi) / r / sin(theta) * spherical[2];
  cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
               + cos(theta) * sin(phi) / r * spherical[1]
               + cos(phi) / r / sin(theta) * spherical[2];
  cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
               - sin(theta) / r * spherical[1];
}

}

template<>
void Kernel<Laplace>::initialize() {}

template<>
void Kernel<Laplace>::P2M(C_iter Cj) {
  real Rmax = 0;
  complex Ynm[4*P*P], YnmTheta[4*P*P];
  for( B_iter B=Cj->LEAF; B!=Cj->LEAF+Cj->NDLEAF; ++B ) {
    vect dist = B->X - Cj->X;
    real R = std::sqrt(norm(dist));
    if( R > Rmax ) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        int nm  = n * n + n + m;
        int nms = n * (n + 1) / 2 + m;
        Cj->M[nms] += B->SRC * Ynm[nm];
      }
    }
  }
  Cj->RMAX = Rmax;
  Cj->RCRIT = std::min(Cj->R,Rmax);
}

template<>
void Kernel<Laplace>::M2M(C_iter Ci) {
  const complex I(0.,1.);
  complex Ynm[4*P*P], YnmTheta[4*P*P];
  real Rmax = Ci->RMAX;
  for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {
    vect dist = Ci->X - Cj->X;
    real R = std::sqrt(norm(dist)) + Cj->RCRIT;
    if( R > Rmax ) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        complex M = 0;
        for( int n=0; n<=j; ++n ) {
          for( int m=-n; m<=std::min(k-1,n); ++m ) {
            if( j-n >= k-m ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n + m;
              M += Cj->M[jnkms] * std::pow(I,real(m-abs(m))) * Ynm[nm]
                 * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          for( int m=k; m<=n; ++m ) {
            if( j-n >= m-k ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n + m;
              M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                 * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        Ci->M[jks] += M * EPS;
      }
    }
  }
  Ci->RMAX = Rmax;
  Ci->RCRIT = std::min(Ci->R,Rmax);
}

template<>
void Kernel<Laplace>::M2L(C_iter Ci, C_iter Cj) const {
  complex Ynm[4*P*P], YnmTheta[4*P*P];
  vect dist = Ci->X - Cj->X - Xperiodic;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalLocal(rho,alpha,beta,Ynm,YnmTheta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex L = 0;
      for( int n=0; n!=P; ++n ) {
        for( int m=-n; m<0; ++m ) {
          int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 - m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynm[jnkm];
        }
        for( int m=0; m<=n; ++m ) {
          int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 + m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += Cj->M[nms] * Cnm[jknm] * Ynm[jnkm];
        }
      }
      Ci->L[jks] += L;
    }
  }
}

template<>
void Kernel<Laplace>::M2P(C_iter Ci, C_iter Cj) const {
  const complex I(0.,1.);                                       // Imaginary unit
  complex Ynm[4*P*P], YnmTheta[4*P*P];
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Cj->X - Xperiodic;
    vect spherical = 0;
    vect cartesian = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->TRG[0] += std::real(Cj->M[nms] * Ynm[nm]);
      spherical[0] -= std::real(Cj->M[nms] * Ynm[nm]) / r * (n+1);
      spherical[1] += std::real(Cj->M[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->TRG[0] += 2 * std::real(Cj->M[nms] * Ynm[nm]);
        spherical[0] -= 2 * std::real(Cj->M[nms] *Ynm[nm]) / r * (n+1);
        spherical[1] += 2 * std::real(Cj->M[nms] *YnmTheta[nm]);
        spherical[2] += 2 * std::real(Cj->M[nms] *Ynm[nm] * I) * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    B->TRG[1] += cartesian[0];
    B->TRG[2] += cartesian[1];
    B->TRG[3] += cartesian[2];
  }
}

template<>
void Kernel<Laplace>::L2L(C_iter Ci) const {
  const complex I(0.,1.);
  complex Ynm[P*P], YnmTheta[P*P];
  C_iter Cj = Ci0 + Ci->PARENT;
  vect dist = Ci->X - Cj->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex L = 0;
      for( int n=j; n!=P; ++n ) {
        for( int m=j+k-n; m<0; ++m ) {
          int jnkm = (n - j) * (n - j) + n - j + m - k;
          int nm   = n * n + n - m;
          int nms  = n * (n + 1) / 2 - m;
          L += std::conj(Cj->L[nms]) * Ynm[jnkm]
             * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nm   = n * n + n + m;
            int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[nms] * std::pow(I,real(m-k-abs(m-k)))
               * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
          }
        }
      }
      Ci->L[jks] += L * EPS;
    }
  }
}

template<>
void Kernel<Laplace>::L2P(C_iter Ci) const {
  const complex I(0.,1.);                                       // Imaginary unit
  complex Ynm[4*P*P], YnmTheta[4*P*P];
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Ci->X;
    vect spherical = 0;
    vect cartesian = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
      spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
      spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->TRG[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]);
        spherical[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) / r * n;
        spherical[1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
        spherical[2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;
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

#include "../kernel/CPUEwaldLaplace.cxx"
