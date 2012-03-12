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
  r = sqrt(norm(dist))+EPS;                                   // r = sqrt(x^2 + y^2 + z^2) + eps
  theta = acos(dist[2] / r);                                  // theta = acos(z / r)
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

template<int n, int m>
struct Index {
  static const int npm = Index<n,m-1>::npm + 1;
  static const int nmm = Index<n,m-1>::nmm - 1;
  static const int nms = Index<n,m-1>::nms + 1;
};

template<int n>
struct Index<n,0> {
  static const int npm = Index<n-1,n-1>::npm + n + 1;
  static const int nmm = Index<n-1,n-1>::npm + n + 1;
  static const int nms = Index<n-1,n-1>::nms + 1;
};

template<>
struct Index<0,0> {
  static const int npm = 0;
  static const int nmm = 0;
  static const int nms = 0;
};

template<int p, int n=p-1, int m=p-1>
struct Expansion {
  static inline void getYnm(real &x, real &y, real &rho,
                            real &rhom, real &rhon, real &beta, complex &eim,
                            real &Pn, real &P0, real &P1, real &P2,
                            real *prefactor, complex *Ynm) {
    Expansion<p,n-1,m>::getYnm(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
    Ynm[Index<n,m>::npm] = rhon * P0 * prefactor[Index<n,m>::npm] * eim;
    Ynm[Index<n,m>::nmm] = std::conj(Ynm[Index<n,m>::npm]);
    P2 = P1;
    P1 = P0;
    P0 = (x * (2 * n + 1) * P1 - (n + m) * P2) / (n - m + 1);
    rhon *= rho;
  }
  static inline void getYnmTheta(real &x, real &y, real &rho,
                                 real &rhom, real &rhon, real &beta, complex &eim,
                                 real &Pn, real &P0, real &P1, real &P2,
                                 real *prefactor, complex *Ynm, complex *YnmTheta) {
    Expansion<p,n-1,m>::getYnmTheta(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm,YnmTheta);
    Ynm[Index<n,m>::npm] = rhon * P0 * prefactor[Index<n,m>::npm] * eim;
    Ynm[Index<n,m>::nmm] = std::conj(Ynm[Index<n,m>::npm]);
    P2 = P1;
    P1 = P0;
    P0 = (x * (2 * n + 1) * P1 - (n + m) * P2) / (n - m + 1);
    YnmTheta[Index<n,m>::npm] = rhon * ((n - m + 1) * P0 - (n + 1) * x * P1) / y * prefactor[Index<n,m>::npm] * eim;
    rhon *= rho;
  }
};

template<int p, int m>
struct Expansion<p,m,m> {
  static inline void getYnm(real &x, real &y, real &rho,
                            real &rhom, real &rhon, real &beta, complex &eim,
                            real &Pn, real &P0, real &P1, real &P2,
                            real *prefactor, complex *Ynm) {
    Expansion<p,p-1,m-1>::getYnm(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
    const complex I(0.,1.);
    eim = std::exp(I * real(m * beta));
    Pn = -Pn * (2 * m - 1) * y;
    P0 = Pn;
    Ynm[Index<m,m>::npm] = rhom * P0 * prefactor[Index<m,m>::npm] * eim;
    Ynm[Index<m,m>::nmm] = std::conj(Ynm[Index<m,m>::npm]);
    P1 = P0;
    P0 = x * (2*m+1) * P0;
    rhom *= rho;
    rhon = rhom;
  }
  static inline void getYnmTheta(real &x, real &y, real &rho,
                                 real &rhom, real &rhon, real &beta, complex &eim,
                                 real &Pn, real &P0, real &P1, real &P2,
                                 real *prefactor, complex *Ynm, complex *YnmTheta) {
    Expansion<p,p-1,m-1>::getYnmTheta(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm,YnmTheta);
    const complex I(0.,1.);
    eim = std::exp(I * real(m * beta));
    Pn = -Pn * (2 * m - 1) * y;
    P0 = Pn;
    Ynm[Index<m,m>::npm] = rhom * P0 * prefactor[Index<m,m>::npm] * eim;
    Ynm[Index<m,m>::nmm] = std::conj(Ynm[Index<m,m>::npm]);
    P1 = P0;
    P0 = x * (2*m+1) * P0;
    YnmTheta[Index<m,m>::npm] = rhom * (P0 - (m + 1) * x * P1) / y * prefactor[Index<m,m>::npm] * eim;
    rhom *= rho;
    rhon = rhom;
  }
};

template<int p>
struct Expansion<p,0,0> {
  static inline void getYnm(real &, real &, real &rho,
                            real &rhom, real &rhon, real&, complex&,
                            real&, real &, real &, real&,
                            real*, complex *Ynm) {
    Ynm[0] = rhom;
    rhom *= rho;
    rhon = rhom;
  }
  static inline void getYnmTheta(real &, real &, real &rho,
                                 real &rhom, real &rhon, real&, complex&,
                                 real&, real &, real &, real&,
                                 real*, complex *Ynm, complex *YnmTheta) {
    Ynm[0] = rhom;
    YnmTheta[0] = 0;
    rhom *= rho;
    rhon = rhom;
  }
};

template<int n, int m>
struct Terms {
  static inline void P2M(Mset &M, const real &C, complex *Ynm) {
    Terms<n,m-1>::P2M(M,C,Ynm);
    M[Index<n,m>::nms] += C * Ynm[Index<n,m>::npm];
  }
};

template<int n>
struct Terms<n,0> {
  static inline void P2M(Mset &M, const real &C, complex *Ynm) {
    Terms<n-1,n-1>::P2M(M,C,Ynm);
    M[Index<n,0>::nms] += C * Ynm[Index<n,0>::npm];
  }
};

template<>
struct Terms<0,0> {
  static inline void P2M(Mset &M, const real &C, complex *Ynm) {
    M[Index<0,0>::nms] += C * Ynm[Index<0,0>::npm];
  }
};

template<int p, int j=p-1, int k=p-1, int n=p-1, int m=p-1>
struct M2Ltemplate {
  static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm) {
    M2Ltemplate<p,j,k,n,m-1>::loop(M,L,Cnm,Ynm);
    int nm   = n * n + n + m;
    int nms  = n * (n + 1) / 2 + m;
    int jk = j * j + j + k;
    int jks = j * (j + 1) / 2 + k;
    int jknm = jk * P * P + nm;
    int jnkm = (j + n) * (j + n) + j + n + m - k;
    L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
    nm   = n * n + n - m;
    jknm = jk * P * P + nm;
    jnkm = (j + n) * (j + n) + j + n - m - k;
    L[jks] += std::conj(M[nms]) * Cnm[jknm] * Ynm[jnkm];
  }
};

template<int p, int j, int k, int n>
struct M2Ltemplate<p,j,k,n,1> {
  static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm) {
    M2Ltemplate<p,j,k,n-1,n-1>::loop(M,L,Cnm,Ynm);
    int nm   = n * n + n;
    int nms  = n * (n + 1) / 2;
    int jk = j * j + j + k;
    int jks = j * (j + 1) / 2 + k;
    int jknm = jk * P * P + nm;
    int jnkm = (j + n) * (j + n) + j + n - k;
    L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
    nm   = n * n + n + 1;
    nms  = n * (n + 1) / 2 + 1;
    jknm = jk * P * P + nm;
    jnkm = (j + n) * (j + n) + j + n + 1 - k;
    L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
    nm   = n * n + n - 1;
    jknm = jk * P * P + nm;
    jnkm = (j + n) * (j + n) + j + n - 1 - k;
    L[jks] += std::conj(M[nms]) * Cnm[jknm] * Ynm[jnkm];
  }
};

template<int p, int j, int k>
struct M2Ltemplate<p,j,k,2,1> {
  static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm) {
    M2Ltemplate<p,j,k-1,p-1,p-1>::loop(M,L,Cnm,Ynm);
    int jk = j * j + j + k;
    int jks = j * (j + 1) / 2 + k;
    L[jks] += M[0] * Cnm[jk*P*P] * Ynm[j*j+j-k];
    int jknm = jk * P * P + 6;
    int jnkm = (j + 2) * (j + 2) + j + 2 - k;
    L[jks] += M[3] * Cnm[jknm] * Ynm[jnkm];
    jknm = jk * P * P + 7;
    jnkm = (j + 2) * (j + 2) + j + 3 - k;
    L[jks] += M[4] * Cnm[jknm] * Ynm[jnkm];
    jknm = jk * P * P + 5;
    jnkm = (j + 2) * (j + 2) + j + 1 - k;
    L[jks] += std::conj(M[4]) * Cnm[jknm] * Ynm[jnkm];
  }
};

template<int p, int j>
struct M2Ltemplate<p,j,0,2,1> {
  static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm) {
    M2Ltemplate<p,j-1,j-1,p-1,p-1>::loop(M,L,Cnm,Ynm);
    int jk = j * j + j;
    int jks = j * (j + 1) / 2;
    L[jks] += M[0] * Cnm[jk*P*P] * Ynm[j*j+j];
    int jknm = jk * P * P + 6;
    int jnkm = (j + 2) * (j + 2) + j + 2;
    L[jks] += M[3] * Cnm[jknm] * Ynm[jnkm];
    jknm = jk * P * P + 7;
    jnkm = (j + 2) * (j + 2) + j + 3;
    L[jks] += M[4] * Cnm[jknm] * Ynm[jnkm];
    jknm = jk * P * P + 5;
    jnkm = (j + 2) * (j + 2) + j + 1;
    L[jks] += std::conj(M[4]) * Cnm[jknm] * Ynm[jnkm];
  }
};

template<int p>
struct M2Ltemplate<p,0,0,2,1> {
  static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm) {
    L[0] += M[0] * Cnm[0] * Ynm[0];
    L[0] += M[3] * Cnm[6] * Ynm[6];
    L[0] += M[4] * Cnm[7] * Ynm[7];
    L[0] += std::conj(M[4]) * Cnm[5] * Ynm[5];
  }
};

template<>
void Kernel<Laplace>::initialize() {}

template<>
void Kernel<Laplace>::evalMultipole(real rho, real alpha, real beta, complex *Ynm) const {
  real x = std::cos(alpha);
  real y = std::sin(alpha);
  real rhom=1,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
  complex eim = 1;
  Expansion<P>::getYnm(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
}

template<>
void Kernel<Laplace>::evalMultipoleTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
  real x = std::cos(alpha);
  real y = std::sin(alpha);
  real rhom=1,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
  complex eim = 1;
  Expansion<P>::getYnmTheta(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm,YnmTheta);
}

template<>
void Kernel<Laplace>::evalLocal(real rho, real alpha, real beta, complex *Ynm) const {
  real x = std::cos(alpha);
  real y = std::sin(alpha);
  real invR = 1 / rho;
  real rhom=invR,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
  complex eim = 1;
  Expansion<2*P>::getYnm(x,y,invR,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
}

template<>
void Kernel<Laplace>::evalLocalTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
  real x = std::cos(alpha);
  real y = std::sin(alpha);
  real invR = 1 / rho;
  real rhom=invR,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
  complex eim = 1;
  Expansion<2*P>::getYnmTheta(x,y,invR,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm,YnmTheta);
}

template<>
void Kernel<Laplace>::P2M(C_iter Ci) const {
  real Rmax = 0;
  complex Ynm[4*P*P];
  setCenter(Ci);
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B ) {
    vect dist = B->X - Ci->X;
    real R = std::sqrt(norm(dist));
    if( R > Rmax ) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm);
    Terms<P-1,P-1>::P2M(Ci->M,B->SRC,Ynm);
  }
  Ci->RMAX = Rmax;
  Ci->RCRIT = std::min(Ci->R,Rmax);
}

template<>
void Kernel<Laplace>::M2M(C_iter Ci) const {
  real Rmax = Ci->RMAX;
  const complex I(0.,1.);
  complex Ynm[4*P*P];
  setCenter(Ci);
  for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {
    vect dist = Ci->X - Cj->X;
    real R = std::sqrt(norm(dist)) + Cj->RCRIT;
    if( R > Rmax ) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm);
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
              M += Cj->M[jnkms] * std::pow(I,real(m-abs(m))) * Ynm[nm]
                 * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          for( int m=k; m<=n; ++m ) {
            if( j-n >= m-k ) {
              const int jnkm  = (j - n) * (j - n) + j - n + k - m;
              const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              const int nm    = n * n + n + m;
              M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                 * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        Ci->M[jks] += M * EPS;
      }
    }
  }
  Ci->RCRIT = std::min(Ci->R,Rmax);
}

template<>
void Kernel<Laplace>::M2L(C_iter Ci, C_iter Cj) const {
  complex Ynm[4*P*P];
  vect dist = Ci->X - Cj->X - Xperiodic;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalLocal(rho,alpha,beta,Ynm);
  M2Ltemplate<P>::loop(Cj->M,Ci->L,Cnm,Ynm);
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
    evalLocalTheta(r,theta,phi,Ynm,YnmTheta);
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
  complex Ynm[4*P*P];
  C_iter Cj = Ci0 + Ci->PARENT;
  vect dist = Ci->X - Cj->X;
  real rho, alpha, beta;
  cart2sph(rho,alpha,beta,dist);
  evalMultipole(rho,alpha,beta,Ynm);
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
          L += std::conj(Cj->L[nms]) * Ynm[jnkm]
             * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
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
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B ) {
    vect dist = B->X - Ci->X;
    vect spherical = 0;
    vect cartesian = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipoleTheta(r,theta,phi,Ynm,YnmTheta);
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
