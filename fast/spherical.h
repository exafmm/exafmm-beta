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
#ifndef kernel_h
#define kernel_h
#include "../include/sort.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

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

class Kernel : public Sort {
private:
  real *factorial, *prefactor, *Anm;
  complex *Ynm, *YnmTheta, *Cnm;

protected:
  vect   X0;
  real   R0;
  C_iter Ci0;
  C_iter Cj0;

private:
  void cart2sph(real& r, real& theta, real& phi, vect dist) const {
    r = std::sqrt(norm(dist))+EPS;
    theta = std::acos(dist[2] / r);
    if( std::abs(dist[0]) + std::abs(dist[1]) < EPS ) {
      phi = 0;
    } else if( std::abs(dist[0]) < EPS ) {
      phi = dist[1] / std::abs(dist[1]) * M_PI * 0.5;
    } else if( dist[0] > 0 ) {
      phi = std::atan(dist[1] / dist[0]);
    } else {
      phi = std::atan(dist[1] / dist[0]) + M_PI;
    }
  }

  template<typename T>
  void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) const {
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]
                 + cos(theta) * cos(phi) / r * spherical[1]
                 - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]
                 + cos(theta) * sin(phi) / r * spherical[1]
                 + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]
                 - sin(theta) / r * spherical[1];
  }

  void evalMultipole(real rho, real alpha, real beta) const {
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real rhom=1,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
    complex eim = 1;
    Expansion<P>::getYnm(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
  }

  void evalMultipoleTheta(real rho, real alpha, real beta) const {
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real rhom=1,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
    complex eim = 1;
    Expansion<P>::getYnmTheta(x,y,rho,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm,YnmTheta);
  }

  void evalLocal(real rho, real alpha, real beta) const {
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real invR = 1 / rho;
    real rhom=invR,rhon=rhom,P0=x,P1=1,P2=1,Pn=1;
    complex eim = 1;
    Expansion<2*P>::getYnm(x,y,invR,rhom,rhon,beta,eim,Pn,P0,P1,P2,prefactor,Ynm);
  }

public:
  Kernel() : X0(0), R0(0) {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new real  [P];
    prefactor = new real  [4*P*P];
    Anm       = new real  [4*P*P];
    Ynm       = new complex [4*P*P];
    YnmTheta  = new complex [4*P*P];
    Cnm       = new complex [P*P*P*P];

    factorial[0] = 1;
    for( int n=1; n!=P; ++n ) {
      factorial[n] = factorial[n-1] * n;
    }

    for( int n=0; n!=2*P; ++n ) {
      for( int m=-n; m<=n; ++m ) {
        int nm = n*n+n+m;
        int nabsm = abs(m);
        real fnmm = 1.0;
        for( int i=1; i<=n-m; ++i ) fnmm *= i;
        real fnpm = 1.0;
        for( int i=1; i<=n+m; ++i ) fnpm *= i;
        real fnma = 1.0;
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;
        real fnpa = 1.0;
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;
        prefactor[nm] = std::sqrt(fnma/fnpa);
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);
      }
    }

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {
      for( int k=-j; k<=j; ++k, ++jk ){
        for( int n=0, nm=0; n!=P; ++n ) {
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {
            int jnkm = (j+n)*(j+n)+j+n+m-k;
            Cnm[jknm] = std::pow(I,abs(k-m)-abs(k)-abs(m))*(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]);
          }
        }
      }
    }
  }

  ~Kernel() {
    delete[] factorial;
    delete[] prefactor;
    delete[] Anm;
    delete[] Ynm;
    delete[] YnmTheta;
    delete[] Cnm;
  }

  void P2P(C_iter Ci, C_iter Cj, bool mutual=true) const {
    for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {
        vect dR = Bi->X - Bj->X;
        real D1 = norm(dR) + EPS2;
        real D0 = Bi->SRC * Bj->SRC;
        real XX = 1.0 / D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        Bj->TRG[0] -= D0 * mutual;
        Bj->TRG[1] += dR[0] * mutual;
        Bj->TRG[2] += dR[1] * mutual;
        Bj->TRG[3] += dR[2] * mutual;
      }
      Bi->TRG[0] += P0;
      Bi->TRG[1] += F0[0];
      Bi->TRG[2] += F0[1];
      Bi->TRG[3] += F0[2];
    }
  }

  void P2P(C_iter C) const {
    unsigned NJ = C->NDLEAF;
    for( B_iter Bi=C->LEAF; Bi!=C->LEAF+C->NDLEAF; ++Bi, --NJ ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter Bj=Bi+1; Bj!=Bi+NJ; ++Bj ) {
        vect dR = Bi->X - Bj->X;
        real D1 = norm(dR) + EPS2;
        real D0 = Bi->SRC * Bj->SRC;
        real XX = 1.0 / D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        Bj->TRG[0] -= D0;
        Bj->TRG[1] += dR[0];
        Bj->TRG[2] += dR[1];
        Bj->TRG[3] += dR[2];
      }
      Bi->TRG[0] += P0;
      Bi->TRG[1] += F0[0];
      Bi->TRG[2] += F0[1];
      Bi->TRG[3] += F0[2];
    }
  }

  void P2M(C_iter C, real &Rmax) const {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vect dist = C->X - B->X;
      real R = std::sqrt(norm(dist));
      if( R > Rmax ) Rmax = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta);
      Terms<P-1,P-1>::P2M(C->M,B->SRC,Ynm);
    }
    C->RCRIT = std::min(C->R,Rmax);
  }

  void M2M(C_iter Ci, real &Rmax) const {
    for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {
      vect dist = Ci->X - Cj->X;
      real R = std::sqrt(norm(dist)) + Cj->RCRIT;
      if( R > Rmax ) Rmax = R;
      const complex I(0.,1.);
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta);
      for( int j=0; j!=P; ++j ) {
        for( int k=0; k<=j; ++k ) {
          int jk = j * j + j + k;
          int jks = j * (j + 1) / 2 + k;
          complex M = 0;
          for( int n=0; n<=j; ++n ) {
            for( int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); ++m ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n + m;
              M += Cj->M[jnkms] * std::pow(I,real(m-abs(m))) * Ynm[nm]
                 * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
            for( int m=k; m<=std::min(n,j+k-n); ++m ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n + m;
              M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                 * real(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          Ci->M[jks] += M;
        }
      }
    }
    Ci->RCRIT = std::min(Ci->R,Rmax);
  }

  void M2L(C_iter Ci, C_iter Cj, bool mutual=true) const {
    vect dist = Ci->X - Cj->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalLocal(rho,alpha,beta);
    M2Ltemplate<P>::loop(Cj->M,Ci->L,Cnm,Ynm);
    if( mutual ) {
      dist = Cj->X - Ci->X;
      cart2sph(rho,alpha,beta,dist);
      evalLocal(rho,alpha,beta);
      M2Ltemplate<P>::loop(Ci->M,Cj->L,Cnm,Ynm);
    }
  }

  void M2P(C_iter Ci, C_iter Cj, bool mutual=true) const {
    const complex I(0.,1.);
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
      vect dist = B->X - Cj->X;
      vect spherical = 0;
      vect cartesian = 0;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalLocal(r,theta,phi);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->TRG[0] -= B->SRC * std::real(Cj->M[nms] * Ynm[nm]);
        spherical[0] -= B->SRC * std::real(Cj->M[nms] * Ynm[nm]) / r * (n+1);
        spherical[1] += B->SRC * std::real(Cj->M[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->TRG[0] -= 2 * B->SRC * std::real(Cj->M[nms] * Ynm[nm]);
          spherical[0] -= 2 * B->SRC * std::real(Cj->M[nms] *Ynm[nm]) / r * (n+1);
          spherical[1] += 2 * B->SRC * std::real(Cj->M[nms] *YnmTheta[nm]);
          spherical[2] += 2 * B->SRC * std::real(Cj->M[nms] *Ynm[nm] * I) * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      B->TRG[1] += cartesian[0];
      B->TRG[2] += cartesian[1];
      B->TRG[3] += cartesian[2];
    }
    if( mutual ) {
      for( B_iter B=Cj->LEAF; B!=Cj->LEAF+Cj->NDLEAF; ++B ) {
        vect dist = B->X - Ci->X;
        vect spherical = 0;
        vect cartesian = 0;
        real r, theta, phi;
        cart2sph(r,theta,phi,dist);
        evalLocal(r,theta,phi);
        for( int n=0; n!=P; ++n ) {
          int nm  = n * n + n;
          int nms = n * (n + 1) / 2;
          B->TRG[0] -= B->SRC * std::real(Ci->M[nms] * Ynm[nm]);
          spherical[0] -= B->SRC * std::real(Ci->M[nms] * Ynm[nm]) / r * (n+1);
          spherical[1] += B->SRC * std::real(Ci->M[nms] * YnmTheta[nm]);
          for( int m=1; m<=n; ++m ) {
            nm  = n * n + n + m;
            nms = n * (n + 1) / 2 + m;
            B->TRG[0] -= 2 * B->SRC * std::real(Ci->M[nms] * Ynm[nm]);
            spherical[0] -= 2 * B->SRC * std::real(Ci->M[nms] *Ynm[nm]) / r * (n+1);
            spherical[1] += 2 * B->SRC * std::real(Ci->M[nms] *YnmTheta[nm]);
            spherical[2] += 2 * B->SRC * std::real(Ci->M[nms] *Ynm[nm] * I) * m;
          }
        }
        sph2cart(r,theta,phi,spherical,cartesian);
        B->TRG[1] += cartesian[0];
        B->TRG[2] += cartesian[1];
        B->TRG[3] += cartesian[2];
      }
    }
  }

  void L2L(C_iter Ci) const {
    C_iter Cj = Ci0 + Ci->PARENT;
    vect dist = Ci->X - Cj->X;
    const complex I(0.,1.);
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,beta);
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
          for( int m=std::max(0,j+k-n); m<=std::min(n,-j+k+n); ++m ) {
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nm   = n * n + n + m;
            int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[nms] * std::pow(I,real(m-k-abs(m-k)))
               * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
          }
        }
        Ci->L[jks] += L;
      }
    }
  }

  void L2P(C_iter Ci) const {
    const complex I(0.,1.);
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B ) {
      vect dist = B->X - Ci->X;
      B->TRG /= B->SRC;
      vect spherical = 0;
      vect cartesian = 0;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipoleTheta(r,theta,phi);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->TRG[0] -= std::real(Ci->L[nms] * Ynm[nm]);
        spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
        spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->TRG[0] -= 2 * std::real(Ci->L[nms] * Ynm[nm]);
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
};

#endif
