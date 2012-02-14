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

template<int nx, int ny, int nz>
struct Index {
  static const int  I = Index<nx,ny+1,nz-1>::I + 1;
  static const real F = Index<nx,ny,nz-1>::F * nz;
};

template<int nx, int ny>
struct Index<nx,ny,0> {
  static const int  I = Index<nx+1,0,ny-1>::I + 1;
  static const real F = Index<nx,ny-1,0>::F * ny;
};

template<int nx>
struct Index<nx,0,0> {
  static const int  I = Index<0,0,nx-1>::I + 1;
  static const real F = Index<nx-1,0,0>::F * nx;
};

template<>
struct Index<0,0,0> {
  static const int  I = 0;
  static const real F = 1;
};


template<int n, int kx, int ky , int kz, int d>
struct DerivativeTerm {
  static const int coef = 1 - 2 * n;
  static inline real kernel(const Lset &C, const vect &dist) {
    return coef * dist[d] * C[Index<kx,ky,kz>::I];
  }
};

template<int n, int kx, int ky , int kz>
struct DerivativeTerm<n,kx,ky,kz,-1> {
  static const int coef = 1 - n;
  static inline real kernel(const Lset &C, const vect&) {
    return coef * C[Index<kx,ky,kz>::I];
  }
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz, int flag=5>
struct DerivativeSum {
  static const int nextflag = 5 - (kz < nz || kz == 1);
  static const int dim = kz == (nz-1) ? -1 : 2;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,nx,ny,kz-1,nextflag>::loop(C,dist)
         + DerivativeTerm<n,nx,ny,kz-1,dim>::kernel(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,4> {
  static const int nextflag = 3 - (ny == 0);
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,3> {
  static const int nextflag = 3 - (ky < ny || ky == 1);
  static const int dim = ky == (ny-1) ? -1 : 1;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,nx,ky-1,nz,nextflag>::loop(C,dist)
         + DerivativeTerm<n,nx,ky-1,nz,dim>::kernel(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,2> {
  static const int nextflag = 1 - (nx == 0);
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,1> {
  static const int nextflag = 1 - (kx < nx || kx == 1);
  static const int dim = kx == (nx-1) ? -1 : 0;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,kx-1,ny,nz,nextflag>::loop(C,dist)
         + DerivativeTerm<n,kx-1,ny,nz,dim>::kernel(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,0> {
  static inline real loop(const Lset&, const vect&) {
    return 0;
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct DerivativeSum<nx,ny,nz,kx,ky,0,5> {
  static inline real loop(const Lset &C, const vect &dist) {
    return DerivativeSum<nx,ny,nz,nx,ny,0,4>::loop(C,dist);
  }
};


template<int nx, int ny, int nz>
struct Terms {
  static inline void power(Lset &C, const vect &dist) {
    Terms<nx,ny+1,nz-1>::power(C,dist);
    C[Index<nx,ny,nz>::I] = C[Index<nx,ny,nz-1>::I] * dist[2] / nz;
  }
  static inline void derivative(Lset &C, const vect &dist, const real &invR2) {
    static const int n = nx + ny + nz;
    Terms<nx,ny+1,nz-1>::derivative(C,dist,invR2);
    C[Index<nx,ny,nz>::I] = DerivativeSum<nx,ny,nz>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Terms<nx,ny+1,nz-1>::scale(C);
    C[Index<nx,ny,nz>::I] *= Index<nx,ny,nz>::F;
  }
};

template<int nx, int ny>
struct Terms<nx,ny,0> {
  static inline void power(Lset &C, const vect &dist) {
    Terms<nx+1,0,ny-1>::power(C,dist);
    C[Index<nx,ny,0>::I] = C[Index<nx,ny-1,0>::I] * dist[1] / ny;
  }
  static inline void derivative(Lset &C, const vect &dist, const real &invR2) {
    static const int n = nx + ny;
    Terms<nx+1,0,ny-1>::derivative(C,dist,invR2);
    C[Index<nx,ny,0>::I] = DerivativeSum<nx,ny,0>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Terms<nx+1,0,ny-1>::scale(C);
    C[Index<nx,ny,0>::I] *= Index<nx,ny,0>::F;
  }
};

template<int nx>
struct Terms<nx,0,0> {
  static inline void power(Lset &C, const vect &dist) {
    Terms<0,0,nx-1>::power(C,dist);
    C[Index<nx,0,0>::I] = C[Index<nx-1,0,0>::I] * dist[0] / nx;
  }
  static inline void derivative(Lset &C, const vect &dist, const real &invR2) {
    static const int n = nx;
    Terms<0,0,nx-1>::derivative(C,dist,invR2);
    C[Index<nx,0,0>::I] = DerivativeSum<nx,0,0>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Terms<0,0,nx-1>::scale(C);
    C[Index<nx,0,0>::I] *= Index<nx,0,0>::F;
  }
};

template<>
struct Terms<0,0,0> {
  static inline void power(Lset&, const vect&) {}
  static inline void derivative(Lset&, const vect&, const real&) {}
  static inline void scale(Lset&) {}
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz>
struct M2MSum {
  static inline real kernel(const Lset &C, const Mset &M) {
    return M2MSum<nx,ny,nz,kx,ky,kz-1>::kernel(C,M)
         + C[Index<nx-kx,ny-ky,nz-kz>::I]*M[Index<kx,ky,kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct M2MSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return M2MSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M)
         + C[Index<nx-kx,ny-ky,nz>::I]*M[Index<kx,ky,0>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct M2MSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return M2MSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M)
         + C[Index<nx-kx,ny,nz>::I]*M[Index<kx,0,0>::I];
  }
};

template<int nx, int ny, int nz>
struct M2MSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};


template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-nx-ny-nz>
struct M2LSum {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(L,M)
         + M[Index<kx,ky,kz>::I] * L[Index<nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct M2LSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,kx+1,0,ky-1>::kernel(L,M)
         + M[Index<kx,ky,0>::I] * L[Index<nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct M2LSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,0,0,kx-1>::kernel(L,M)
         + M[Index<kx,0,0>::I] * L[Index<nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct M2LSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};


template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-nx-ny-nz>
struct LocalSum {
  static inline real kernel(const Lset &C, const Lset &L) {
    return LocalSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(C,L)
         + C[Index<kx,ky,kz>::I] * L[Index<nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct LocalSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &C, const Lset &L) {
    return LocalSum<nx,ny,nz,kx+1,0,ky-1>::kernel(C,L)
         + C[Index<kx,ky,0>::I] * L[Index<nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct LocalSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &C, const Lset &L) {
    return LocalSum<nx,ny,nz,0,0,kx-1>::kernel(C,L)
         + C[Index<kx,0,0>::I] * L[Index<nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Lset&) { return 0; }
};


template<int nx, int ny, int nz>
struct Upward {
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Upward<nx,ny+1,nz-1>::M2M(MI,C,MJ);
    MI[Index<nx,ny,nz>::I] += M2MSum<nx,ny,nz>::kernel(C,MJ);
  }
};

template<int nx, int ny>
struct Upward<nx,ny,0> {
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Upward<nx+1,0,ny-1>::M2M(MI,C,MJ);
    MI[Index<nx,ny,0>::I] += M2MSum<nx,ny,0>::kernel(C,MJ);
  }
};

template<int nx>
struct Upward<nx,0,0> {
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Upward<0,0,nx-1>::M2M(MI,C,MJ);
    MI[Index<nx,0,0>::I] += M2MSum<nx,0,0>::kernel(C,MJ);
  }
};

template<>
struct Upward<0,0,0> {
  static inline void M2M(Mset&, const Lset&, const Mset&) {}
};


template<int nx, int ny, int nz>
struct Downward {
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Downward<nx,ny+1,nz-1>::M2L(L,C,M);
    L[Index<nx,ny,nz>::I] += M2LSum<nx,ny,nz>::kernel(C,M);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Downward<nx,ny+1,nz-1>::M2P(B,C,M);
    B->TRG[Index<nx,ny,nz>::I] += M2LSum<nx,ny,nz>::kernel(C,M);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Downward<nx,ny+1,nz-1>::L2L(LI,C,LJ);
    LI[Index<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Downward<nx,ny+1,nz-1>::L2P(B,C,L);
    B->TRG[Index<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(C,L);
  }
};

template<int nx, int ny>
struct Downward<nx,ny,0> {
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Downward<nx+1,0,ny-1>::M2L(L,C,M);
    L[Index<nx,ny,0>::I] += M2LSum<nx,ny,0>::kernel(C,M);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Downward<nx+1,0,ny-1>::M2P(B,C,M);
    B->TRG[Index<nx,ny,0>::I] += M2LSum<nx,ny,0>::kernel(C,M);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Downward<nx+1,0,ny-1>::L2L(LI,C,LJ);
    LI[Index<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Downward<nx+1,0,ny-1>::L2P(B,C,L);
    B->TRG[Index<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(C,L);
  }
};

template<int nx>
struct Downward<nx,0,0> {
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Downward<0,0,nx-1>::M2L(L,C,M);
    L[Index<nx,0,0>::I] += M2LSum<nx,0,0>::kernel(C,M);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Downward<0,0,nx-1>::M2P(B,C,M);
    B->TRG[Index<nx,0,0>::I] += M2LSum<nx,0,0>::kernel(C,M);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Downward<0,0,nx-1>::L2L(LI,C,LJ);
    LI[Index<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Downward<0,0,nx-1>::L2P(B,C,L);
    B->TRG[Index<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(C,L);
  }
};

template<>
struct Downward<0,0,0> {
  static inline void M2L(Lset&, const Lset&, const Mset&) {}
  static inline void M2P(B_iter, const Lset&, const Mset&) {}
  static inline void L2L(Lset&, const Lset&, const Lset&) {}
  static inline void L2P(B_iter, const Lset&, const Lset&) {}
};

inline void getCoef(Lset &C, const vect &dist, real &invR2, const real &invR) {
  C[0] = invR;
  Terms<0,0,P>::derivative(C,dist,invR2);
  Terms<0,0,P>::scale(C);
}

inline void sumM2L(Lset &L, const Lset &C, const Mset &M) {
  L += C * M[0];
  for( int i=1; i<MTERM; ++i ) L[0] += M[i] * C[i];
  Downward<0,0,P-1>::M2L(L,C,M);
}

inline void sumM2P(B_iter B, const Lset &C, const Mset &M) {
  B->TRG[0] += C[0] * M[0];
  B->TRG[1] += C[1] * M[0];
  B->TRG[2] += C[2] * M[0];
  B->TRG[3] += C[3] * M[0];
  for( int i=1; i<MTERM; ++i ) B->TRG[0] += M[i] * C[i];
  Downward<0,0,1>::M2P(B,C,M);
}

template<>
void Kernel<Laplace>::initialize() {}

template<>
void Kernel<Laplace>::P2M(C_iter Ci) const {
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = Ci->X - B->X;
    Lset M;
    M[0] = B->SRC;
    Terms<0,0,P-1>::power(M,dist);
    for( int i=0; i<MTERM; ++i ) Ci->M[i] += M[i];
  }
}

template<>
void Kernel<Laplace>::M2M(C_iter Ci, C_iter Cj) const {
  vect dist = Ci->X - Cj->X;
  Mset M;
  Lset C;
  C[0] = 1;
  Terms<0,0,P-1>::power(C,dist);
  M = Cj->M;
  for( int i=0; i<MTERM; ++i ) Ci->M[i] += C[i] * M[0];
  Upward<0,0,P-1>::M2M(Ci->M,C,M);
}

template<>
void Kernel<Laplace>::M2L(C_iter Ci, C_iter Cj) const {
  vect dist = Ci->X - Cj->X - Xperiodic;
  real invR2 = 1 / norm(dist);
  real invR  = std::sqrt(invR2);
  Lset C;
  getCoef(C,dist,invR2,invR);
  sumM2L(Ci->L,C,Cj->M);
}

template<>
void Kernel<Laplace>::M2P(C_iter Ci, C_iter Cj) const {
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
    vect dist = B->X - Cj->X - Xperiodic;
    real invR2 = 1 / norm(dist);
    real invR  = std::sqrt(invR2);
    Lset C;
    getCoef(C,dist,invR2,invR);
    sumM2P(B,C,Cj->M);
  }
}

template<>
void Kernel<Laplace>::L2L(C_iter Ci, C_iter Cj) const {
  vect dist = Ci->X - Cj->X;
  Lset C;
  C[0] = 1;
  Terms<0,0,P>::power(C,dist);

  Ci->L += Cj->L;
  for( int i=1; i<LTERM; ++i ) Ci->L[0] += C[i] * Cj->L[i];
  Downward<0,0,P-1>::L2L(Ci->L,C,Cj->L);
}

template<>
void Kernel<Laplace>::L2P(C_iter Ci) const {
  for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
      vect dist = B->X - Ci->X;
      Lset C, L;
      C[0] = 1;
      Terms<0,0,P>::power(C,dist);

      L = Ci->L;
      B->TRG[0] += L[0];
      B->TRG[1] += L[1];
      B->TRG[2] += L[2];
      B->TRG[3] += L[3];
      for( int i=1; i<LTERM; ++i ) B->TRG[0] += C[i]*L[i];
      Downward<0,0,1>::L2P(B,C,L);
  }
}

template<>
void Kernel<Laplace>::finalize() {}

#include "../kernel/cpuEwaldLaplace.cxx"
