#ifndef kernel_h
#define kernel_h
#include "sort.h"

template<int nx, int ny, int nz>
struct Index {
  static const int  M = Index<nx,ny+1,nz-1>::M + 1;
  static const int  I = Index<nx,ny+1,nz-1>::I + 1;
  static const real F = Index<nx,ny,nz-1>::F * nz;
};

template<int nx, int ny>
struct Index<nx,ny,0> {
  static const int  M = Index<nx+1,0,ny-1>::M + 1;
  static const int  I = Index<nx+1,0,ny-1>::I + 1;
  static const real F = Index<nx,ny-1,0>::F * ny;
};

template<int nx>
struct Index<nx,0,0> {
  static const int  M = Index<0,0,nx-1>::M + 1;
  static const int  I = Index<0,0,nx-1>::I + 1;
  static const real F = Index<nx-1,0,0>::F * nx;
};

template<>
struct Index<2,0,0> {
  static const int  M = 1;
  static const int  I = 4;
  static const real F = 2;
};

template<>
struct Index<0,0,1> {
  static const int  M = -1;
  static const int  I = 3;
  static const real F = 1;
};

template<>
struct Index<0,1,0> {
  static const int  M = -1;
  static const int  I = 2;
  static const real F = 1;
};

template<>
struct Index<1,0,0> {
  static const int  M = -1;
  static const int  I = 1;
  static const real F = 1;
};

template<>
struct Index<0,0,0> {
  static const int  M = 0;
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
         + C[Index<nx-kx,ny-ky,nz-kz>::I]*M[Index<kx,ky,kz>::M];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct M2MSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return M2MSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M)
         + C[Index<nx-kx,ny-ky,nz>::I]*M[Index<kx,ky,0>::M];
  }
};

template<int nx, int ny, int nz, int kx>
struct M2MSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return M2MSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M)
         + C[Index<nx-kx,ny,nz>::I]*M[Index<kx,0,0>::M];
  }
};

template<int nx, int ny, int nz>
struct M2MSum<nx,ny,nz,0,0,1> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2MSum<nx,ny,nz,0,1,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2MSum<nx,ny,nz,1,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2MSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};


template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-1-nx-ny-nz>
struct M2LSum {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(L,M)
         + M[Index<kx,ky,kz>::M] * L[Index<nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct M2LSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,kx+1,0,ky-1>::kernel(L,M)
         + M[Index<kx,ky,0>::M] * L[Index<nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct M2LSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &L, const Mset &M) {
    return M2LSum<nx,ny,nz,0,0,kx-1>::kernel(L,M)
         + M[Index<kx,0,0>::M] * L[Index<nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct M2LSum<nx,ny,nz,0,0,1> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2LSum<nx,ny,nz,0,1,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2LSum<nx,ny,nz,1,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct M2LSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};


template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-1-nx-ny-nz>
struct LocalSum {
  static inline real kernel(const Lset &L, const Lset &M) {
    return LocalSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(L,M)
         + M[Index<kx,ky,kz>::I] * L[Index<nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct LocalSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &L, const Lset &M) {
    return LocalSum<nx,ny,nz,kx+1,0,ky-1>::kernel(L,M)
         + M[Index<kx,ky,0>::I] * L[Index<nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct LocalSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &L, const Lset &M) {
    return LocalSum<nx,ny,nz,0,0,kx-1>::kernel(L,M)
         + M[Index<kx,0,0>::I] * L[Index<nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Lset&) { return 0; }
};


template<int nx, int ny, int nz>
struct Upward {
  static inline void M2M(C_iter C, const Lset &L, const Mset &M) {
    Upward<nx,ny+1,nz-1>::M2M(C,L,M);
    C->M[Index<nx,ny,nz>::M] += M2MSum<nx,ny,nz>::kernel(L,M);
  }
};

template<int nx, int ny>
struct Upward<nx,ny,0> {
  static inline void M2M(C_iter C, const Lset &L, const Mset &M) {
    Upward<nx+1,0,ny-1>::M2M(C,L,M);
    C->M[Index<nx,ny,0>::M] += M2MSum<nx,ny,0>::kernel(L,M);
  }
};

template<int nx>
struct Upward<nx,0,0> {
  static inline void M2M(C_iter C, const Lset &L, const Mset &M) {
    Upward<0,0,nx-1>::M2M(C,L,M);
    C->M[Index<nx,0,0>::M] += M2MSum<nx,0,0>::kernel(L,M);
  }
};

template<>
struct Upward<0,0,1> {
  static inline void M2M(C_iter, const Lset&, const Mset&) {}
};

template<>
struct Upward<0,1,0> {
  static inline void M2M(C_iter, const Lset&, const Mset&) {}
};

template<>
struct Upward<1,0,0> {
  static inline void M2M(C_iter, const Lset&, const Mset&) {}
};

template<>
struct Upward<0,0,0> {
  static inline void M2M(C_iter, const Lset&, const Mset&) {}
};


template<int nx, int ny, int nz>
struct Downward {
  static inline void M2L(C_iter C, const Lset &L, const Mset &M) {
    Downward<nx,ny+1,nz-1>::M2L(C,L,M);
    C->L[Index<nx,ny,nz>::I] += M2LSum<nx,ny,nz>::kernel(L,M);
  }
  static inline void L2L(C_iter C, const Lset &L, const Lset &M) {
    Downward<nx,ny+1,nz-1>::L2L(C,L,M);
    C->L[Index<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(L,M);
  }
  static inline void L2P(B_iter B, const Lset &L, const Lset &M) {
    Downward<nx,ny+1,nz-1>::L2P(B,L,M);
    B->TRG[Index<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(L,M);
  }
};

template<int nx, int ny>
struct Downward<nx,ny,0> {
  static inline void M2L(C_iter C, const Lset &L, const Mset &M) {
    Downward<nx+1,0,ny-1>::M2L(C,L,M);
    C->L[Index<nx,ny,0>::I] += M2LSum<nx,ny,0>::kernel(L,M);
  }
  static inline void L2L(C_iter C, const Lset &L, const Lset &M) {
    Downward<nx+1,0,ny-1>::L2L(C,L,M);
    C->L[Index<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(L,M);
  }
  static inline void L2P(B_iter B, const Lset &L, const Lset &M) {
    Downward<nx+1,0,ny-1>::L2P(B,L,M);
    B->TRG[Index<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(L,M);
  }
};

template<int nx>
struct Downward<nx,0,0> {
  static inline void M2L(C_iter C, const Lset &L, const Mset &M) {
    Downward<0,0,nx-1>::M2L(C,L,M);
    C->L[Index<nx,0,0>::I] += M2LSum<nx,0,0>::kernel(L,M);
  }
  static inline void L2L(C_iter C, const Lset &L, const Lset &M) {
    Downward<0,0,nx-1>::L2L(C,L,M);
    C->L[Index<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(L,M);
  }
  static inline void L2P(B_iter B, const Lset &L, const Lset &M) {
    Downward<0,0,nx-1>::L2P(B,L,M);
    B->TRG[Index<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(L,M);
  }
};

template<>
struct Downward<0,0,0> {
  static inline void M2L(C_iter, const Lset&, const Mset&) {}
  static inline void L2L(C_iter, const Lset&, const Lset&) {}
  static inline void L2P(B_iter, const Lset&, const Lset&) {}
};


class Kernel : public Sort {
private:
  real DMAX;

protected:
  vect   X0;
  real   R0;
  C_iter C0;

private:
  real getBmax(vect const&X, C_iter C) {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

public:
  Kernel() : X0(0), R0(0) {}

  ~Kernel() {}

  void setCenter(C_iter C) {
    DMAX = 0;
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += B->SRC[0];
      X += B->X * B->SRC[0];
    }
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      m += c->M[0];
      X += c->X * c->M[0];
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  void P2P(C_iter CI, C_iter CJ, bool mutual=true) const {
    for( B_iter BI=CI->LEAF; BI!=CI->LEAF+CI->NDLEAF; ++BI ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=CJ->LEAF; BJ!=CJ->LEAF+CJ->NDLEAF; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0 * mutual;
        BJ->TRG[1] += dR[0] * mutual;
        BJ->TRG[2] += dR[1] * mutual;
        BJ->TRG[3] += dR[2] * mutual;
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void P2P(C_iter C) const {
    unsigned NJ = C->NDLEAF;
    for( B_iter BI=C->LEAF; BI!=C->LEAF+C->NDLEAF; ++BI, --NJ ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=BI+1; BJ!=BI+NJ; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0;
        BJ->TRG[1] += dR[0];
        BJ->TRG[2] += dR[1];
        BJ->TRG[3] += dR[2];
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void P2M(C_iter C) {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vect dist = C->X - B->X;
      real R = std::sqrt(norm(dist));
      if( R > DMAX ) DMAX = R;
      Lset M;
      M[0] = B->SRC[0];
      Terms<0,0,P-1>::power(M,dist);
      C->M[0] += M[0];
      for( int i=4; i<NCOEF; ++i ) C->M[i-3] += M[i];
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2M(C_iter CI) {
    for( C_iter CJ=C0+CI->CHILD; CJ!=C0+CI->CHILD+CI->NCHILD; ++CJ ) {
      vect dist = CI->X - CJ->X;
      real R = std::sqrt(norm(dist)) + CJ->RCRIT;
      if( R > DMAX ) DMAX = R;
      Mset M;
      Lset C;
      C[0] = 1;
      Terms<0,0,P-1>::power(C,dist);

      M = CJ->M;
      CI->M[0] += C[0] * M[0];
      for( int i=4; i<NCOEF; ++i ) CI->M[i-3] += C[i] * M[0];
      Upward<0,0,P-1>::M2M(CI,C,M);
    }
    CI->RCRIT = std::min(CI->R,DMAX);
  }

  void M2L(C_iter CI, C_iter CJ, bool mutual=true) const {
    vect dist = CI->X - CJ->X;
    real invR2 = 1 / norm(dist);
    real invR  = std::sqrt(invR2);
    Lset C;
    C[0] = invR;

#if 0
    Terms<0,0,P-1>::derivative(C,dist,invR2);
    Terms<0,0,P-1>::scale(C);
#else
    invR2 = -invR2;
    real x = dist[0], y = dist[1], z = dist[2];

    real invR3 = invR * invR2;
    C[1] = x * invR3;
    C[2] = y * invR3;
    C[3] = z * invR3;

    real invR5 = 3 * invR3 * invR2;
    real t = x * invR5;
    C[4] = x * t + invR3;
    C[5] = y * t;
    C[6] = z * t;
    t = y * invR5;
    C[7] = y * t + invR3;
    C[8] = z * t;
    C[9] = z * z * invR5 + invR3;

    real invR7 = 5 * invR5 * invR2;
    t = x * x * invR7;
    C[10] = x * (t + 3 * invR5);
    C[11] = y * (t +     invR5);
    C[12] = z * (t +     invR5);
    t = y * y * invR7;
    C[13] = x * (t +     invR5);
    C[16] = y * (t + 3 * invR5);
    C[17] = z * (t +     invR5);
    t = z * z * invR7;
    C[15] = x * (t +     invR5);
    C[18] = y * (t +     invR5);
    C[19] = z * (t + 3 * invR5);
    C[14] = x * y * z * invR7;

/*
    real invR9 = 7 * invR7 * invR2;
    t = x * x * invR9;
    C[20] = x * x * (t + 6 * invR7) + 3 * invR5;
    C[21] = x * y * (t + 3 * invR7);
    C[22] = x * z * (t + 3 * invR7);
    C[23] = y * y * (t +     invR7) + x * x * invR7 + invR5;
    C[24] = y * z * (t +     invR7);
    C[25] = z * z * (t +     invR7) + x * x * invR7 + invR5;
    t = y * y * invR9;
    C[26] = x * y * (t + 3 * invR7);
    C[27] = x * z * (t +     invR7);
    C[30] = y * y * (t + 6 * invR7) + 3 * invR5;
    C[31] = y * z * (t + 3 * invR7);
    C[32] = z * z * (t +     invR7) + y * y * invR7 + invR5;
    t = z * z * invR9;
    C[28] = x * y * (t +     invR7);
    C[29] = x * z * (t + 3 * invR7);
    C[33] = y * z * (t + 3 * invR7);
    C[34] = z * z * (t + 6 * invR7) + 3 * invR5;

    real invR11 = 9 * invR9 * invR2;
    t = x * x * invR11;
    C[35] = x * x * x * (t + 10 * invR9) + 15 * x * invR7;
    C[36] = x * x * y * (t +  6 * invR9) +  3 * y * invR7;
    C[37] = x * x * z * (t +  6 * invR9) +  3 * z * invR7;
    C[38] = x * y * y * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
    C[39] = x * y * z * (t +  3 * invR9);
    C[40] = x * z * z * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
    C[41] = y * y * y * (t +      invR9) + 3 * x * x * y * invR9 + 3 * y * invR7;
    C[42] = y * y * z * (t +      invR9) + x * x * z * invR9 + z * invR7;
    C[43] = y * z * z * (t +      invR9) + x * x * y * invR9 + y * invR7;
    C[44] = z * z * z * (t +      invR9) + 3 * x * x * z * invR9 + 3 * z * invR7;
    t = y * y * invR11;
    C[45] = x * y * y * (t +  6 * invR9) +  3 * x * invR7;
    C[46] = x * y * z * (t +  3 * invR9);
    C[47] = x * z * z * (t +      invR9) + x * y * y * invR9 + x * invR7;
    C[50] = y * y * y * (t + 10 * invR9) + 15 * y * invR7;
    C[51] = y * y * z * (t +  6 * invR9) + 3 * z * invR7;
    C[52] = y * z * z * (t +  3 * invR9) + y * y * y * invR9 + 3 * y * invR7;
    C[53] = z * z * z * (t +      invR9) + 3 * y * y * z * invR9 + 3 * z * invR7;
    t = z * z * invR11;
    C[48] = x * y * z * (t +  3 * invR9);
    C[49] = x * z * z * (t +  6 * invR9) +  3 * x * invR7;
    C[54] = y * z * z * (t +  6 * invR9) +  3 * y * invR7;
    C[55] = z * z * z * (t + 10 * invR9) + 15 * z * invR7;
*/
#endif

    Mset M = CJ->M;
    CI->L += C * M[0];
    for( int i=4; i<NCOEF; ++i ) CI->L[0] += M[i-3]*C[i];
    Downward<0,0,P-2>::M2L(CI,C,M);

    invR = mutual;
  }

  void L2L(C_iter CI) const {
    C_iter CJ = C0 + CI->PARENT;
    vect dist = CI->X - CJ->X;
    Lset C, L;
    C[0] = 1;
    Terms<0,0,P-1>::power(C,dist);

    L = CJ->L;
    CI->L += L;
    for( int i=1; i<NCOEF; ++i ) CI->L[0] += C[i]*L[i];
    Downward<0,0,P-2>::L2L(CI,L,C);

  }

  void L2P(C_iter CI) const {
    for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NCLEAF; ++B ) {
      vect dist = B->X - CI->X;
      Lset C, L;
      C[0] = 1;
      Terms<0,0,P-1>::power(C,dist);

      L = CI->L;
      B->TRG /= B->SRC[0];
      B->TRG[0] -= L[0];
      B->TRG[1] += L[1];
      B->TRG[2] += L[2];
      B->TRG[3] += L[3];
      for( int i=1; i<NCOEF; ++i ) B->TRG[0] -= C[i]*L[i];
      Downward<0,0,1>::L2P(B,L,C);
    }
  }
};

#endif
