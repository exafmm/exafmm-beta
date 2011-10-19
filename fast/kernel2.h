#ifndef kernel_h
#define kernel_h
#include "sort.h"

template<int nx, int ny, int nz>
struct TaylorIndex {
  static const int  I = TaylorIndex<nx,ny+1,nz-1>::I + 1;
  static const real F = TaylorIndex<nx,ny,nz-1>::F * nz;
};

template<int nx, int ny>
struct TaylorIndex<nx,ny,0> {
  static const int  I = TaylorIndex<nx+1,0,ny-1>::I + 1;
  static const real F = TaylorIndex<nx,ny-1,0>::F * ny;
};

template<int nx>
struct TaylorIndex<nx,0,0> {
  static const int  I = TaylorIndex<0,0,nx-1>::I + 1;
  static const real F = TaylorIndex<nx-1,0,0>::F * nx;
};

template<>
struct TaylorIndex<0,0,0> {
  static const int  I = 0;
  static const real F = 1;
};


template<int nx, int ny, int nz>
struct MultipoleCoef {
  static inline void loop(Mset &M, const vect &dist) {
    MultipoleCoef<nx,ny+1,nz-1>::loop(M,dist);
    M[TaylorIndex<nx,ny,nz>::I] = M[TaylorIndex<nx,ny,nz-1>::I] * dist[2] / nz;
  }
};

template<int nx, int ny>
struct MultipoleCoef<nx,ny,0> {
  static inline void loop(Mset &M, const vect &dist) {
    MultipoleCoef<nx+1,0,ny-1>::loop(M,dist);
    M[TaylorIndex<nx,ny,0>::I] = M[TaylorIndex<nx,ny-1,0>::I] * dist[1] / ny;
  }
};

template<int nx>
struct MultipoleCoef<nx,0,0> {
  static inline void loop(Mset &M, const vect &dist) {
    MultipoleCoef<0,0,nx-1>::loop(M,dist);
    M[TaylorIndex<nx,0,0>::I] = M[TaylorIndex<nx-1,0,0>::I] * dist[0] / nx;
  }
};

template<>
struct MultipoleCoef<0,0,0> {
  static inline void loop(Mset&, const vect&) {}
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz>
struct MultipoleSum {
  static inline real kernel(const Mset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx,ky,kz-1>::kernel(C,M) + C[TaylorIndex<kx,ky,kz>::I]*M[TaylorIndex<nx-kx,ny-ky,nz-kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct MultipoleSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Mset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M) + C[TaylorIndex<kx,ky,0>::I]*M[TaylorIndex<nx-kx,ny-ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct MultipoleSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Mset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M) + C[TaylorIndex<kx,0,0>::I]*M[TaylorIndex<nx-kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Mset&, const Mset&) { return 0; }
};


template<int nx, int ny, int nz>
struct MultipoleShift {
  static inline void loop(C_iter CI, const Mset &C, const Mset &M) {
    MultipoleShift<nx,ny+1,nz-1>::loop(CI,C,M);
    CI->M[TaylorIndex<nx,ny,nz>::I] += MultipoleSum<nx,ny,nz>::kernel(C,M);
  }
};

template<int nx, int ny>
struct MultipoleShift<nx,ny,0> {
  static inline void loop(C_iter CI, const Mset &C, const Mset &M) {
    MultipoleShift<nx+1,0,ny-1>::loop(CI,C,M);
    CI->M[TaylorIndex<nx,ny,0>::I] += MultipoleSum<nx,ny,0>::kernel(C,M);
  }
};

template<int nx>
struct MultipoleShift<nx,0,0> {
  static inline void loop(C_iter CI, const Mset &C, const Mset &M) {
    MultipoleShift<0,0,nx-1>::loop(CI,C,M);
    CI->M[TaylorIndex<nx,0,0>::I] += MultipoleSum<nx,0,0>::kernel(C,M);
  }
};

template<>
struct MultipoleShift<0,0,0> {
  static inline void loop(C_iter, const Mset&, const Mset&) {}
};


template<int n, int kx, int ky , int kz, int d>
struct LocalCoefTerm {
  static const int coef = 1 - 2 * n;
  static inline real kernel(const Mset &C, const vect &dist) {
    return coef * dist[d] * C[TaylorIndex<kx,ky,kz>::I];
  }
};

template<int n, int kx, int ky , int kz>
struct LocalCoefTerm<n,kx,ky,kz,-1> {
  static const int coef = 1 - n;
  static inline real kernel(const Mset &C, const vect&) {
    return coef * C[TaylorIndex<kx,ky,kz>::I];
  }
};
  

template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz, int flag=5>
struct LocalCoefSum { 
  static const int nextflag = 5 - (kz < nz || kz == 1);
  static const int dim = kz == (nz-1) ? -1 : 2;
  static const int n = nx + ny + nz;
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,nx,ny,kz-1,nextflag>::loop(C,dist) + LocalCoefTerm<n,nx,ny,kz-1,dim>::kernel(C,dist);
  }
};
  
template<int nx, int ny, int nz, int kx, int ky, int kz>
struct LocalCoefSum<nx,ny,nz,kx,ky,kz,4> {
  static const int nextflag = 3 - (ny == 0);
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct LocalCoefSum<nx,ny,nz,kx,ky,kz,3> {
  static const int nextflag = 3 - (ky < ny || ky == 1);
  static const int dim = ky == (ny-1) ? -1 : 1;
  static const int n = nx + ny + nz;
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,nx,ky-1,nz,nextflag>::loop(C,dist) + LocalCoefTerm<n,nx,ky-1,nz,dim>::kernel(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct LocalCoefSum<nx,ny,nz,kx,ky,kz,2> {
  static const int nextflag = 1 - (nx == 0);
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct LocalCoefSum<nx,ny,nz,kx,ky,kz,1> {
  static const int nextflag = 1 - (kx < nx || kx == 1);
  static const int dim = kx == (nx-1) ? -1 : 0;
  static const int n = nx + ny + nz;
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,kx-1,ny,nz,nextflag>::loop(C,dist) + LocalCoefTerm<n,kx-1,ny,nz,dim>::kernel(C,dist);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct LocalCoefSum<nx,ny,nz,kx,ky,kz,0> {
  static inline real loop(const Mset&, const vect&) { 
    return 0;
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct LocalCoefSum<nx,ny,nz,kx,ky,0,5> {
  static inline real loop(const Mset &C, const vect &dist) {
    return LocalCoefSum<nx,ny,nz,nx,ny,0,4>::loop(C,dist);
  }
};

template<int nx, int ny, int nz>
struct LocalCoef {
  static inline void loop(Mset &C, const vect &dist, const real &invR2) {
    static const int n = nx + ny + nz;
    LocalCoef<nx,ny+1,nz-1>::loop(C,dist,invR2);
    C[TaylorIndex<nx,ny,nz>::I] = LocalCoefSum<nx,ny,nz>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Mset &C) {
    LocalCoef<nx,ny+1,nz-1>::scale(C);
    C[TaylorIndex<nx,ny,nz>::I] *= TaylorIndex<nx,ny,nz>::F;
  }
};

template<int nx, int ny>
struct LocalCoef<nx,ny,0> {
  static inline void loop(Mset &C, const vect &dist, const real &invR2) {
    static const int n = nx + ny;
    LocalCoef<nx+1,0,ny-1>::loop(C,dist,invR2);
    C[TaylorIndex<nx,ny,0>::I] = LocalCoefSum<nx,ny,0>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Mset &C) {
    LocalCoef<nx+1,0,ny-1>::scale(C);
    C[TaylorIndex<nx,ny,0>::I] *= TaylorIndex<nx,ny,0>::F;
  }
};

template<int nx>
struct LocalCoef<nx,0,0> {
  static inline void loop(Mset &C, const vect &dist, const real &invR2) {
    static const int n = nx;
    LocalCoef<0,0,nx-1>::loop(C,dist,invR2);
    C[TaylorIndex<nx,0,0>::I] = LocalCoefSum<nx,0,0>::loop(C,dist) / n * invR2;
  }
  static inline void scale(Mset &C) {
    LocalCoef<0,0,nx-1>::scale(C);
    C[TaylorIndex<nx,0,0>::I] *= TaylorIndex<nx,0,0>::F;
  }
};

template<>
struct LocalCoef<0,0,0> {
  static inline void loop(Mset&, const vect&, const real&) {}
  static inline void scale(Mset&) {}
};


template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-1-nx-ny-nz>
struct LocalSum {
  static inline real kernel(const Mset &L, const Mset &M) {
    return LocalSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(L,M)
         + M[TaylorIndex<kx,ky,kz>::I] * L[TaylorIndex<nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct LocalSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Mset &L, const Mset &M) {
    return LocalSum<nx,ny,nz,kx+1,0,ky-1>::kernel(L,M)
         + M[TaylorIndex<kx,ky,0>::I] * L[TaylorIndex<nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct LocalSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Mset &L, const Mset &M) {
    return LocalSum<nx,ny,nz,0,0,kx-1>::kernel(L,M)
         + M[TaylorIndex<kx,0,0>::I] * L[TaylorIndex<nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Mset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct LocalShift {
  static inline void loopB(B_iter B, const Mset &L, const Mset &M) {
    LocalShift<nx,ny+1,nz-1>::loopB(B,L,M);
    B->TRG[TaylorIndex<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(L,M);
  }
  static inline void loopC(C_iter C, const Mset &L, const Mset &M) {
    LocalShift<nx,ny+1,nz-1>::loopC(C,L,M);
    C->L[TaylorIndex<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(L,M);
  }
};

template<int nx, int ny>
struct LocalShift<nx,ny,0> {
  static inline void loopB(B_iter B, const Mset &L, const Mset &M) {
    LocalShift<nx+1,0,ny-1>::loopB(B,L,M);
    B->TRG[TaylorIndex<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(L,M);
  }
  static inline void loopC(C_iter C, const Mset &L, const Mset &M) {
    LocalShift<nx+1,0,ny-1>::loopC(C,L,M);
    C->L[TaylorIndex<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(L,M);
  }
};

template<int nx>
struct LocalShift<nx,0,0> {
  static inline void loopB(B_iter B, const Mset &L, const Mset &M) {
    LocalShift<0,0,nx-1>::loopB(B,L,M);
    B->TRG[TaylorIndex<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(L,M);
  }
  static inline void loopC(C_iter C, const Mset &L, const Mset &M) {
    LocalShift<0,0,nx-1>::loopC(C,L,M);
    C->L[TaylorIndex<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(L,M);
  }
};

template<>
struct LocalShift<0,0,0> {
  static inline void loopB(B_iter, const Mset&, const Mset&) {}
  static inline void loopC(C_iter, const Mset&, const Mset&) {}
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
      Mset M;
      M[0] = B->SRC[0];
      MultipoleCoef<0,0,P-1>::loop(M,dist);
      C->M += M;
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2M(C_iter CI) {
    for( C_iter CJ=C0+CI->CHILD; CJ!=C0+CI->CHILD+CI->NCHILD; ++CJ ) {
      vect dist = CI->X - CJ->X;
      real R = std::sqrt(norm(dist)) + CJ->RCRIT;
      if( R > DMAX ) DMAX = R;
      Mset C, M;
      C[0] = 1;

      MultipoleCoef<0,0,P-1>::loop(C,dist);

      M = CJ->M;
      CI->M += M;

      MultipoleShift<0,0,P-1>::loop(CI,C,M);

    }
    CI->RCRIT = std::min(CI->R,DMAX);
  }

  void M2L(C_iter CI, C_iter CJ, bool mutual=true) const {
    vect dist = CI->X - CJ->X;
    real invR2 = 1 / norm(dist);
    real invR  = std::sqrt(invR2);
    Mset C, M;

    C[0] = invR;

#if 0
    LocalCoef<0,0,P-1>::loop(C,dist,invR2);
    LocalCoef<0,0,P-1>::scale(C);
#else
    invR2 = -invR2;
    real x = dist[0], y = dist[1], z = dist[2];

    real invR3 = invR * invR2;
    C[1] = x * invR3;
    C[2] = y * invR3;
    C[3] = z * invR3;

    real invR5 = 3 * invR3 * invR2;
    C[4] = x * x * invR5 + invR3;
    C[5] = x * y * invR5;
    C[6] = x * z * invR5;
    C[7] = y * y * invR5 + invR3;
    C[8] = y * z * invR5;
    C[9] = z * z * invR5 + invR3;

    real invR7 = 5 * invR5 * invR2;
    C[10] = x * x * x * invR7 + 3 * x * invR5;
    C[11] = x * x * y * invR7 +     y * invR5;
    C[12] = x * x * z * invR7 +     z * invR5;
    C[13] = x * y * y * invR7 +     x * invR5;
    C[14] = x * y * z * invR7;
    C[15] = x * z * z * invR7 +     x * invR5;
    C[16] = y * y * y * invR7 + 3 * y * invR5;
    C[17] = y * y * z * invR7 +     z * invR5;
    C[18] = y * z * z * invR7 +     y * invR5;
    C[19] = z * z * z * invR7 + 3 * z * invR5;

/*
    real invR9 = 7 * invR7 * invR2;
    C[20] = x * x * x * x * invR9 + 6 * x * x * invR7 + 3 * invR5;
    C[21] = x * x * x * y * invR9 + 3 * x * y * invR7;
    C[22] = x * x * x * z * invR9 + 3 * x * z * invR7;
    C[23] = x * x * y * y * invR9 + (x * x + y * y) * invR7 + invR5;
    C[24] = x * x * y * z * invR9 +     y * z * invR7;
    C[25] = x * x * z * z * invR9 + (x * x + z * z) * invR7 + invR5;
    C[26] = x * y * y * y * invR9 + 3 * x * y * invR7;
    C[27] = x * y * y * z * invR9 +     x * z * invR7;
    C[28] = x * y * z * z * invR9 +     x * y * invR7;
    C[29] = x * z * z * z * invR9 + 3 * x * z * invR7;
    C[30] = y * y * y * y * invR9 + 6 * y * y * invR7 + 3 * invR5;
    C[31] = y * y * y * z * invR9 + 3 * y * z * invR7;
    C[32] = y * y * z * z * invR9 + (y * y + z * z) * invR7 + invR5;
    C[33] = y * z * z * z * invR9 + 3 * y * z * invR7;
    C[34] = z * z * z * z * invR9 + 6 * z * z * invR7 + 3 * invR5;

    real invR11 = 9 * invR9 * invR2;
    C[35] = x * x * x * x * x * invR11 + 10 * x * x * x * invR9 + 15 * x * invR7;
    C[36] = x * x * x * x * y * invR11 +  6 * x * x * y * invR9 +  3 * y * invR7;
    C[37] = x * x * x * x * z * invR11 +  6 * x * x * z * invR9 +  3 * z * invR7;
    C[38] = x * x * x * y * y * invR11 + (x * x * x + 3 * x * y * y) * invR9 + 3 * x * invR7;
    C[39] = x * x * x * y * z * invR11 +  3 * x * y * z * invR9;
    C[40] = x * x * x * z * z * invR11 + (x * x * x + 3 * x * z * z) * invR9 + 3 * x * invR7;
    C[41] = x * x * y * y * y * invR11 + (3 * x * x * y + y * y * y) * invR9 + 3 * y * invR7;
    C[42] = x * x * y * y * z * invR11 + (x * x * z + y * y * z) * invR9 + z * invR7;
    C[43] = x * x * y * z * z * invR11 + (x * x * y + y * z * z) * invR9 + y * invR7;
    C[44] = x * x * z * z * z * invR11 + (3 * x * x * z + z * z * z) * invR9 + 3 * z * invR7;
    C[45] = x * y * y * y * y * invR11 +  6 * x * y * y * invR9 +  3 * x * invR7;
    C[46] = x * y * y * y * z * invR11 +  3 * x * y * z * invR9;
    C[47] = x * y * y * z * z * invR11 + (x * y * y + x * z * z) * invR9 + x * invR7;
    C[48] = x * y * z * z * z * invR11 +  3 * x * y * z * invR9;
    C[49] = x * z * z * z * z * invR11 +  6 * x * z * z * invR9 +  3 * x * invR7;
    C[50] = y * y * y * y * y * invR11 + 10 * y * y * y * invR9 + 15 * y * invR7;
    C[51] = y * y * y * y * z * invR11 +  6 * y * y * z * invR9 +  3 * z * invR7;
    C[52] = y * y * y * z * z * invR11 + (y * y * y + 3 * y * z * z) * invR9 + 3 * y * invR7;
    C[53] = y * y * z * z * z * invR11 + (3 * y * y * z + z * z * z) * invR9 + 3 * z * invR7;
    C[54] = y * z * z * z * z * invR11 +  6 * y * z * z * invR9 +  3 * y * invR7;
    C[55] = z * z * z * z * z * invR11 + 10 * z * z * z * invR9 + 15 * z * invR7;
*/
#endif

    M = CJ->M;
    CI->L += C*M[0];
    for( int i=1; i<10; ++i ) CI->L[0] += M[i]*C[i];
    LocalShift<0,0,P-2>::loopC(CI,C,M);

    invR = mutual;
  }

  void L2L(C_iter CI) const {
    C_iter CJ = C0 + CI->PARENT;
    vect dist = CI->X - CJ->X;
    Lset C, L;
    C[0] = 1;

    MultipoleCoef<0,0,P-1>::loop(C,dist);

    L = CJ->L;
    CI->L += L;
    for( int i=1; i<NCOEF; ++i ) CI->L[0] += C[i]*L[i];

    LocalShift<0,0,P-2>::loopC(CI,L,C);

  }

  void L2P(C_iter CI) const {
    for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NCLEAF; ++B ) {
      vect dist = B->X - CI->X;
      Lset C, L;
      C[0] = 1;
      MultipoleCoef<0,0,P-1>::loop(C,dist);

      L = CI->L;
      B->TRG /= B->SRC[0];

      B->TRG[0] -= L[0];
      B->TRG[1] += L[1];
      B->TRG[2] += L[2];
      B->TRG[3] += L[3];
      for( int i=1; i<NCOEF; ++i ) B->TRG[0] -= C[i]*L[i];

      LocalShift<0,0,1>::loopB(B,L,C);
    }
  }
};

#endif
