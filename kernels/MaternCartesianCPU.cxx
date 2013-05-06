#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "kernel.h"
using boost::math::cyl_bessel_k;
using boost::math::tgamma;

template<typename T, int nx, int ny, int nz>
struct Index {
  static const int                I = Index<T,nx,ny+1,nz-1>::I + 1;
  static const unsigned long long F = Index<T,nx,ny,nz-1>::F * nz;
};

template<typename T, int nx, int ny>
struct Index<T,nx,ny,0> {
  static const int                I = Index<T,nx+1,0,ny-1>::I + 1;
  static const unsigned long long F = Index<T,nx,ny-1,0>::F * ny;
};

template<typename T, int nx>
struct Index<T,nx,0,0> {
  static const int                I = Index<T,0,0,nx-1>::I + 1;
  static const unsigned long long F = Index<T,nx-1,0,0>::F * nx;
};

template<typename T>
struct Index<T,0,0,0> {
  static const int                I = 0;
  static const unsigned long long F = 1;
};

#if COMkernel
template<>
struct Index<vecM,2,0,0> {
  static const int                I = 1;
  static const unsigned long long F = 2;
};
#endif


template<int kx, int ky , int kz, int d>
struct DerivativeTerm {
  static inline real_t kernel(const vecL &C, const vec3 &dX) {
    return dX[d] * C[Index<vecL,kx,ky,kz>::I];
  }
};

template<int kx, int ky , int kz>
struct DerivativeTerm<kx,ky,kz,-1> {
  static inline real_t kernel(const vecL &C, const vec3&) {
    return C[Index<vecL,kx,ky,kz>::I];
  }
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz, int flag=5>
struct DerivativeSum {
  static const int nextflag = 5 - (kz < nz || kz == 1);
  static const int dim = kz == (nz-1) ? -1 : 2;
  static const int n = nx + ny + nz;
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,kz-1,nextflag>::loop(C,dX)
      - DerivativeTerm<nx,ny,kz-1,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,4> {
  static const int nextflag = 3 - (ny == 0);
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,3> {
  static const int nextflag = 3 - (ky < ny || ky == 1);
  static const int dim = ky == (ny-1) ? -1 : 1;
  static const int n = nx + ny + nz;
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ky-1,nz,nextflag>::loop(C,dX)
      - DerivativeTerm<nx,ky-1,nz,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,2> {
  static const int nextflag = 1 - (nx == 0);
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,1> {
  static const int nextflag = 1 - (kx < nx || kx == 1);
  static const int dim = kx == (nx-1) ? -1 : 0;
  static const int n = nx + ny + nz;
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,kx-1,ny,nz,nextflag>::loop(C,dX)
      - DerivativeTerm<kx-1,ny,nz,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,0> {
  static inline real_t loop(const vecL&, const vec3&) {
    return 0;
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct DerivativeSum<nx,ny,nz,kx,ky,0,5> {
  static inline real_t loop(const vecL &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,0,4>::loop(C,dX);
  }
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz>
struct MultipoleSum {
  static inline real_t kernel(const vecL &C, const vecM &M) {
    return MultipoleSum<nx,ny,nz,kx,ky,kz-1>::kernel(C,M)
         + C[Index<vecL,nx-kx,ny-ky,nz-kz>::I]*M[Index<vecM,kx,ky,kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct MultipoleSum<nx,ny,nz,kx,ky,0> {
  static inline real_t kernel(const vecL &C, const vecM &M) {
    return MultipoleSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M)
         + C[Index<vecL,nx-kx,ny-ky,nz>::I]*M[Index<vecM,kx,ky,0>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct MultipoleSum<nx,ny,nz,kx,0,0> {
  static inline real_t kernel(const vecL &C, const vecM &M) {
    return MultipoleSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M)
         + C[Index<vecL,nx-kx,ny,nz>::I]*M[Index<vecM,kx,0,0>::I];
  }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,0,0> {
  static inline real_t kernel(const vecL&, const vecM&) { return 0; }
};

#if COMkernel
template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,0,1> {
  static inline real_t kernel(const vecL&, const vecM&) { return 0; }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,1,0> {
  static inline real_t kernel(const vecL&, const vecM&) { return 0; }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,1,0,0> {
  static inline real_t kernel(const vecL&, const vecM&) { return 0; }
};
#endif

template<int nx, int ny, int nz, typename T, int kx=0, int ky=0, int kz=P-nx-ny-nz>
struct LocalSum {
  static inline real_t kernel(const T &M, const vecL &L) {
    return LocalSum<nx,ny,nz,T,kx,ky+1,kz-1>::kernel(M,L)
         + M[Index<T,kx,ky,kz>::I] * L[Index<vecL,nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, typename T, int kx, int ky>
struct LocalSum<nx,ny,nz,T,kx,ky,0> {
  static inline real_t kernel(const T &M, const vecL &L) {
    return LocalSum<nx,ny,nz,T,kx+1,0,ky-1>::kernel(M,L)
         + M[Index<T,kx,ky,0>::I] * L[Index<vecL,nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, typename T, int kx>
struct LocalSum<nx,ny,nz,T,kx,0,0> {
  static inline real_t kernel(const T &M, const vecL &L) {
    return LocalSum<nx,ny,nz,T,0,0,kx-1>::kernel(M,L)
         + M[Index<T,kx,0,0>::I] * L[Index<vecL,nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz, typename T>
struct LocalSum<nx,ny,nz,T,0,0,0> {
  static inline real_t kernel(const T&, const vecL&) { return 0; }
};

#if COMkernel
template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,vecM,0,0,1> {
  static inline real_t kernel(const vecM&, const vecL&) { return 0; }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,vecM,0,1,0> {
  static inline real_t kernel(const vecM&, const vecL&) { return 0; }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,vecM,1,0,0> {
  static inline real_t kernel(const vecM&, const vecL&) { return 0; }
};
#endif


template<int nx, int ny, int nz>
struct Kernels {
  static inline void power(vecL &C, const vec3 &dX) {
    Kernels<nx,ny+1,nz-1>::power(C,dX);
    C[Index<vecL,nx,ny,nz>::I] = C[Index<vecL,nx,ny,nz-1>::I] * dX[2] / nz;
  }
  static inline void scale(vecL &C) {
    Kernels<nx,ny+1,nz-1>::scale(C);
    C[Index<vecL,nx,ny,nz>::I] *= Index<vecL,nx,ny,nz>::F;
  }
  static inline void M2M(vecM &MI, const vecL &C, const vecM &MJ) {
    Kernels<nx,ny+1,nz-1>::M2M(MI,C,MJ);
    MI[Index<vecM,nx,ny,nz>::I] += MultipoleSum<nx,ny,nz>::kernel(C,MJ);
  }
  static inline void M2L(vecL &L, const vecL &C, const vecM &M) {
    Kernels<nx,ny+1,nz-1>::M2L(L,C,M);
    L[Index<vecL,nx,ny,nz>::I] += LocalSum<nx,ny,nz,vecM>::kernel(M,C);
  }
  static inline void L2L(vecL &LI, const vecL &C, const vecL &LJ) {
    Kernels<nx,ny+1,nz-1>::L2L(LI,C,LJ);
    LI[Index<vecL,nx,ny,nz>::I] += LocalSum<nx,ny,nz,vecL>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const vecL &C, const vecL &L) {
    Kernels<nx,ny+1,nz-1>::L2P(B,C,L);
    B->TRG[Index<vecL,nx,ny,nz>::I] += LocalSum<nx,ny,nz,vecL>::kernel(C,L);
  }
};

template<int nx, int ny>
struct Kernels<nx,ny,0> {
  static inline void power(vecL &C, const vec3 &dX) {
    Kernels<nx+1,0,ny-1>::power(C,dX);
    C[Index<vecL,nx,ny,0>::I] = C[Index<vecL,nx,ny-1,0>::I] * dX[1] / ny;
  }
  static inline void scale(vecL &C) {
    Kernels<nx+1,0,ny-1>::scale(C);
    C[Index<vecL,nx,ny,0>::I] *= Index<vecL,nx,ny,0>::F;
  }
  static inline void M2M(vecM &MI, const vecL &C, const vecM &MJ) {
    Kernels<nx+1,0,ny-1>::M2M(MI,C,MJ);
    MI[Index<vecM,nx,ny,0>::I] += MultipoleSum<nx,ny,0>::kernel(C,MJ);
  }
  static inline void M2L(vecL &L, const vecL &C, const vecM &M) {
    Kernels<nx+1,0,ny-1>::M2L(L,C,M);
    L[Index<vecL,nx,ny,0>::I] += LocalSum<nx,ny,0,vecM>::kernel(M,C);
  }
  static inline void L2L(vecL &LI, const vecL &C, const vecL &LJ) {
    Kernels<nx+1,0,ny-1>::L2L(LI,C,LJ);
    LI[Index<vecL,nx,ny,0>::I] += LocalSum<nx,ny,0,vecL>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const vecL &C, const vecL &L) {
    Kernels<nx+1,0,ny-1>::L2P(B,C,L);
    B->TRG[Index<vecL,nx,ny,0>::I] += LocalSum<nx,ny,0,vecL>::kernel(C,L);
  }
};

template<int nx>
struct Kernels<nx,0,0> {
  static inline void power(vecL &C, const vec3 &dX) {
    Kernels<0,0,nx-1>::power(C,dX);
    C[Index<vecL,nx,0,0>::I] = C[Index<vecL,nx-1,0,0>::I] * dX[0] / nx;
  }
  static inline void scale(vecL &C) {
    Kernels<0,0,nx-1>::scale(C);
    C[Index<vecL,nx,0,0>::I] *= Index<vecL,nx,0,0>::F;
  }
  static inline void M2M(vecM &MI, const vecL &C, const vecM &MJ) {
    Kernels<0,0,nx-1>::M2M(MI,C,MJ);
    MI[Index<vecM,nx,0,0>::I] += MultipoleSum<nx,0,0>::kernel(C,MJ);
  }
  static inline void M2L(vecL &L, const vecL &C, const vecM &M) {
    Kernels<0,0,nx-1>::M2L(L,C,M);
    L[Index<vecL,nx,0,0>::I] += LocalSum<nx,0,0,vecM>::kernel(M,C);
  }
  static inline void L2L(vecL &LI, const vecL &C, const vecL &LJ) {
    Kernels<0,0,nx-1>::L2L(LI,C,LJ);
    LI[Index<vecL,nx,0,0>::I] += LocalSum<nx,0,0,vecL>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const vecL &C, const vecL &L) {
    Kernels<0,0,nx-1>::L2P(B,C,L);
    B->TRG[Index<vecL,nx,0,0>::I] += LocalSum<nx,0,0,vecL>::kernel(C,L);
  }
};

template<>
struct Kernels<0,0,0> {
  static inline void power(vecL&, const vec3&) {}
  static inline void scale(vecL&) {}
  static inline void M2M(vecM&, const vecL&, const vecM&) {}
  static inline void M2L(vecL&, const vecL&, const vecM&) {}
  static inline void L2L(vecL&, const vecL&, const vecL&) {}
  static inline void L2P(B_iter, const vecL&, const vecL&) {}
};


template<int np, int nx, int ny, int nz>
struct Kernels2 {
  static inline void derivative(vecL &C, vecL &G, const vec3 &dX, const real_t &nu, real_t &coef) {
    static const int n = nx + ny + nz;
    Kernels2<np,nx,ny+1,nz-1>::derivative(C,G,dX,nu,coef);
    C[Index<vecL,nx,ny,nz>::I] = DerivativeSum<nx,ny,nz>::loop(G,dX) / n * coef;
  }
};

template<int np, int nx, int ny>
struct Kernels2<np,nx,ny,0> {
  static inline void derivative(vecL &C, vecL &G, const vec3 &dX, const real_t &nu, real_t &coef) {
    static const int n = nx + ny;
    Kernels2<np,nx+1,0,ny-1>::derivative(C,G,dX,nu,coef);
    C[Index<vecL,nx,ny,0>::I] = DerivativeSum<nx,ny,0>::loop(G,dX) / n * coef;
  }
};

template<int np, int nx>
struct Kernels2<np,nx,0,0> {
  static inline void derivative(vecL &C, vecL &G, const vec3 &dX, const real_t &nu, real_t &coef) {
    static const int n = nx;
    Kernels2<np,0,0,nx-1>::derivative(C,G,dX,nu,coef);
    C[Index<vecL,nx,0,0>::I] = DerivativeSum<nx,0,0>::loop(G,dX) / n * coef;
  }
};

template<int np>
struct Kernels2<np,0,0,0> {
  static inline void derivative(vecL &C, vecL &G, const vec3 &dX, const real_t &nu, real_t &coef) {
    Kernels2<np-1,0,0,np-1>::derivative(G,C,dX,nu,coef);
    static const real_t c = std::sqrt(2 * nu);
    real_t R = c * std::sqrt(norm(dX));
    real_t zR = (-0.577216-log(R/2)) * (R<0.413) + 1 * (R>=0.413);
    static const real_t u = nu - P + np;
    static const real_t gu = tgamma(1-u) / tgamma(u);
    static const real_t aum = std::abs(u-1);
    static const real_t gaum = 1 / tgamma(aum);
    static const real_t au = std::abs(u);
    static const real_t gau = 1 / tgamma(au);
    if (aum < 1e-12) {
      G[0] = cyl_bessel_k(0,R) / zR;
    } else {
      G[0] = std::pow(R/2,aum) * 2 * cyl_bessel_k(aum,R) * gaum;
    }
    if (au < 1e-12) {
      C[0] = cyl_bessel_k(0,R) / zR;
    } else {
      C[0] = std::pow(R/2,au) * 2 * cyl_bessel_k(au,R) * gau;
    }
    real_t hu = 0;
    if (u > 1) {
      hu = 0.5 / (u-1);
    } else if (nu == 0) {
      hu = zR;
    } else if (u > 0 && u < 1) {
      hu = std::pow(R/2,2*u-2) / 2 * gu;
    } else if (u == 0) {
      hu = 1 / (R * R * zR);
    } else {
      hu = -2 * u / (R * R);
    }
    coef = c * c * hu;
  }
};

template<>
struct Kernels2<0,0,0,0> {
  static inline void derivative(vecL, vecL, const vec3, const real_t, real_t) {}
};


template<int PP>
inline void getCoef(vecL &C, const vec3 &dX, const real_t &nu) {
  real_t coef;
  vecL G;
  Kernels2<PP,0,0,PP>::derivative(C,G,dX,nu,coef);
  Kernels<0,0,PP>::scale(C);
}


template<int PP>
inline void sumM2L(vecL &L, const vecL &C, const vecM &M) {
  L += C;
#if COMkernel
  for (int i=1; i<MTERM; i++) L[0] += M[i] * C[i+3];
#else
  for (int i=1; i<MTERM; i++) L[0] += M[i] * C[i];
#endif
  Kernels<0,0,PP-1>::M2L(L,C,M);
}


template<int PP, bool odd>
struct Coefs {
  static const int begin = PP*(PP+1)*(PP+2)/6;
  static const int end = (PP+1)*(PP+2)*(PP+3)/6;
  static inline void negate(vecL &C) {
    for (int i=begin; i<end; i++) C[i] = -C[i];
    Coefs<PP-1,1-odd>::negate(C);
  }
};

template<int PP>
struct Coefs<PP,0> {
  static inline void negate(vecL &C) {
    Coefs<PP-1,1>::negate(C);
  }
};

template<>
struct Coefs<0,0> {
  static inline void negate(vecL){}
};

void Kernel::P2M(C_iter C) const {
  for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {
    vec3 dX = C->X - B->X;
    real_t R = std::sqrt(norm(dX));
    if (R > C->RMAX) C->RMAX = R;
    vecL M;
    M[0] = B->SRC;
    Kernels<0,0,P-1>::power(M,dX/RHO);
#if COMkernel
    C->M[0] += M[0];
    for (int i=1; i<MTERM; i++) C->M[i] += M[i+3];
#else
    for (int i=0; i<MTERM; i++) C->M[i] += M[i];
#endif
  }
#if USE_RMAX
  C->RCRIT = std::min(C->R,C->RMAX);
#else
  C->RCRIT = C->R;
#endif
}

void Kernel::M2M(C_iter Ci, C_iter C0) const {
  for (C_iter Cj=C0+Ci->CHILD; Cj!=C0+Ci->CHILD+Ci->NCHILD; Cj++) {
    vec3 dX = Ci->X - Cj->X;
    real_t R = std::sqrt(norm(dX)) + Cj->RCRIT;
    if (R > Ci->RMAX) Ci->RMAX = R;
    vecM M;
    vecL C;
    C[0] = 1;
    Kernels<0,0,P-1>::power(C,dX/RHO);
    M = Cj->M;
#if COMkernel
    Ci->M[0] += C[0] * M[0];
    for (int i=1; i<MTERM; i++) Ci->M[i] += C[i+3] * M[0];
#else
    for (int i=0; i<MTERM; i++) Ci->M[i] += C[i] * M[0];
#endif
    Kernels<0,0,P-1>::M2M(Ci->M,C,M);
  }
#if USE_RMAX
  Ci->RCRIT = std::min(Ci->R,Ci->RMAX);
#else
  Ci->RCRIT = Ci->R;
#endif
}

void Kernel::M2L(C_iter Ci, C_iter Cj, bool mutual) const {
  vec3 dX = Ci->X - Cj->X - Xperiodic;
  vecL C;
  getCoef<P>(C,dX/RHO,NU);
  C *= Ci->M[0] * Cj->M[0];
  sumM2L<P>(Ci->L,C,Cj->M);
  if (mutual) {
    Coefs<P,P&1>::negate(C);
    sumM2L<P>(Cj->L,C,Ci->M);
  }
}

void Kernel::L2L(C_iter Ci, C_iter Ci0) const {
  C_iter Cj = Ci0 + Ci->PARENT;
  vec3 dX = Ci->X - Cj->X;
  vecL C;
  C[0] = 1;
  Kernels<0,0,P>::power(C,dX/RHO);
  Ci->L /= Ci->M[0];
  Ci->L += Cj->L;
  for (int i=1; i<LTERM; i++) Ci->L[0] += C[i] * Cj->L[i];
  Kernels<0,0,P-1>::L2L(Ci->L,C,Cj->L);
}

void Kernel::L2P(C_iter Ci) const {
  for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {
    vec3 dX = B->X - Ci->X;
    vecL C;
    C[0] = 1;
    Kernels<0,0,P>::power(C,dX/RHO);
    vecL L = Ci->L;
    B->TRG /= B->SRC;
    for (int i=0; i<LTERM; i++) B->TRG[0] += C[i] * L[i];
  }
}
