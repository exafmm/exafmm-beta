#ifndef kernel_h
#define kernel_h
#include "sort.h"

struct float4 {
  float x;
  float y;
  float z;
  float w;
};

float4 make_float4(B_iter B) {
  float4 a;
  a.x = B->X[0];
  a.y = B->X[1];
  a.z = B->X[2];
  a.w = B->SRC;
  return a;
}

struct float16 {
  float x[4];
  float y[4];
  float z[4];
  float w[4];
};


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
  static const unsigned long long F = 1.;
};

#if COMkernel
template<>
struct Index<Mset,2,0,0> {
  static const int                I = 1;
  static const unsigned long long F = 2.;
};
#endif


template<int n, int kx, int ky , int kz, int d>
struct DerivativeTerm {
  static const int coef = 1 - 2 * n;
  static inline real kernel(const Lset &C, const vec3 &dX) {
    return coef * dX[d] * C[Index<Lset,kx,ky,kz>::I];
  }
};

template<int n, int kx, int ky , int kz>
struct DerivativeTerm<n,kx,ky,kz,-1> {
  static const int coef = 1 - n;
  static inline real kernel(const Lset &C, const vec3&) {
    return coef * C[Index<Lset,kx,ky,kz>::I];
  }
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz, int flag=5>
struct DerivativeSum {
  static const int nextflag = 5 - (kz < nz || kz == 1);
  static const int dim = kz == (nz-1) ? -1 : 2;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,kz-1,nextflag>::loop(C,dX)
         + DerivativeTerm<n,nx,ny,kz-1,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,4> {
  static const int nextflag = 3 - (ny == 0);
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,3> {
  static const int nextflag = 3 - (ky < ny || ky == 1);
  static const int dim = ky == (ny-1) ? -1 : 1;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ky-1,nz,nextflag>::loop(C,dX)
         + DerivativeTerm<n,nx,ky-1,nz,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,2> {
  static const int nextflag = 1 - (nx == 0);
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,1> {
  static const int nextflag = 1 - (kx < nx || kx == 1);
  static const int dim = kx == (nx-1) ? -1 : 0;
  static const int n = nx + ny + nz;
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,kx-1,ny,nz,nextflag>::loop(C,dX)
         + DerivativeTerm<n,kx-1,ny,nz,dim>::kernel(C,dX);
  }
};

template<int nx, int ny, int nz, int kx, int ky, int kz>
struct DerivativeSum<nx,ny,nz,kx,ky,kz,0> {
  static inline real loop(const Lset&, const vec3&) {
    return 0;
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct DerivativeSum<nx,ny,nz,kx,ky,0,5> {
  static inline real loop(const Lset &C, const vec3 &dX) {
    return DerivativeSum<nx,ny,nz,nx,ny,0,4>::loop(C,dX);
  }
};


template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz>
struct MultipoleSum {
  static inline real kernel(const Lset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx,ky,kz-1>::kernel(C,M)
         + C[Index<Lset,nx-kx,ny-ky,nz-kz>::I]*M[Index<Mset,kx,ky,kz>::I];
  }
};

template<int nx, int ny, int nz, int kx, int ky>
struct MultipoleSum<nx,ny,nz,kx,ky,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M)
         + C[Index<Lset,nx-kx,ny-ky,nz>::I]*M[Index<Mset,kx,ky,0>::I];
  }
};

template<int nx, int ny, int nz, int kx>
struct MultipoleSum<nx,ny,nz,kx,0,0> {
  static inline real kernel(const Lset &C, const Mset &M) {
    return MultipoleSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M)
         + C[Index<Lset,nx-kx,ny,nz>::I]*M[Index<Mset,kx,0,0>::I];
  }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

#if COMkernel
template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,0,1> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,0,1,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};

template<int nx, int ny, int nz>
struct MultipoleSum<nx,ny,nz,1,0,0> {
  static inline real kernel(const Lset&, const Mset&) { return 0; }
};
#endif

template<int nx, int ny, int nz, typename T, int kx=0, int ky=0, int kz=P-nx-ny-nz>
struct LocalSum {
  static inline real kernel(const T &M, const Lset &L) {
    return LocalSum<nx,ny,nz,T,kx,ky+1,kz-1>::kernel(M,L)
         + M[Index<T,kx,ky,kz>::I] * L[Index<Lset,nx+kx,ny+ky,nz+kz>::I];
  }
};

template<int nx, int ny, int nz, typename T, int kx, int ky>
struct LocalSum<nx,ny,nz,T,kx,ky,0> {
  static inline real kernel(const T &M, const Lset &L) {
    return LocalSum<nx,ny,nz,T,kx+1,0,ky-1>::kernel(M,L)
         + M[Index<T,kx,ky,0>::I] * L[Index<Lset,nx+kx,ny+ky,nz>::I];
  }
};

template<int nx, int ny, int nz, typename T, int kx>
struct LocalSum<nx,ny,nz,T,kx,0,0> {
  static inline real kernel(const T &M, const Lset &L) {
    return LocalSum<nx,ny,nz,T,0,0,kx-1>::kernel(M,L)
         + M[Index<T,kx,0,0>::I] * L[Index<Lset,nx+kx,ny,nz>::I];
  }
};

template<int nx, int ny, int nz, typename T>
struct LocalSum<nx,ny,nz,T,0,0,0> {
  static inline real kernel(const T&, const Lset&) { return 0; }
};

#if COMkernel
template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,Mset,0,0,1> {
  static inline real kernel(const Mset&, const Lset&) { return 0; }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,Mset,0,1,0> {
  static inline real kernel(const Mset&, const Lset&) { return 0; }
};

template<int nx, int ny, int nz>
struct LocalSum<nx,ny,nz,Mset,1,0,0> {
  static inline real kernel(const Mset&, const Lset&) { return 0; }
};
#endif


template<int nx, int ny, int nz>
struct Kernels {
  static inline void power(Lset &C, const vec3 &dX) {
    Kernels<nx,ny+1,nz-1>::power(C,dX);
    C[Index<Lset,nx,ny,nz>::I] = C[Index<Lset,nx,ny,nz-1>::I] * dX[2] / nz;
  }
  static inline void derivative(Lset &C, const vec3 &dX, const real &invR2) {
    static const int n = nx + ny + nz;
    Kernels<nx,ny+1,nz-1>::derivative(C,dX,invR2);
    C[Index<Lset,nx,ny,nz>::I] = DerivativeSum<nx,ny,nz>::loop(C,dX) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Kernels<nx,ny+1,nz-1>::scale(C);
    C[Index<Lset,nx,ny,nz>::I] *= Index<Lset,nx,ny,nz>::F;
  }
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Kernels<nx,ny+1,nz-1>::M2M(MI,C,MJ);
    MI[Index<Mset,nx,ny,nz>::I] += MultipoleSum<nx,ny,nz>::kernel(C,MJ);
  }
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Kernels<nx,ny+1,nz-1>::M2L(L,C,M);
    L[Index<Lset,nx,ny,nz>::I] += LocalSum<nx,ny,nz,Mset>::kernel(M,C);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Kernels<nx,ny+1,nz-1>::M2P(B,C,M);
    B->TRG[Index<Lset,nx,ny,nz>::I] += LocalSum<nx,ny,nz,Mset>::kernel(M,C);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Kernels<nx,ny+1,nz-1>::L2L(LI,C,LJ);
    LI[Index<Lset,nx,ny,nz>::I] += LocalSum<nx,ny,nz,Lset>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Kernels<nx,ny+1,nz-1>::L2P(B,C,L);
    B->TRG[Index<Lset,nx,ny,nz>::I] += LocalSum<nx,ny,nz,Lset>::kernel(C,L);
  }
};

template<int nx, int ny>
struct Kernels<nx,ny,0> {
  static inline void power(Lset &C, const vec3 &dX) {
    Kernels<nx+1,0,ny-1>::power(C,dX);
    C[Index<Lset,nx,ny,0>::I] = C[Index<Lset,nx,ny-1,0>::I] * dX[1] / ny;
  }
  static inline void derivative(Lset &C, const vec3 &dX, const real &invR2) {
    static const int n = nx + ny;
    Kernels<nx+1,0,ny-1>::derivative(C,dX,invR2);
    C[Index<Lset,nx,ny,0>::I] = DerivativeSum<nx,ny,0>::loop(C,dX) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Kernels<nx+1,0,ny-1>::scale(C);
    C[Index<Lset,nx,ny,0>::I] *= Index<Lset,nx,ny,0>::F;
  }
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Kernels<nx+1,0,ny-1>::M2M(MI,C,MJ);
    MI[Index<Mset,nx,ny,0>::I] += MultipoleSum<nx,ny,0>::kernel(C,MJ);
  }
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Kernels<nx+1,0,ny-1>::M2L(L,C,M);
    L[Index<Lset,nx,ny,0>::I] += LocalSum<nx,ny,0,Mset>::kernel(M,C);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Kernels<nx+1,0,ny-1>::M2P(B,C,M);
    B->TRG[Index<Lset,nx,ny,0>::I] += LocalSum<nx,ny,0,Mset>::kernel(M,C);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Kernels<nx+1,0,ny-1>::L2L(LI,C,LJ);
    LI[Index<Lset,nx,ny,0>::I] += LocalSum<nx,ny,0,Lset>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Kernels<nx+1,0,ny-1>::L2P(B,C,L);
    B->TRG[Index<Lset,nx,ny,0>::I] += LocalSum<nx,ny,0,Lset>::kernel(C,L);
  }
};

template<int nx>
struct Kernels<nx,0,0> {
  static inline void power(Lset &C, const vec3 &dX) {
    Kernels<0,0,nx-1>::power(C,dX);
    C[Index<Lset,nx,0,0>::I] = C[Index<Lset,nx-1,0,0>::I] * dX[0] / nx;
  }
  static inline void derivative(Lset &C, const vec3 &dX, const real &invR2) {
    static const int n = nx;
    Kernels<0,0,nx-1>::derivative(C,dX,invR2);
    C[Index<Lset,nx,0,0>::I] = DerivativeSum<nx,0,0>::loop(C,dX) / n * invR2;
  }
  static inline void scale(Lset &C) {
    Kernels<0,0,nx-1>::scale(C);
    C[Index<Lset,nx,0,0>::I] *= Index<Lset,nx,0,0>::F;
  }
  static inline void M2M(Mset &MI, const Lset &C, const Mset &MJ) {
    Kernels<0,0,nx-1>::M2M(MI,C,MJ);
    MI[Index<Mset,nx,0,0>::I] += MultipoleSum<nx,0,0>::kernel(C,MJ);
  }
  static inline void M2L(Lset &L, const Lset &C, const Mset &M) {
    Kernels<0,0,nx-1>::M2L(L,C,M);
    L[Index<Lset,nx,0,0>::I] += LocalSum<nx,0,0,Mset>::kernel(M,C);
  }
  static inline void M2P(B_iter B, const Lset &C, const Mset &M) {
    Kernels<0,0,nx-1>::M2P(B,C,M);
    B->TRG[Index<Lset,nx,0,0>::I] += LocalSum<nx,0,0,Mset>::kernel(M,C);
  }
  static inline void L2L(Lset &LI, const Lset &C, const Lset &LJ) {
    Kernels<0,0,nx-1>::L2L(LI,C,LJ);
    LI[Index<Lset,nx,0,0>::I] += LocalSum<nx,0,0,Lset>::kernel(C,LJ);
  }
  static inline void L2P(B_iter B, const Lset &C, const Lset &L) {
    Kernels<0,0,nx-1>::L2P(B,C,L);
    B->TRG[Index<Lset,nx,0,0>::I] += LocalSum<nx,0,0,Lset>::kernel(C,L);
  }
};

template<>
struct Kernels<0,0,0> {
  static inline void power(Lset&, const vec3&) {}
  static inline void derivative(Lset&, const vec3&, const real&) {}
  static inline void scale(Lset&) {}
  static inline void M2M(Mset&, const Lset&, const Mset&) {}
  static inline void M2L(Lset&, const Lset&, const Mset&) {}
  static inline void M2P(B_iter, const Lset&, const Mset&) {}
  static inline void L2L(Lset&, const Lset&, const Lset&) {}
  static inline void L2P(B_iter, const Lset&, const Lset&) {}
};


template<int PP>
inline void getCoef(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  C[0] = invR;
  Kernels<0,0,PP>::derivative(C,dX,invR2);
  Kernels<0,0,PP>::scale(C);
}

template<>
inline void getCoef<1>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  C[0] = invR;
  invR2 = -invR2;
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  C[1] = x * invR3;
  C[2] = y * invR3;
  C[3] = z * invR3;
}

template<>
inline void getCoef<2>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  getCoef<1>(C,dX,invR2,invR);
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real t = x * invR5;
  C[4] = x * t + invR3;
  C[5] = y * t;
  C[6] = z * t;
  t = y * invR5;
  C[7] = y * t + invR3;
  C[8] = z * t;
  C[9] = z * z * invR5 + invR3;
}
template<>
inline void getCoef<3>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
#if 1
  getCoef<2>(C,dX,invR2,invR);
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real invR7 = 5 * invR5 * invR2;
  real t = x * x * invR7;
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
#else
  real* C_arr = (real*) C;

  __m128 result[5];
  __m128 term0;
  __m128 term1;

  invR2 = -invR2;
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real t = x * invR5;
  real invR7 = 5 * invR5 * invR2;
  real t1 = y * invR5;
  real t2 = x * x * invR7;
  real t3 = y * y * invR7;
  real t4 = z * z * invR7;

  term0 = _mm_set_ps(z, y, x, invR);
  term1 = _mm_set_ps(invR3, invR3, invR3, 1);
  result[0] = _mm_mul_ps(term0, term1);

  term0 = _mm_set_ps(y, z, y, x);
  term1 = _mm_set_ps(t1, t, t, t);
  result[1] = _mm_mul_ps(term0, term1);
  term1 = _mm_set_ps(invR3, 0, 0, invR3);
  result[1] = _mm_add_ps(result[1], term1);

  term0 = _mm_set_ps(y, x, z, z);
  term1 = _mm_set_ps(t2 + invR5, t2 + 3 * invR5, z * invR5, t1);
  result[2] = _mm_mul_ps(term0, term1);
  term1 = _mm_set_ps(0, 0, invR3, 0);
  result[2] = _mm_add_ps(result[2], term1);

  term0 = _mm_set_ps(x, x, x, z);
  term1 = _mm_set_ps(t4 + invR5, y * z * invR7, t3 + invR5, t2 + invR5);
  result[3] = _mm_mul_ps(term0, term1);

  term0 = _mm_set_ps(z, y, z, y);
  term1 = _mm_set_ps(t4 + 3 * invR5, t4 + invR5, t3 + invR5, t3 + 3 * invR5);
  result[4] = _mm_mul_ps(term0, term1);

  _mm_store_ps(C_arr, result[0]);
  _mm_store_ps(C_arr + 4, result[1]);
  _mm_store_ps(C_arr + 8, result[2]);
  _mm_store_ps(C_arr + 12, result[3]);
  _mm_store_ps(C_arr + 16, result[4]);
#endif
}

template<>
inline void getCoef<4>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  getCoef<3>(C,dX,invR2,invR);
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real invR7 = 5 * invR5 * invR2;
  real invR9 = 7 * invR7 * invR2;
  real t = x * x * invR9;
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
}

template<>
inline void getCoef<5>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  getCoef<4>(C,dX,invR2,invR);
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real invR7 = 5 * invR5 * invR2;
  real invR9 = 7 * invR7 * invR2;
  real invR11 = 9 * invR9 * invR2;
  real t = x * x * invR11;
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
}

template<>
inline void getCoef<6>(Lset &C, const vec3 &dX, real &invR2, const real &invR) {
  getCoef<5>(C,dX,invR2,invR);
  real x = dX[0], y = dX[1], z = dX[2];
  real invR3 = invR * invR2;
  real invR5 = 3 * invR3 * invR2;
  real invR7 = 5 * invR5 * invR2;
  real invR9 = 7 * invR7 * invR2;
  real invR11 = 9 * invR9 * invR2;
  real invR13 = 11 * invR11 * invR2;
  real t = x * x * invR13;
  C[56] = x * x * x * x * (t + 15 * invR11) + 45 * x * x * invR9 + 15 * invR7;
  C[57] = x * x * x * y * (t + 10 * invR11) + 15 * x * y * invR9;
  C[58] = x * x * x * z * (t + 10 * invR11) + 15 * x * z * invR9;
  C[59] = x * x * y * y * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * y * y) * invR9 + 3 * invR7;
  C[60] = x * x * y * z * (t +  6 * invR11) + 3 * y * z * invR9;
  C[61] = x * x * z * z * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * z * z) * invR9 + 3 * invR7;
  C[62] = x * y * y * y * (t +  3 * invR11) + 3 * x * x * x * y * invR11 + 9 * x * y * invR9;
  C[63] = x * y * y * z * (t +  3 * invR11) + x * x * x * z * invR11 + 3 * x * z * invR9;
  C[64] = x * y * z * z * (t +  3 * invR11) + x * x * x * y * invR11 + 3 * x * y * invR9;
  C[65] = x * z * z * z * (t +  3 * invR11) + 3 * x * x * x * z * invR11 + 9 * x * z * invR9;
  C[66] = y * y * y * y * (t +      invR11) + 6 * x * x * y * y * invR11 + (3 * x * x + 6 * y * y) * invR9 + 3 * invR7;
  C[67] = y * y * y * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
  C[68] = y * y * z * z * (t +      invR11) + (x * x * y * y + x * x * z * z) * invR11 + (x * x + y * y + z * z) * invR9 + invR7;
  C[69] = y * z * z * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
  C[70] = z * z * z * z * (t +      invR11) + 6 * x * x * z * z * invR11 + (3 * x * x + 6 * z * z) * invR9 + 3 * invR7;
  t = y * y * invR13;
  C[71] = x * y * y * y * (t + 10 * invR11) + 15 * x * y * invR9;
  C[72] = x * y * y * z * (t +  6 * invR11) + 3 * x * z * invR9;
  C[73] = x * y * z * z * (t +  3 * invR11) + x * y * y * y * invR11 + 3 * x * y * invR9;
  C[74] = x * z * z * z * (t +      invR11) + 3 * x * y * y * z * invR11 + 3 * x * z * invR9;
  C[77] = y * y * y * y * (t + 15 * invR11) + 45 * y * y * invR9 + 15 * invR7;
  C[78] = y * y * y * z * (t + 10 * invR11) + 15 * y * z * invR9;
  C[79] = y * y * z * z * (t +  6 * invR11) + y * y * y * y * invR11 + (6 * y * y + 3 * z * z) * invR9 + 3 * invR7;
  C[80] = y * z * z * z * (t +  3 * invR11) + 3 * y * y * y * z * invR11 + 9 * y * z * invR9;
  C[81] = z * z * z * z * (t +      invR11) + 6 * y * y * z * z * invR11 + (3 * y * y + 6 * z * z) * invR9 + 3 * invR7;
  t = z * z * invR13;
  C[75] = x * y * z * z * (t +  6 * invR11) + 3 * x * y * invR9;
  C[76] = x * z * z * z * (t + 10 * invR11) + 15 * x * z * invR9;
  C[82] = y * z * z * z * (t + 10 * invR11) + 15 * y * z * invR9;
  C[83] = z * z * z * z * (t + 15 * invR11) + 45 * z * z * invR9 + 15 * invR7;
}


template<int PP>
inline void sumM2L(Lset &L, const Lset &C, const Mset &M) {
  L = C;
#if COMkernel
  for( int i=1; i<MTERM; ++i ) L[0] += M[i] * C[i+3];
#else
  for( int i=1; i<MTERM; ++i ) L[0] += M[i] * C[i];
#endif
  Kernels<0,0,PP-1>::M2L(L,C,M);
}

template<>
inline void sumM2L<1>(Lset &L, const Lset &C, const Mset&) {
  L = C;
}

template<>
inline void sumM2L<2>(Lset &L, const Lset &C, const Mset &M) {
  sumM2L<1>(L,C,M);
#if COMkernel
#else
  L[0] += M[1]*C[1]+M[2]*C[2]+M[3]*C[3];
  L[1] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6];
  L[2] += M[1]*C[5]+M[2]*C[7]+M[3]*C[8];
  L[3] += M[1]*C[6]+M[2]*C[8]+M[3]*C[9];
#endif
}

template<>
inline void sumM2L<3>(Lset &L, const Lset &C, const Mset &M) {
  sumM2L<2>(L,C,M);
#if COMkernel
  L[0] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6]+M[4]*C[7]+M[5]*C[8]+M[6]*C[9];
  L[1] += M[1]*C[10]+M[2]*C[11]+M[3]*C[12]+M[4]*C[13]+M[5]*C[14]+M[6]*C[15];
  L[2] += M[1]*C[11]+M[2]*C[13]+M[3]*C[14]+M[4]*C[16]+M[5]*C[17]+M[6]*C[18];
  L[3] += M[1]*C[12]+M[2]*C[14]+M[3]*C[15]+M[4]*C[17]+M[5]*C[18]+M[6]*C[19];
#else
  L[0] += M[4]*C[4]+M[5]*C[5]+M[6]*C[6]+M[7]*C[7]+M[8]*C[8]+M[9]*C[9];
  L[1] += M[4]*C[10]+M[5]*C[11]+M[6]*C[12]+M[7]*C[13]+M[8]*C[14]+M[9]*C[15];
  L[2] += M[4]*C[11]+M[5]*C[13]+M[6]*C[14]+M[7]*C[16]+M[8]*C[17]+M[9]*C[18];
  L[3] += M[4]*C[12]+M[5]*C[14]+M[6]*C[15]+M[7]*C[17]+M[8]*C[18]+M[9]*C[19];
  L[4] += M[1]*C[10]+M[2]*C[11]+M[3]*C[12];
  L[5] += M[1]*C[11]+M[2]*C[13]+M[3]*C[14];
  L[6] += M[1]*C[12]+M[2]*C[14]+M[3]*C[15];
  L[7] += M[1]*C[13]+M[2]*C[16]+M[3]*C[17];
  L[8] += M[1]*C[14]+M[2]*C[17]+M[3]*C[18];
  L[9] += M[1]*C[15]+M[2]*C[18]+M[3]*C[19];
#endif
}

template<>
inline void sumM2L<4>(Lset &L, const Lset &C, const Mset &M) {
  sumM2L<3>(L,C,M);
#if COMkernel
  L[0] += M[7]*C[10]+M[8]*C[11]+M[9]*C[12]+M[10]*C[13]+M[11]*C[14]+M[12]*C[15]+M[13]*C[16]+M[14]*C[17]+M[15]*C[18]+M[16]*C[19];
  L[1] += M[7]*C[20]+M[8]*C[21]+M[9]*C[22]+M[10]*C[23]+M[11]*C[24]+M[12]*C[25]+M[13]*C[26]+M[14]*C[27]+M[15]*C[28]+M[16]*C[29];
  L[2] += M[7]*C[21]+M[8]*C[23]+M[9]*C[24]+M[10]*C[26]+M[11]*C[27]+M[12]*C[28]+M[13]*C[30]+M[14]*C[31]+M[15]*C[32]+M[16]*C[33];
  L[3] += M[7]*C[22]+M[8]*C[24]+M[9]*C[25]+M[10]*C[27]+M[11]*C[28]+M[12]*C[29]+M[13]*C[31]+M[14]*C[32]+M[15]*C[33]+M[16]*C[34];
  L[4] += M[1]*C[20]+M[2]*C[21]+M[3]*C[22]+M[4]*C[23]+M[5]*C[24]+M[6]*C[25];
  L[5] += M[1]*C[21]+M[2]*C[23]+M[3]*C[24]+M[4]*C[26]+M[5]*C[27]+M[6]*C[28];
  L[6] += M[1]*C[22]+M[2]*C[24]+M[3]*C[25]+M[4]*C[27]+M[5]*C[28]+M[6]*C[29];
  L[7] += M[1]*C[23]+M[2]*C[26]+M[3]*C[27]+M[4]*C[30]+M[5]*C[31]+M[6]*C[32];
  L[8] += M[1]*C[24]+M[2]*C[27]+M[3]*C[28]+M[4]*C[31]+M[5]*C[32]+M[6]*C[33];
  L[9] += M[1]*C[25]+M[2]*C[28]+M[3]*C[29]+M[4]*C[32]+M[5]*C[33]+M[6]*C[34];
#else
  L[0] += M[10]*C[10]+M[11]*C[11]+M[12]*C[12]+M[13]*C[13]+M[14]*C[14]+M[15]*C[15]+M[16]*C[16]+M[17]*C[17]+M[18]*C[18]+M[19]*C[19];
  L[1] += M[10]*C[20]+M[11]*C[21]+M[12]*C[22]+M[13]*C[23]+M[14]*C[24]+M[15]*C[25]+M[16]*C[26]+M[17]*C[27]+M[18]*C[28]+M[19]*C[29];
  L[2] += M[10]*C[21]+M[11]*C[23]+M[12]*C[24]+M[13]*C[26]+M[14]*C[27]+M[15]*C[28]+M[16]*C[30]+M[17]*C[31]+M[18]*C[32]+M[19]*C[33];
  L[3] += M[10]*C[22]+M[11]*C[24]+M[12]*C[25]+M[13]*C[27]+M[14]*C[28]+M[15]*C[29]+M[16]*C[31]+M[17]*C[32]+M[18]*C[33]+M[19]*C[34];
  L[4] += M[4]*C[20]+M[5]*C[21]+M[6]*C[22]+M[7]*C[23]+M[8]*C[24]+M[9]*C[25];
  L[5] += M[4]*C[21]+M[5]*C[23]+M[6]*C[24]+M[7]*C[26]+M[8]*C[27]+M[9]*C[28];
  L[6] += M[4]*C[22]+M[5]*C[24]+M[6]*C[25]+M[7]*C[27]+M[8]*C[28]+M[9]*C[29];
  L[7] += M[4]*C[23]+M[5]*C[26]+M[6]*C[27]+M[7]*C[30]+M[8]*C[31]+M[9]*C[32];
  L[8] += M[4]*C[24]+M[5]*C[27]+M[6]*C[28]+M[7]*C[31]+M[8]*C[32]+M[9]*C[33];
  L[9] += M[4]*C[25]+M[5]*C[28]+M[6]*C[29]+M[7]*C[32]+M[8]*C[33]+M[9]*C[34];
  L[10] += M[1]*C[20]+M[2]*C[21]+M[3]*C[22];
  L[11] += M[1]*C[21]+M[2]*C[23]+M[3]*C[24];
  L[12] += M[1]*C[22]+M[2]*C[24]+M[3]*C[25];
  L[13] += M[1]*C[23]+M[2]*C[26]+M[3]*C[27];
  L[14] += M[1]*C[24]+M[2]*C[27]+M[3]*C[28];
  L[15] += M[1]*C[25]+M[2]*C[28]+M[3]*C[29];
  L[16] += M[1]*C[26]+M[2]*C[30]+M[3]*C[31];
  L[17] += M[1]*C[27]+M[2]*C[31]+M[3]*C[32];
  L[18] += M[1]*C[28]+M[2]*C[32]+M[3]*C[33];
  L[19] += M[1]*C[29]+M[2]*C[33]+M[3]*C[34];
#endif
}

template<>
inline void sumM2L<5>(Lset &L, const Lset &C, const Mset &M) {
  sumM2L<4>(L,C,M);
#if COMkernel
  L[0] += M[17]*C[20]+M[18]*C[21]+M[19]*C[22]+M[20]*C[23]+M[21]*C[24]+M[22]*C[25]+M[23]*C[26]+M[24]*C[27]+M[25]*C[28]+M[26]*C[29]+M[27]*C[30]+M[28]*C[31]+M[29]*C[32]+M[30]*C[33]+M[31]*C[34];
  L[1] += M[17]*C[35]+M[18]*C[36]+M[19]*C[37]+M[20]*C[38]+M[21]*C[39]+M[22]*C[40]+M[23]*C[41]+M[24]*C[42]+M[25]*C[43]+M[26]*C[44]+M[27]*C[45]+M[28]*C[46]+M[29]*C[47]+M[30]*C[48]+M[31]*C[49];
  L[2] += M[17]*C[36]+M[18]*C[38]+M[19]*C[39]+M[20]*C[41]+M[21]*C[42]+M[22]*C[43]+M[23]*C[45]+M[24]*C[46]+M[25]*C[47]+M[26]*C[48]+M[27]*C[50]+M[28]*C[51]+M[29]*C[52]+M[30]*C[53]+M[31]*C[54];
  L[3] += M[17]*C[37]+M[18]*C[39]+M[19]*C[40]+M[20]*C[42]+M[21]*C[43]+M[22]*C[44]+M[23]*C[46]+M[24]*C[47]+M[25]*C[48]+M[26]*C[49]+M[27]*C[51]+M[28]*C[52]+M[29]*C[53]+M[30]*C[54]+M[31]*C[55];
  L[4] += M[7]*C[35]+M[8]*C[36]+M[9]*C[37]+M[10]*C[38]+M[11]*C[39]+M[12]*C[40]+M[13]*C[41]+M[14]*C[42]+M[15]*C[43]+M[16]*C[44];
  L[5] += M[7]*C[36]+M[8]*C[38]+M[9]*C[39]+M[10]*C[41]+M[11]*C[42]+M[12]*C[43]+M[13]*C[45]+M[14]*C[46]+M[15]*C[47]+M[16]*C[48];
  L[6] += M[7]*C[37]+M[8]*C[39]+M[9]*C[40]+M[10]*C[42]+M[11]*C[43]+M[12]*C[44]+M[13]*C[46]+M[14]*C[47]+M[15]*C[48]+M[16]*C[49];
  L[7] += M[7]*C[38]+M[8]*C[41]+M[9]*C[42]+M[10]*C[45]+M[11]*C[46]+M[12]*C[47]+M[13]*C[50]+M[14]*C[51]+M[15]*C[52]+M[16]*C[53];
  L[8] += M[7]*C[39]+M[8]*C[42]+M[9]*C[43]+M[10]*C[46]+M[11]*C[47]+M[12]*C[48]+M[13]*C[51]+M[14]*C[52]+M[15]*C[53]+M[16]*C[54];
  L[9] += M[7]*C[40]+M[8]*C[43]+M[9]*C[44]+M[10]*C[47]+M[11]*C[48]+M[12]*C[49]+M[13]*C[52]+M[14]*C[53]+M[15]*C[54]+M[16]*C[55];
  L[10] += M[1]*C[35]+M[2]*C[36]+M[3]*C[37]+M[4]*C[38]+M[5]*C[39]+M[6]*C[40];
  L[11] += M[1]*C[36]+M[2]*C[38]+M[3]*C[39]+M[4]*C[41]+M[5]*C[42]+M[6]*C[43];
  L[12] += M[1]*C[37]+M[2]*C[39]+M[3]*C[40]+M[4]*C[42]+M[5]*C[43]+M[6]*C[44];
  L[13] += M[1]*C[38]+M[2]*C[41]+M[3]*C[42]+M[4]*C[45]+M[5]*C[46]+M[6]*C[47];
  L[14] += M[1]*C[39]+M[2]*C[42]+M[3]*C[43]+M[4]*C[46]+M[5]*C[47]+M[6]*C[48];
  L[15] += M[1]*C[40]+M[2]*C[43]+M[3]*C[44]+M[4]*C[47]+M[5]*C[48]+M[6]*C[49];
  L[16] += M[1]*C[41]+M[2]*C[45]+M[3]*C[46]+M[4]*C[50]+M[5]*C[51]+M[6]*C[52];
  L[17] += M[1]*C[42]+M[2]*C[46]+M[3]*C[47]+M[4]*C[51]+M[5]*C[52]+M[6]*C[53];
  L[18] += M[1]*C[43]+M[2]*C[47]+M[3]*C[48]+M[4]*C[52]+M[5]*C[53]+M[6]*C[54];
  L[19] += M[1]*C[44]+M[2]*C[48]+M[3]*C[49]+M[4]*C[53]+M[5]*C[54]+M[6]*C[55];
#else
  L[0] += M[20]*C[20]+M[21]*C[21]+M[22]*C[22]+M[23]*C[23]+M[24]*C[24]+M[25]*C[25]+M[26]*C[26]+M[27]*C[27]+M[28]*C[28]+M[29]*C[29]+M[30]*C[30]+M[31]*C[31]+M[32]*C[32]+M[33]*C[33]+M[34]*C[34];
  L[1] += M[20]*C[35]+M[21]*C[36]+M[22]*C[37]+M[23]*C[38]+M[24]*C[39]+M[25]*C[40]+M[26]*C[41]+M[27]*C[42]+M[28]*C[43]+M[29]*C[44]+M[30]*C[45]+M[31]*C[46]+M[32]*C[47]+M[33]*C[48]+M[34]*C[49];
  L[2] += M[20]*C[36]+M[21]*C[38]+M[22]*C[39]+M[23]*C[41]+M[24]*C[42]+M[25]*C[43]+M[26]*C[45]+M[27]*C[46]+M[28]*C[47]+M[29]*C[48]+M[30]*C[50]+M[31]*C[51]+M[32]*C[52]+M[33]*C[53]+M[34]*C[54];
  L[3] += M[20]*C[37]+M[21]*C[39]+M[22]*C[40]+M[23]*C[42]+M[24]*C[43]+M[25]*C[44]+M[26]*C[46]+M[27]*C[47]+M[28]*C[48]+M[29]*C[49]+M[30]*C[51]+M[31]*C[52]+M[32]*C[53]+M[33]*C[54]+M[34]*C[55];
  L[4] += M[10]*C[35]+M[11]*C[36]+M[12]*C[37]+M[13]*C[38]+M[14]*C[39]+M[15]*C[40]+M[16]*C[41]+M[17]*C[42]+M[18]*C[43]+M[19]*C[44];
  L[5] += M[10]*C[36]+M[11]*C[38]+M[12]*C[39]+M[13]*C[41]+M[14]*C[42]+M[15]*C[43]+M[16]*C[45]+M[17]*C[46]+M[18]*C[47]+M[19]*C[48];
  L[6] += M[10]*C[37]+M[11]*C[39]+M[12]*C[40]+M[13]*C[42]+M[14]*C[43]+M[15]*C[44]+M[16]*C[46]+M[17]*C[47]+M[18]*C[48]+M[19]*C[49];
  L[7] += M[10]*C[38]+M[11]*C[41]+M[12]*C[42]+M[13]*C[45]+M[14]*C[46]+M[15]*C[47]+M[16]*C[50]+M[17]*C[51]+M[18]*C[52]+M[19]*C[53];
  L[8] += M[10]*C[39]+M[11]*C[42]+M[12]*C[43]+M[13]*C[46]+M[14]*C[47]+M[15]*C[48]+M[16]*C[51]+M[17]*C[52]+M[18]*C[53]+M[19]*C[54];
  L[9] += M[10]*C[40]+M[11]*C[43]+M[12]*C[44]+M[13]*C[47]+M[14]*C[48]+M[15]*C[49]+M[16]*C[52]+M[17]*C[53]+M[18]*C[54]+M[19]*C[55];
  L[10] += M[4]*C[35]+M[5]*C[36]+M[6]*C[37]+M[7]*C[38]+M[8]*C[39]+M[9]*C[40];
  L[11] += M[4]*C[36]+M[5]*C[38]+M[6]*C[39]+M[7]*C[41]+M[8]*C[42]+M[9]*C[43];
  L[12] += M[4]*C[37]+M[5]*C[39]+M[6]*C[40]+M[7]*C[42]+M[8]*C[43]+M[9]*C[44];
  L[13] += M[4]*C[38]+M[5]*C[41]+M[6]*C[42]+M[7]*C[45]+M[8]*C[46]+M[9]*C[47];
  L[14] += M[4]*C[39]+M[5]*C[42]+M[6]*C[43]+M[7]*C[46]+M[8]*C[47]+M[9]*C[48];
  L[15] += M[4]*C[40]+M[5]*C[43]+M[6]*C[44]+M[7]*C[47]+M[8]*C[48]+M[9]*C[49];
  L[16] += M[4]*C[41]+M[5]*C[45]+M[6]*C[46]+M[7]*C[50]+M[8]*C[51]+M[9]*C[52];
  L[17] += M[4]*C[42]+M[5]*C[46]+M[6]*C[47]+M[7]*C[51]+M[8]*C[52]+M[9]*C[53];
  L[18] += M[4]*C[43]+M[5]*C[47]+M[6]*C[48]+M[7]*C[52]+M[8]*C[53]+M[9]*C[54];
  L[19] += M[4]*C[44]+M[5]*C[48]+M[6]*C[49]+M[7]*C[53]+M[8]*C[54]+M[9]*C[55];
  L[20] += M[1]*C[35]+M[2]*C[36]+M[3]*C[37];
  L[21] += M[1]*C[36]+M[2]*C[38]+M[3]*C[39];
  L[22] += M[1]*C[37]+M[2]*C[39]+M[3]*C[40];
  L[23] += M[1]*C[38]+M[2]*C[41]+M[3]*C[42];
  L[24] += M[1]*C[39]+M[2]*C[42]+M[3]*C[43];
  L[25] += M[1]*C[40]+M[2]*C[43]+M[3]*C[44];
  L[26] += M[1]*C[41]+M[2]*C[45]+M[3]*C[46];
  L[27] += M[1]*C[42]+M[2]*C[46]+M[3]*C[47];
  L[28] += M[1]*C[43]+M[2]*C[47]+M[3]*C[48];
  L[29] += M[1]*C[44]+M[2]*C[48]+M[3]*C[49];
  L[30] += M[1]*C[45]+M[2]*C[50]+M[3]*C[51];
  L[31] += M[1]*C[46]+M[2]*C[51]+M[3]*C[52];
  L[32] += M[1]*C[47]+M[2]*C[52]+M[3]*C[53];
  L[33] += M[1]*C[48]+M[2]*C[53]+M[3]*C[54];
  L[34] += M[1]*C[49]+M[2]*C[54]+M[3]*C[55];
#endif
}

template<>
inline void sumM2L<6>(Lset &L, const Lset &C, const Mset &M) {
  sumM2L<5>(L,C,M);
#if COMkernel
  L[0] += M[32]*C[35]+M[33]*C[36]+M[34]*C[37]+M[35]*C[38]+M[36]*C[39]+M[37]*C[40]+M[38]*C[41]+M[39]*C[42]+M[40]*C[43]+M[41]*C[44]+M[42]*C[45]+M[43]*C[46]+M[44]*C[47]+M[45]*C[48]+M[46]*C[49]+M[47]*C[50]+M[48]*C[51]+M[49]*C[52]+M[50]*C[53]+M[51]*C[54]+M[52]*C[55];
  L[1] += M[32]*C[56]+M[33]*C[57]+M[34]*C[58]+M[35]*C[59]+M[36]*C[60]+M[37]*C[61]+M[38]*C[62]+M[39]*C[63]+M[40]*C[64]+M[41]*C[65]+M[42]*C[66]+M[43]*C[67]+M[44]*C[68]+M[45]*C[69]+M[46]*C[70]+M[47]*C[71]+M[48]*C[72]+M[49]*C[73]+M[50]*C[74]+M[51]*C[75]+M[52]*C[76];
  L[2] += M[32]*C[57]+M[33]*C[59]+M[34]*C[60]+M[35]*C[62]+M[36]*C[63]+M[37]*C[64]+M[38]*C[66]+M[39]*C[67]+M[40]*C[68]+M[41]*C[69]+M[42]*C[71]+M[43]*C[72]+M[44]*C[73]+M[45]*C[74]+M[46]*C[75]+M[47]*C[77]+M[48]*C[78]+M[49]*C[79]+M[50]*C[80]+M[51]*C[81]+M[52]*C[82];
  L[3] += M[32]*C[58]+M[33]*C[60]+M[34]*C[61]+M[35]*C[63]+M[36]*C[64]+M[37]*C[65]+M[38]*C[67]+M[39]*C[68]+M[40]*C[69]+M[41]*C[70]+M[42]*C[72]+M[43]*C[73]+M[44]*C[74]+M[45]*C[75]+M[46]*C[76]+M[47]*C[78]+M[48]*C[79]+M[49]*C[80]+M[50]*C[81]+M[51]*C[82]+M[52]*C[83];
  L[4] += M[17]*C[56]+M[18]*C[57]+M[19]*C[58]+M[20]*C[59]+M[21]*C[60]+M[22]*C[61]+M[23]*C[62]+M[24]*C[63]+M[25]*C[64]+M[26]*C[65]+M[27]*C[66]+M[28]*C[67]+M[29]*C[68]+M[30]*C[69]+M[31]*C[70];
  L[5] += M[17]*C[57]+M[18]*C[59]+M[19]*C[60]+M[20]*C[62]+M[21]*C[63]+M[22]*C[64]+M[23]*C[66]+M[24]*C[67]+M[25]*C[68]+M[26]*C[69]+M[27]*C[71]+M[28]*C[72]+M[29]*C[73]+M[30]*C[74]+M[31]*C[75];
  L[6] += M[17]*C[58]+M[18]*C[60]+M[19]*C[61]+M[20]*C[63]+M[21]*C[64]+M[22]*C[65]+M[23]*C[67]+M[24]*C[68]+M[25]*C[69]+M[26]*C[70]+M[27]*C[72]+M[28]*C[73]+M[29]*C[74]+M[30]*C[75]+M[31]*C[76];
  L[7] += M[17]*C[59]+M[18]*C[62]+M[19]*C[63]+M[20]*C[66]+M[21]*C[67]+M[22]*C[68]+M[23]*C[71]+M[24]*C[72]+M[25]*C[73]+M[26]*C[74]+M[27]*C[77]+M[28]*C[78]+M[29]*C[79]+M[30]*C[80]+M[31]*C[81];
  L[8] += M[17]*C[60]+M[18]*C[63]+M[19]*C[64]+M[20]*C[67]+M[21]*C[68]+M[22]*C[69]+M[23]*C[72]+M[24]*C[73]+M[25]*C[74]+M[26]*C[75]+M[27]*C[78]+M[28]*C[79]+M[29]*C[80]+M[30]*C[81]+M[31]*C[82];
  L[9] += M[17]*C[61]+M[18]*C[64]+M[19]*C[65]+M[20]*C[68]+M[21]*C[69]+M[22]*C[70]+M[23]*C[73]+M[24]*C[74]+M[25]*C[75]+M[26]*C[76]+M[27]*C[79]+M[28]*C[80]+M[29]*C[81]+M[30]*C[82]+M[31]*C[83];
  L[7] += M[7]*C[56]+M[8]*C[57]+M[9]*C[58]+M[10]*C[59]+M[11]*C[60]+M[12]*C[61]+M[13]*C[62]+M[14]*C[63]+M[15]*C[64]+M[16]*C[65];
  L[8] += M[7]*C[57]+M[8]*C[59]+M[9]*C[60]+M[10]*C[62]+M[11]*C[63]+M[12]*C[64]+M[13]*C[66]+M[14]*C[67]+M[15]*C[68]+M[16]*C[69];
  L[9] += M[7]*C[58]+M[8]*C[60]+M[9]*C[61]+M[10]*C[63]+M[11]*C[64]+M[12]*C[65]+M[13]*C[67]+M[14]*C[68]+M[15]*C[69]+M[16]*C[70];
  L[10] += M[7]*C[59]+M[8]*C[62]+M[9]*C[63]+M[10]*C[66]+M[11]*C[67]+M[12]*C[68]+M[13]*C[71]+M[14]*C[72]+M[15]*C[73]+M[16]*C[74];
  L[11] += M[7]*C[60]+M[8]*C[63]+M[9]*C[64]+M[10]*C[67]+M[11]*C[68]+M[12]*C[69]+M[13]*C[72]+M[14]*C[73]+M[15]*C[74]+M[16]*C[75];
  L[12] += M[7]*C[61]+M[8]*C[64]+M[9]*C[65]+M[10]*C[68]+M[11]*C[69]+M[12]*C[70]+M[13]*C[73]+M[14]*C[74]+M[15]*C[75]+M[16]*C[76];
  L[13] += M[7]*C[62]+M[8]*C[66]+M[9]*C[67]+M[10]*C[71]+M[11]*C[72]+M[12]*C[73]+M[13]*C[77]+M[14]*C[78]+M[15]*C[79]+M[16]*C[80];
  L[14] += M[7]*C[63]+M[8]*C[67]+M[9]*C[68]+M[10]*C[72]+M[11]*C[73]+M[12]*C[74]+M[13]*C[78]+M[14]*C[79]+M[15]*C[80]+M[16]*C[81];
  L[15] += M[7]*C[64]+M[8]*C[68]+M[9]*C[69]+M[10]*C[73]+M[11]*C[74]+M[12]*C[75]+M[13]*C[79]+M[14]*C[80]+M[15]*C[81]+M[16]*C[82];
  L[16] += M[7]*C[65]+M[8]*C[69]+M[9]*C[70]+M[10]*C[74]+M[11]*C[75]+M[12]*C[76]+M[13]*C[80]+M[14]*C[81]+M[15]*C[82]+M[16]*C[83];
  L[20] += M[1]*C[56]+M[2]*C[57]+M[3]*C[58]+M[4]*C[59]+M[5]*C[60]+M[6]*C[61];
  L[21] += M[1]*C[57]+M[2]*C[59]+M[3]*C[60]+M[4]*C[62]+M[5]*C[63]+M[6]*C[64];
  L[22] += M[1]*C[58]+M[2]*C[60]+M[3]*C[61]+M[4]*C[63]+M[5]*C[64]+M[6]*C[65];
  L[23] += M[1]*C[59]+M[2]*C[62]+M[3]*C[63]+M[4]*C[66]+M[5]*C[67]+M[6]*C[68];
  L[24] += M[1]*C[60]+M[2]*C[63]+M[3]*C[64]+M[4]*C[67]+M[5]*C[68]+M[6]*C[69];
  L[25] += M[1]*C[61]+M[2]*C[64]+M[3]*C[65]+M[4]*C[68]+M[5]*C[69]+M[6]*C[70];
  L[26] += M[1]*C[62]+M[2]*C[66]+M[3]*C[67]+M[4]*C[71]+M[5]*C[72]+M[6]*C[73];
  L[27] += M[1]*C[63]+M[2]*C[67]+M[3]*C[68]+M[4]*C[72]+M[5]*C[73]+M[6]*C[74];
  L[28] += M[1]*C[64]+M[2]*C[68]+M[3]*C[69]+M[4]*C[73]+M[5]*C[74]+M[6]*C[75];
  L[29] += M[1]*C[65]+M[2]*C[69]+M[3]*C[70]+M[4]*C[74]+M[5]*C[75]+M[6]*C[76];
  L[30] += M[1]*C[66]+M[2]*C[71]+M[3]*C[72]+M[4]*C[77]+M[5]*C[78]+M[6]*C[79];
  L[31] += M[1]*C[67]+M[2]*C[72]+M[3]*C[73]+M[4]*C[78]+M[5]*C[79]+M[6]*C[80];
  L[32] += M[1]*C[68]+M[2]*C[73]+M[3]*C[74]+M[4]*C[79]+M[5]*C[80]+M[6]*C[81];
  L[33] += M[1]*C[69]+M[2]*C[74]+M[3]*C[75]+M[4]*C[80]+M[5]*C[81]+M[6]*C[82];
  L[34] += M[1]*C[70]+M[2]*C[75]+M[3]*C[76]+M[4]*C[81]+M[5]*C[82]+M[6]*C[83];
#else
  L[0] += M[35]*C[35]+M[36]*C[36]+M[37]*C[37]+M[38]*C[38]+M[39]*C[39]+M[40]*C[40]+M[41]*C[41]+M[42]*C[42]+M[43]*C[43]+M[44]*C[44]+M[45]*C[45]+M[46]*C[46]+M[47]*C[47]+M[48]*C[48]+M[49]*C[49]+M[50]*C[50]+M[51]*C[51]+M[52]*C[52]+M[53]*C[53]+M[54]*C[54]+M[55]*C[55];
  L[1] += M[35]*C[56]+M[36]*C[57]+M[37]*C[58]+M[38]*C[59]+M[39]*C[60]+M[40]*C[61]+M[41]*C[62]+M[42]*C[63]+M[43]*C[64]+M[44]*C[65]+M[45]*C[66]+M[46]*C[67]+M[47]*C[68]+M[48]*C[69]+M[49]*C[70]+M[50]*C[71]+M[51]*C[72]+M[52]*C[73]+M[53]*C[74]+M[54]*C[75]+M[55]*C[76];
  L[2] += M[35]*C[57]+M[36]*C[59]+M[37]*C[60]+M[38]*C[62]+M[39]*C[63]+M[40]*C[64]+M[41]*C[66]+M[42]*C[67]+M[43]*C[68]+M[44]*C[69]+M[45]*C[71]+M[46]*C[72]+M[47]*C[73]+M[48]*C[74]+M[49]*C[75]+M[50]*C[77]+M[51]*C[78]+M[52]*C[79]+M[53]*C[80]+M[54]*C[81]+M[55]*C[82];
  L[3] += M[35]*C[58]+M[36]*C[60]+M[37]*C[61]+M[38]*C[63]+M[39]*C[64]+M[40]*C[65]+M[41]*C[67]+M[42]*C[68]+M[43]*C[69]+M[44]*C[70]+M[45]*C[72]+M[46]*C[73]+M[47]*C[74]+M[48]*C[75]+M[49]*C[76]+M[50]*C[78]+M[51]*C[79]+M[52]*C[80]+M[53]*C[81]+M[54]*C[82]+M[55]*C[83];
  L[4] += M[20]*C[56]+M[21]*C[57]+M[22]*C[58]+M[23]*C[59]+M[24]*C[60]+M[25]*C[61]+M[26]*C[62]+M[27]*C[63]+M[28]*C[64]+M[29]*C[65]+M[30]*C[66]+M[31]*C[67]+M[32]*C[68]+M[33]*C[69]+M[34]*C[70];
  L[5] += M[20]*C[57]+M[21]*C[59]+M[22]*C[60]+M[23]*C[62]+M[24]*C[63]+M[25]*C[64]+M[26]*C[66]+M[27]*C[67]+M[28]*C[68]+M[29]*C[69]+M[30]*C[71]+M[31]*C[72]+M[32]*C[73]+M[33]*C[74]+M[34]*C[75];
  L[6] += M[20]*C[58]+M[21]*C[60]+M[22]*C[61]+M[23]*C[63]+M[24]*C[64]+M[25]*C[65]+M[26]*C[67]+M[27]*C[68]+M[28]*C[69]+M[29]*C[70]+M[30]*C[72]+M[31]*C[73]+M[32]*C[74]+M[33]*C[75]+M[34]*C[76];
  L[7] += M[20]*C[59]+M[21]*C[62]+M[22]*C[63]+M[23]*C[66]+M[24]*C[67]+M[25]*C[68]+M[26]*C[71]+M[27]*C[72]+M[28]*C[73]+M[29]*C[74]+M[30]*C[77]+M[31]*C[78]+M[32]*C[79]+M[33]*C[80]+M[34]*C[81];
  L[8] += M[20]*C[60]+M[21]*C[63]+M[22]*C[64]+M[23]*C[67]+M[24]*C[68]+M[25]*C[69]+M[26]*C[72]+M[27]*C[73]+M[28]*C[74]+M[29]*C[75]+M[30]*C[78]+M[31]*C[79]+M[32]*C[80]+M[33]*C[81]+M[34]*C[82];
  L[9] += M[20]*C[61]+M[21]*C[64]+M[22]*C[65]+M[23]*C[68]+M[24]*C[69]+M[25]*C[70]+M[26]*C[73]+M[27]*C[74]+M[28]*C[75]+M[29]*C[76]+M[30]*C[79]+M[31]*C[80]+M[32]*C[81]+M[33]*C[82]+M[34]*C[83];
  L[10] += M[10]*C[56]+M[11]*C[57]+M[12]*C[58]+M[13]*C[59]+M[14]*C[60]+M[15]*C[61]+M[16]*C[62]+M[17]*C[63]+M[18]*C[64]+M[19]*C[65];
  L[11] += M[10]*C[57]+M[11]*C[59]+M[12]*C[60]+M[13]*C[62]+M[14]*C[63]+M[15]*C[64]+M[16]*C[66]+M[17]*C[67]+M[18]*C[68]+M[19]*C[69];
  L[12] += M[10]*C[58]+M[11]*C[60]+M[12]*C[61]+M[13]*C[63]+M[14]*C[64]+M[15]*C[65]+M[16]*C[67]+M[17]*C[68]+M[18]*C[69]+M[19]*C[70];
  L[13] += M[10]*C[59]+M[11]*C[62]+M[12]*C[63]+M[13]*C[66]+M[14]*C[67]+M[15]*C[68]+M[16]*C[71]+M[17]*C[72]+M[18]*C[73]+M[19]*C[74];
  L[14] += M[10]*C[60]+M[11]*C[63]+M[12]*C[64]+M[13]*C[67]+M[14]*C[68]+M[15]*C[69]+M[16]*C[72]+M[17]*C[73]+M[18]*C[74]+M[19]*C[75];
  L[15] += M[10]*C[61]+M[11]*C[64]+M[12]*C[65]+M[13]*C[68]+M[14]*C[69]+M[15]*C[70]+M[16]*C[73]+M[17]*C[74]+M[18]*C[75]+M[19]*C[76];
  L[16] += M[10]*C[62]+M[11]*C[66]+M[12]*C[67]+M[13]*C[71]+M[14]*C[72]+M[15]*C[73]+M[16]*C[77]+M[17]*C[78]+M[18]*C[79]+M[19]*C[80];
  L[17] += M[10]*C[63]+M[11]*C[67]+M[12]*C[68]+M[13]*C[72]+M[14]*C[73]+M[15]*C[74]+M[16]*C[78]+M[17]*C[79]+M[18]*C[80]+M[19]*C[81];
  L[18] += M[10]*C[64]+M[11]*C[68]+M[12]*C[69]+M[13]*C[73]+M[14]*C[74]+M[15]*C[75]+M[16]*C[79]+M[17]*C[80]+M[18]*C[81]+M[19]*C[82];
  L[19] += M[10]*C[65]+M[11]*C[69]+M[12]*C[70]+M[13]*C[74]+M[14]*C[75]+M[15]*C[76]+M[16]*C[80]+M[17]*C[81]+M[18]*C[82]+M[19]*C[83];
  L[20] += M[4]*C[56]+M[5]*C[57]+M[6]*C[58]+M[7]*C[59]+M[8]*C[60]+M[9]*C[61];
  L[21] += M[4]*C[57]+M[5]*C[59]+M[6]*C[60]+M[7]*C[62]+M[8]*C[63]+M[9]*C[64];
  L[22] += M[4]*C[58]+M[5]*C[60]+M[6]*C[61]+M[7]*C[63]+M[8]*C[64]+M[9]*C[65];
  L[23] += M[4]*C[59]+M[5]*C[62]+M[6]*C[63]+M[7]*C[66]+M[8]*C[67]+M[9]*C[68];
  L[24] += M[4]*C[60]+M[5]*C[63]+M[6]*C[64]+M[7]*C[67]+M[8]*C[68]+M[9]*C[69];
  L[25] += M[4]*C[61]+M[5]*C[64]+M[6]*C[65]+M[7]*C[68]+M[8]*C[69]+M[9]*C[70];
  L[26] += M[4]*C[62]+M[5]*C[66]+M[6]*C[67]+M[7]*C[71]+M[8]*C[72]+M[9]*C[73];
  L[27] += M[4]*C[63]+M[5]*C[67]+M[6]*C[68]+M[7]*C[72]+M[8]*C[73]+M[9]*C[74];
  L[28] += M[4]*C[64]+M[5]*C[68]+M[6]*C[69]+M[7]*C[73]+M[8]*C[74]+M[9]*C[75];
  L[29] += M[4]*C[65]+M[5]*C[69]+M[6]*C[70]+M[7]*C[74]+M[8]*C[75]+M[9]*C[76];
  L[30] += M[4]*C[66]+M[5]*C[71]+M[6]*C[72]+M[7]*C[77]+M[8]*C[78]+M[9]*C[79];
  L[31] += M[4]*C[67]+M[5]*C[72]+M[6]*C[73]+M[7]*C[78]+M[8]*C[79]+M[9]*C[80];
  L[32] += M[4]*C[68]+M[5]*C[73]+M[6]*C[74]+M[7]*C[79]+M[8]*C[80]+M[9]*C[81];
  L[33] += M[4]*C[69]+M[5]*C[74]+M[6]*C[75]+M[7]*C[80]+M[8]*C[81]+M[9]*C[82];
  L[34] += M[4]*C[70]+M[5]*C[75]+M[6]*C[76]+M[7]*C[81]+M[8]*C[82]+M[9]*C[83];
  L[35] += M[1]*C[56]+M[2]*C[57]+M[3]*C[58];
  L[36] += M[1]*C[57]+M[2]*C[59]+M[3]*C[60];
  L[37] += M[1]*C[58]+M[2]*C[60]+M[3]*C[61];
  L[38] += M[1]*C[59]+M[2]*C[62]+M[3]*C[63];
  L[39] += M[1]*C[60]+M[2]*C[63]+M[3]*C[64];
  L[40] += M[1]*C[61]+M[2]*C[64]+M[3]*C[65];
  L[41] += M[1]*C[62]+M[2]*C[66]+M[3]*C[67];
  L[42] += M[1]*C[63]+M[2]*C[67]+M[3]*C[68];
  L[43] += M[1]*C[64]+M[2]*C[68]+M[3]*C[69];
  L[44] += M[1]*C[65]+M[2]*C[69]+M[3]*C[70];
  L[45] += M[1]*C[66]+M[2]*C[71]+M[3]*C[72];
  L[46] += M[1]*C[67]+M[2]*C[72]+M[3]*C[73];
  L[47] += M[1]*C[68]+M[2]*C[73]+M[3]*C[74];
  L[48] += M[1]*C[69]+M[2]*C[74]+M[3]*C[75];
  L[49] += M[1]*C[70]+M[2]*C[75]+M[3]*C[76];
  L[50] += M[1]*C[71]+M[2]*C[77]+M[3]*C[78];
  L[51] += M[1]*C[72]+M[2]*C[78]+M[3]*C[79];
  L[52] += M[1]*C[73]+M[2]*C[79]+M[3]*C[80];
  L[53] += M[1]*C[74]+M[2]*C[80]+M[3]*C[81];
  L[54] += M[1]*C[75]+M[2]*C[81]+M[3]*C[82];
  L[55] += M[1]*C[76]+M[2]*C[82]+M[3]*C[83];
#endif
}

template<int PP>
inline void sumM2P(B_iter B, const Lset &C, const Mset &M) {
  B->TRG[0] += C[0];
  B->TRG[1] += C[1];
  B->TRG[2] += C[2];
  B->TRG[3] += C[3];
#if COMkernel
  for( int i=1; i<MTERM; ++i ) B->TRG[0] += M[i] * C[i+3];
#else
  for( int i=1; i<MTERM; ++i ) B->TRG[0] += M[i] * C[i];
#endif
  Kernels<0,0,1>::M2P(B,C,M);
}

template<>
inline void sumM2P<1>(B_iter B, const Lset &C, const Mset&) {
  B->TRG[0] += C[0];
  B->TRG[1] += C[1];
  B->TRG[2] += C[2];
  B->TRG[3] += C[3];
}

template<>
inline void sumM2P<2>(B_iter B, const Lset &C, const Mset &M) {
  sumM2P<1>(B,C,M);
#if COMkernel
#else
  B->TRG[0] += M[1]*C[1]+M[2]*C[2]+M[3]*C[3];
  B->TRG[1] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6];
  B->TRG[2] += M[1]*C[5]+M[2]*C[7]+M[3]*C[8];
  B->TRG[3] += M[1]*C[6]+M[2]*C[8]+M[3]*C[9];
#endif
}

template<>
inline void sumM2P<3>(B_iter B, const Lset &C, const Mset &M) {
  sumM2P<2>(B,C,M);
#if COMkernel
  B->TRG[0] += M[1]*C[4]+M[2]*C[5]+M[3]*C[6]+M[4]*C[7]+M[5]*C[8]+M[6]*C[9];
  B->TRG[1] += M[1]*C[10]+M[2]*C[11]+M[3]*C[12]+M[4]*C[13]+M[5]*C[14]+M[6]*C[15];
  B->TRG[2] += M[1]*C[11]+M[2]*C[13]+M[3]*C[14]+M[4]*C[16]+M[5]*C[17]+M[6]*C[18];
  B->TRG[3] += M[1]*C[12]+M[2]*C[14]+M[3]*C[15]+M[4]*C[17]+M[5]*C[18]+M[6]*C[19];
#else
  B->TRG[0] += M[4]*C[4]+M[5]*C[5]+M[6]*C[6]+M[7]*C[7]+M[8]*C[8]+M[9]*C[9];
  B->TRG[1] += M[4]*C[10]+M[5]*C[11]+M[6]*C[12]+M[7]*C[13]+M[8]*C[14]+M[9]*C[15];
  B->TRG[2] += M[4]*C[11]+M[5]*C[13]+M[6]*C[14]+M[7]*C[16]+M[8]*C[17]+M[9]*C[18];
  B->TRG[3] += M[4]*C[12]+M[5]*C[14]+M[6]*C[15]+M[7]*C[17]+M[8]*C[18]+M[9]*C[19];
#endif
}

template<>
inline void sumM2P<4>(B_iter B, const Lset &C, const Mset &M) {
  sumM2P<3>(B,C,M);
#if COMkernel
  B->TRG[0] += M[7]*C[10]+M[8]*C[11]+M[9]*C[12]+M[10]*C[13]+M[11]*C[14]+M[12]*C[15]+M[13]*C[16]+M[14]*C[17]+M[15]*C[18]+M[16]*C[19];
  B->TRG[1] += M[7]*C[20]+M[8]*C[21]+M[9]*C[22]+M[10]*C[23]+M[11]*C[24]+M[12]*C[25]+M[13]*C[26]+M[14]*C[27]+M[15]*C[28]+M[16]*C[29];
  B->TRG[2] += M[7]*C[21]+M[8]*C[23]+M[9]*C[24]+M[10]*C[26]+M[11]*C[27]+M[12]*C[28]+M[13]*C[30]+M[14]*C[31]+M[15]*C[32]+M[16]*C[33];
  B->TRG[3] += M[7]*C[22]+M[8]*C[24]+M[9]*C[25]+M[10]*C[27]+M[11]*C[28]+M[12]*C[29]+M[13]*C[31]+M[14]*C[32]+M[15]*C[33]+M[16]*C[34];
#else
  B->TRG[0] += M[10]*C[10]+M[11]*C[11]+M[12]*C[12]+M[13]*C[13]+M[14]*C[14]+M[15]*C[15]+M[16]*C[16]+M[17]*C[17]+M[18]*C[18]+M[19]*C[19];
  B->TRG[1] += M[10]*C[20]+M[11]*C[21]+M[12]*C[22]+M[13]*C[23]+M[14]*C[24]+M[15]*C[25]+M[16]*C[26]+M[17]*C[27]+M[18]*C[28]+M[19]*C[29];
  B->TRG[2] += M[10]*C[21]+M[11]*C[23]+M[12]*C[24]+M[13]*C[26]+M[14]*C[27]+M[15]*C[28]+M[16]*C[30]+M[17]*C[31]+M[18]*C[32]+M[19]*C[33];
  B->TRG[3] += M[10]*C[22]+M[11]*C[24]+M[12]*C[25]+M[13]*C[27]+M[14]*C[28]+M[15]*C[29]+M[16]*C[31]+M[17]*C[32]+M[18]*C[33]+M[19]*C[34];
#endif
}

template<>
inline void sumM2P<5>(B_iter B, const Lset &C, const Mset &M) {
  sumM2P<4>(B,C,M);
#if COMkernel
  B->TRG[0] += M[17]*C[20]+M[18]*C[21]+M[19]*C[22]+M[20]*C[23]+M[21]*C[24]+M[22]*C[25]+M[23]*C[26]+M[24]*C[27]+M[25]*C[28]+M[26]*C[29]+M[27]*C[30]+M[28]*C[31]+M[29]*C[32]+M[30]*C[33]+M[31]*C[34];
  B->TRG[1] += M[17]*C[35]+M[18]*C[36]+M[19]*C[37]+M[20]*C[38]+M[21]*C[39]+M[22]*C[40]+M[23]*C[41]+M[24]*C[42]+M[25]*C[43]+M[26]*C[44]+M[27]*C[45]+M[28]*C[46]+M[29]*C[47]+M[30]*C[48]+M[31]*C[49];
  B->TRG[2] += M[17]*C[36]+M[18]*C[38]+M[19]*C[39]+M[20]*C[41]+M[21]*C[42]+M[22]*C[43]+M[23]*C[45]+M[24]*C[46]+M[25]*C[47]+M[26]*C[48]+M[27]*C[50]+M[28]*C[51]+M[29]*C[52]+M[30]*C[53]+M[31]*C[54];
  B->TRG[3] += M[17]*C[37]+M[18]*C[39]+M[19]*C[40]+M[20]*C[42]+M[21]*C[43]+M[22]*C[44]+M[23]*C[46]+M[24]*C[47]+M[25]*C[48]+M[26]*C[49]+M[27]*C[51]+M[28]*C[52]+M[29]*C[53]+M[30]*C[54]+M[31]*C[55];
#else
  B->TRG[0] += M[20]*C[20]+M[21]*C[21]+M[22]*C[22]+M[23]*C[23]+M[24]*C[24]+M[25]*C[25]+M[26]*C[26]+M[27]*C[27]+M[28]*C[28]+M[29]*C[29]+M[30]*C[30]+M[31]*C[31]+M[32]*C[32]+M[33]*C[33]+M[34]*C[34];
  B->TRG[1] += M[20]*C[35]+M[21]*C[36]+M[22]*C[37]+M[23]*C[38]+M[24]*C[39]+M[25]*C[40]+M[26]*C[41]+M[27]*C[42]+M[28]*C[43]+M[29]*C[44]+M[30]*C[45]+M[31]*C[46]+M[32]*C[47]+M[33]*C[48]+M[34]*C[49];
  B->TRG[2] += M[20]*C[36]+M[21]*C[38]+M[22]*C[39]+M[23]*C[41]+M[24]*C[42]+M[25]*C[43]+M[26]*C[45]+M[27]*C[46]+M[28]*C[47]+M[29]*C[48]+M[30]*C[50]+M[31]*C[51]+M[32]*C[52]+M[33]*C[53]+M[34]*C[54];
  B->TRG[3] += M[20]*C[37]+M[21]*C[39]+M[22]*C[40]+M[23]*C[42]+M[24]*C[43]+M[25]*C[44]+M[26]*C[46]+M[27]*C[47]+M[28]*C[48]+M[29]*C[49]+M[30]*C[51]+M[31]*C[52]+M[32]*C[53]+M[33]*C[54]+M[34]*C[55];
#endif
}

template<>
inline void sumM2P<6>(B_iter B, const Lset &C, const Mset &M) {
  sumM2P<5>(B,C,M);
#if COMkernel
  B->TRG[0] += M[32]*C[35]+M[33]*C[36]+M[34]*C[37]+M[35]*C[38]+M[36]*C[39]+M[37]*C[40]+M[38]*C[41]+M[39]*C[42]+M[40]*C[43]+M[41]*C[44]+M[42]*C[45]+M[43]*C[46]+M[44]*C[47]+M[45]*C[48]+M[46]*C[49]+M[47]*C[50]+M[48]*C[51]+M[49]*C[52]+M[50]*C[53]+M[51]*C[54]+M[52]*C[55];
  B->TRG[1] += M[32]*C[56]+M[33]*C[57]+M[34]*C[58]+M[35]*C[59]+M[36]*C[60]+M[37]*C[61]+M[38]*C[62]+M[39]*C[63]+M[40]*C[64]+M[41]*C[65]+M[42]*C[66]+M[43]*C[67]+M[44]*C[68]+M[45]*C[69]+M[46]*C[70]+M[47]*C[71]+M[48]*C[72]+M[49]*C[73]+M[50]*C[74]+M[51]*C[75]+M[52]*C[76];
  B->TRG[2] += M[32]*C[57]+M[33]*C[59]+M[34]*C[60]+M[35]*C[62]+M[36]*C[63]+M[37]*C[64]+M[38]*C[66]+M[39]*C[67]+M[40]*C[68]+M[41]*C[69]+M[42]*C[71]+M[43]*C[72]+M[44]*C[73]+M[45]*C[74]+M[46]*C[75]+M[47]*C[77]+M[48]*C[78]+M[49]*C[79]+M[50]*C[80]+M[51]*C[81]+M[52]*C[82];
  B->TRG[3] += M[32]*C[58]+M[33]*C[60]+M[34]*C[61]+M[35]*C[63]+M[36]*C[64]+M[37]*C[65]+M[38]*C[67]+M[39]*C[68]+M[40]*C[69]+M[41]*C[70]+M[42]*C[72]+M[43]*C[73]+M[44]*C[74]+M[45]*C[75]+M[46]*C[76]+M[47]*C[78]+M[48]*C[79]+M[49]*C[80]+M[50]*C[81]+M[51]*C[82]+M[52]*C[83];
#else
  B->TRG[0] += M[35]*C[35]+M[36]*C[36]+M[37]*C[37]+M[38]*C[38]+M[39]*C[39]+M[40]*C[40]+M[41]*C[41]+M[42]*C[42]+M[43]*C[43]+M[44]*C[44]+M[45]*C[45]+M[46]*C[46]+M[47]*C[47]+M[48]*C[48]+M[49]*C[49]+M[50]*C[50]+M[51]*C[51]+M[52]*C[52]+M[53]*C[53]+M[54]*C[54]+M[55]*C[55];
  B->TRG[1] += M[35]*C[56]+M[36]*C[57]+M[37]*C[58]+M[38]*C[59]+M[39]*C[60]+M[40]*C[61]+M[41]*C[62]+M[42]*C[63]+M[43]*C[64]+M[44]*C[65]+M[45]*C[66]+M[46]*C[67]+M[47]*C[68]+M[48]*C[69]+M[49]*C[70]+M[50]*C[71]+M[51]*C[72]+M[52]*C[73]+M[53]*C[74]+M[54]*C[75]+M[55]*C[76];
  B->TRG[2] += M[35]*C[57]+M[36]*C[59]+M[37]*C[60]+M[38]*C[62]+M[39]*C[63]+M[40]*C[64]+M[41]*C[66]+M[42]*C[67]+M[43]*C[68]+M[44]*C[69]+M[45]*C[71]+M[46]*C[72]+M[47]*C[73]+M[48]*C[74]+M[49]*C[75]+M[50]*C[77]+M[51]*C[78]+M[52]*C[79]+M[53]*C[80]+M[54]*C[81]+M[55]*C[82];
  B->TRG[3] += M[35]*C[58]+M[36]*C[60]+M[37]*C[61]+M[38]*C[63]+M[39]*C[64]+M[40]*C[65]+M[41]*C[67]+M[42]*C[68]+M[43]*C[69]+M[44]*C[70]+M[45]*C[72]+M[46]*C[73]+M[47]*C[74]+M[48]*C[75]+M[49]*C[76]+M[50]*C[78]+M[51]*C[79]+M[52]*C[80]+M[53]*C[81]+M[54]*C[82]+M[55]*C[83];
#endif
}

class Kernel : public Sort {
protected:
  vec3   X0;
  real   R0;
  C_iter Ci0;
  C_iter Cj0;

private:
  inline void flipCoef(Lset &C) const {
    for( int i=1; i!=4; ++i ) C[i] = -C[i];
    for( int i=10; i!=20; ++i ) C[i] = -C[i];
  }

public:
  Kernel() : X0(0), R0(0) {}
  ~Kernel() {}

#ifdef SSE
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const {
    if( mutual ) {
      for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {
        real P0 = 0;
        vec3 F0 = 0;
        for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {
          vec3 dX = Bi->X - Bj->X - Xperiodic;
          real R2 = norm(dX) + EPS2;
          real invR2 = 1.0 / R2;
          if( R2 == 0 ) invR2 = 0;
          real invR = Bi->SRC * Bj->SRC * std::sqrt(invR2);
          dX *= invR2 * invR;
          P0 += invR;
          F0 += dX;
          Bj->TRG[0] += invR * mutual;
          Bj->TRG[1] += dX[0] * mutual;
          Bj->TRG[2] += dX[1] * mutual;
          Bj->TRG[3] += dX[2] * mutual;
        }
        Bi->TRG[0] += P0;
        Bi->TRG[1] -= F0[0];
        Bi->TRG[2] -= F0[1];
        Bi->TRG[3] -= F0[2];
      }
    } else {
      for( int bi=0; bi<Ci->NDLEAF; bi+=4 ) {
        B_iter Bi = Ci->LEAF + bi;
        int nvec = std::min(Ci->NDLEAF-bi,4);
        float16 target4;
        for( int i=0; i<nvec; i++ ) {
          target4.x[i] = Bi[i].X[0] - Xperiodic[0];
          target4.y[i] = Bi[i].X[1] - Xperiodic[1];
          target4.z[i] = Bi[i].X[2] - Xperiodic[2];
          target4.w[i] = EPS2;
        }
        B_iter Bj = Cj->LEAF;
        __m128 mi = _mm_load1_ps(&Bi->SRC);
        __m128 zero = _mm_setzero_ps();
        __m128 ax = _mm_setzero_ps();
        __m128 ay = _mm_setzero_ps();
        __m128 az = _mm_setzero_ps();
        __m128 phi = _mm_setzero_ps();

        __m128 xi = _mm_load_ps(target4.x);
        __m128 yi = _mm_load_ps(target4.y);
        __m128 zi = _mm_load_ps(target4.z);
        __m128 R2 = _mm_load_ps(target4.w);

        float4 source = make_float4(Bj);
        __m128 x2 = _mm_load1_ps(&source.x);
        x2 = _mm_sub_ps(x2, xi);
        __m128 y2 = _mm_load1_ps(&source.y);
        y2 = _mm_sub_ps(y2, yi);
        __m128 z2 = _mm_load1_ps(&source.z);
        z2 = _mm_sub_ps(z2, zi);
        __m128 mj = _mm_load1_ps(&source.w);

        __m128 xj = x2;
        x2 = _mm_mul_ps(x2, x2);
        R2 = _mm_add_ps(R2, x2);
        __m128 yj = y2;
        y2 = _mm_mul_ps(y2, y2);
        R2 = _mm_add_ps(R2, y2);
        __m128 zj = z2;
        z2 = _mm_mul_ps(z2, z2);
        R2 = _mm_add_ps(R2, z2);

        source = make_float4(++Bj);
        x2 = _mm_load_ps(&source.x);
        y2 = x2;
        z2 = x2;
        for( int j=0; j<Cj->NDLEAF; j++ ) {
          __m128 invR = _mm_rsqrt_ps(R2);
          __m128 mask = _mm_cmpgt_ps(R2,zero);
          invR = _mm_and_ps(invR,mask);
          R2 = _mm_load_ps(target4.w);
          x2 = _mm_shuffle_ps(x2, x2, 0x00);
          x2 = _mm_sub_ps(x2, xi);
          y2 = _mm_shuffle_ps(y2, y2, 0x55);
          y2 = _mm_sub_ps(y2, yi);
          z2 = _mm_shuffle_ps(z2, z2, 0xaa);
          z2 = _mm_sub_ps(z2, zi);

          mj = _mm_mul_ps(mj, invR);
          mj = _mm_mul_ps(mj, mi);
          phi = _mm_add_ps(phi, mj);
          invR = _mm_mul_ps(invR, invR);
          invR = _mm_mul_ps(invR, mj);
          mj = _mm_load_ps(&source.x);
          mj = _mm_shuffle_ps(mj, mj, 0xff);

          source = make_float4(++Bj);
          xj = _mm_mul_ps(xj, invR);
          ax = _mm_add_ps(ax, xj);
          xj = x2;
          x2 = _mm_mul_ps(x2, x2);
          R2 = _mm_add_ps(R2, x2);
          x2 = _mm_load_ps(&source.x);

          yj = _mm_mul_ps(yj, invR);
          ay = _mm_add_ps(ay, yj);
          yj = y2;
          y2 = _mm_mul_ps(y2, y2);
          R2 = _mm_add_ps(R2, y2);
          y2 = x2;

          zj = _mm_mul_ps(zj, invR);
          az = _mm_add_ps(az, zj);
          zj = z2;
          z2 = _mm_mul_ps(z2, z2);
          R2 = _mm_add_ps(R2, z2);
          z2 = x2;
        }
        _mm_store_ps(target4.x, ax);
        _mm_store_ps(target4.y, ay);
        _mm_store_ps(target4.z, az);
        _mm_store_ps(target4.w, phi);
        for( int i=0; i<nvec; i++ ) {
          Bi[i].TRG[0] += target4.w[i];
          Bi[i].TRG[1] += target4.x[i];
          Bi[i].TRG[2] += target4.y[i];
          Bi[i].TRG[3] += target4.z[i];
        }
      }
    }
  }
#else
  void P2P(C_iter Ci, C_iter Cj, bool mutual=true) const {
    for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {
      real P0 = 0;
      vec3 F0 = 0;
      for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {
        vec3 dX = Bi->X - Bj->X - Xperiodic;
        real R2 = norm(dX) + EPS2;
        real invR2 = 1.0 / R2;
        if( R2 == 0 ) invR2 = 0;
        real invR = Bi->SRC * Bj->SRC * std::sqrt(invR2);
        dX *= invR2 * invR;
        P0 += invR;
        F0 += dX;
        Bj->TRG[0] += invR * mutual;
        Bj->TRG[1] += dX[0] * mutual;
        Bj->TRG[2] += dX[1] * mutual;
        Bj->TRG[3] += dX[2] * mutual;
      }
      Bi->TRG[0] += P0;
      Bi->TRG[1] -= F0[0];
      Bi->TRG[2] -= F0[1];
      Bi->TRG[3] -= F0[2];
    }
  }
#endif

  void P2P(C_iter C) const {
    int NJ = C->NDLEAF;
    for( B_iter Bi=C->LEAF; Bi!=C->LEAF+C->NDLEAF; ++Bi, --NJ ) {
      real P0 = 0;
      vec3 F0 = 0;
      for( B_iter Bj=Bi+1; Bj!=Bi+NJ; ++Bj ) {
        vec3 dX = Bi->X - Bj->X;
        real R2 = norm(dX) + EPS2;
        real invR2 = 1.0 / R2;
        if( R2 == 0 ) invR2 = 0;
        real invR = Bi->SRC * Bj->SRC * std::sqrt(invR2);
        dX *= invR2 * invR;
        P0 += invR;
        F0 += dX;
        Bj->TRG[0] += invR;
        Bj->TRG[1] += dX[0];
        Bj->TRG[2] += dX[1];
        Bj->TRG[3] += dX[2];
      }
      Bi->TRG[0] += P0;
      Bi->TRG[1] -= F0[0];
      Bi->TRG[2] -= F0[1];
      Bi->TRG[3] -= F0[2];
    }
  }

  void P2M(C_iter C, real &Rmax) const {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vec3 dX = C->X - B->X;
      real R = std::sqrt(norm(dX));
      if( R > Rmax ) Rmax = R;
      Lset M;
      M[0] = B->SRC;
      Kernels<0,0,P-1>::power(M,dX);
#if COMkernel
      C->M[0] += M[0];
      for( int i=1; i<MTERM; ++i ) C->M[i] += M[i+3];
#else
      for( int i=0; i<MTERM; ++i ) C->M[i] += M[i];
#endif
    }
#if USE_RMAX
    C->RCRIT = std::min(C->R,Rmax);
#else
    C->RCRIT = C->R;
#endif
  }

  void M2M(C_iter Ci, real &Rmax) const {
    for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {
      vec3 dX = Ci->X - Cj->X;
      real R = std::sqrt(norm(dX)) + Cj->RCRIT;
      if( R > Rmax ) Rmax = R;
      Mset M;
      Lset C;
      C[0] = 1;
      Kernels<0,0,P-1>::power(C,dX);
      M = Cj->M;
#if COMkernel
      Ci->M[0] += C[0] * M[0];
      for( int i=1; i<MTERM; ++i ) Ci->M[i] += C[i+3] * M[0];
#else
      for( int i=0; i<MTERM; ++i ) Ci->M[i] += C[i] * M[0];
#endif
      Kernels<0,0,P-1>::M2M(Ci->M,C,M);
    }
#if USE_RMAX
    Ci->RCRIT = std::min(Ci->R,Rmax);
#else
    Ci->RCRIT = Ci->R;
#endif
  }

  void M2L(C_iter Ci, C_iter Cj, bool mutual=true) const {
    vec3 dX = Ci->X - Cj->X - Xperiodic;
    real invR2 = 1 / norm(dX);
    real invR  = Ci->M[0] * Cj->M[0] * std::sqrt(invR2);
    Lset C, L;
    getCoef<P>(C,dX,invR2,invR);
    sumM2L<P>(L,C,Cj->M);
    for( int i=0; i<LTERM; ++i ) {
//#pragma omp atomic
      Ci->L[i] += L[i];
    }
    if( mutual ) {
      flipCoef(C);
      sumM2L<P>(L,C,Ci->M);
      for( int i=0; i<LTERM; ++i ) {
//#pragma omp atomic
        Cj->L[i] += L[i];
      }
    }
  }

  void M2P(C_iter Ci, C_iter Cj, bool mutual=true) const {
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {
      vec3 dX = B->X - Cj->X - Xperiodic;
      real invR2 = 1 / norm(dX);
      real invR  = B->SRC * Cj->M[0] * std::sqrt(invR2);
      Lset C;
      getCoef<P>(C,dX,invR2,invR);
      sumM2P<P>(B,C,Cj->M);
    }
    if( mutual ) {
      for( B_iter B=Cj->LEAF; B!=Cj->LEAF+Cj->NDLEAF; ++B ) {
        vec3 dX = B->X - Ci->X + Xperiodic;
        real invR2 = 1 / norm(dX);
        real invR  = B->SRC * Ci->M[0] * std::sqrt(invR2);
        Lset C;
        getCoef<P>(C,dX,invR2,invR);
        sumM2P<P>(B,C,Ci->M);
      }
    }
  }

  void L2L(C_iter Ci) const {
    C_iter Cj = Ci0 + Ci->PARENT;
    vec3 dX = Ci->X - Cj->X;
    Lset C;
    C[0] = 1;
    Kernels<0,0,P>::power(C,dX);
    Ci->L /= Ci->M[0];
    Ci->L += Cj->L;
    for( int i=1; i<LTERM; ++i ) Ci->L[0] += C[i] * Cj->L[i];
    Kernels<0,0,P-1>::L2L(Ci->L,C,Cj->L);
  }

  void L2P(C_iter Ci) const {
    for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NCLEAF; ++B ) {
      vec3 dX = B->X - Ci->X;
      Lset C, L;
      C[0] = 1;
      Kernels<0,0,P>::power(C,dX);
      L = Ci->L;
      B->TRG /= B->SRC;
      B->TRG[0] += L[0];
      B->TRG[1] += L[1];
      B->TRG[2] += L[2];
      B->TRG[3] += L[3];
      for( int i=1; i<LTERM; ++i ) B->TRG[0] += C[i]*L[i];
      Kernels<0,0,1>::L2P(B,C,L);
    }
  }
};

#endif
