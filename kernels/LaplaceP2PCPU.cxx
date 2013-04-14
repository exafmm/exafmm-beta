#include "kernel.h"
template<typename T, int D, int N>
struct SIMD {
  static inline T setBody(B_iter B, int i) {
    T v;
    return v;
  }
  static inline T setIndex(int i) {
    T v;
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,8> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i  ].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D],
        B[i+4].X[D],B[i+5].X[D],B[i+6].X[D],B[i+7].X[D]);
    return v;
  }
  static inline T setIndex(int i) {
    T v(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7);
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,4> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D]);
    return v;
  }
  static inline T setIndex(int i) {
    T v(i,i+1,i+2,i+3);
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,2> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i].X[D],B[i+1].X[D]);
    return v;
  }
  static inline T setIndex(int i) {
    T v(i,i+1);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,8> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i  ].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC,
        B[i+4].SRC,B[i+5].SRC,B[i+6].SRC,B[i+7].SRC);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,4> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,2> {
  static inline T setBody(B_iter B, int i) {
    T v(B[i].SRC,B[i+1].SRC);
    return v;
  }
};

kreal_t transpose(ksimdvec v, int i) {
#if KAHAN
  kreal_t temp;
  temp.s = v.s[i];
  temp.c = v.c[i];
  return temp;
#else
  return v[i];
#endif
}

void Kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) const {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  int i = 0;
#if 1
  for ( ; i<=ni-NSIMD; i+=NSIMD) {
    simdvec zero = 0;
    ksimdvec pot = zero;
    ksimdvec ax = zero;
    ksimdvec ay = zero;
    ksimdvec az = zero;

    simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(Bi,i);
    simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(Bi,i);
    simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(Bi,i);
    simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(Bi,i);
    simdvec R2 = EPS2;

    simdvec xj = Xperiodic[0];
    xi -= xj;
    simdvec yj = Xperiodic[1];
    yi -= yj;
    simdvec zj = Xperiodic[2];
    zi -= zj;

    simdvec x2 = Bj[0].X[0];
    x2 -= xi;
    simdvec y2 = Bj[0].X[1];
    y2 -= yi;
    simdvec z2 = Bj[0].X[2];
    z2 -= zi;
    simdvec mj = Bj[0].SRC;

    xj = x2;
    R2 += x2 * x2;
    yj = y2;
    R2 += y2 * y2;
    zj = z2;
    R2 += z2 * z2;
    simdvec invR;

    x2 = Bj[1].X[0];
    y2 = Bj[1].X[1];
    z2 = Bj[1].X[2];
    for (int j=0; j<nj-2; j++) {
      invR = rsqrt(R2);
      invR &= R2 > zero;
      R2 = EPS2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      if (mutual) Bj[j].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = Bj[j+1].SRC;

      xj *= invR;
      ax += xj;
      if (mutual) Bj[j].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;
      x2 = Bj[j+2].X[0];

      yj *= invR;
      ay += yj;
      if (mutual) Bj[j].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;
      y2 = Bj[j+2].X[1];

      zj *= invR;
      az += zj;
      if (mutual) Bj[j].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
      z2 = Bj[j+2].X[2];
    }
    if ( nj > 1 ) {
      invR = rsqrt(R2);
      invR &= R2 > zero;
      R2 = EPS2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      if (mutual) Bj[nj-2].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = Bj[nj-1].SRC;

      xj *= invR;
      ax += xj;
      if (mutual) Bj[nj-2].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;

      yj *= invR;
      ay += yj;
      if (mutual) Bj[nj-2].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;

      zj *= invR;
      az += zj;
      if (mutual) Bj[nj-2].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
    }
    invR = rsqrt(R2);
    invR &= R2 > zero;
    mj *= invR * mi;
    pot += mj;
    if (mutual) Bj[nj-1].TRG[0] += sum(mj);
    invR = invR * invR * mj;

    xj *= invR;
    ax += xj;
    if (mutual) Bj[nj-1].TRG[1] -= sum(xj);
    yj *= invR;
    ay += yj;
    if (mutual) Bj[nj-1].TRG[2] -= sum(yj);
    zj *= invR;
    az += zj;
    if (mutual) Bj[nj-1].TRG[3] -= sum(zj);
    for (int k=0; k<NSIMD; k++) {
      Bi[i+k].TRG[0] += transpose(pot,k);
      Bi[i+k].TRG[1] += transpose(ax,k);
      Bi[i+k].TRG[2] += transpose(ay,k);
      Bi[i+k].TRG[3] += transpose(az,k);
    }
  }
#endif
  for ( ; i<ni; i++) {
    kreal_t pot = 0; 
    kreal_t ax = 0;
    kreal_t ay = 0;
    kreal_t az = 0;
    for (int j=0; j<nj; j++) {
      vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        ax += dX[0];
        ay += dX[1];
        az += dX[2];
        if (mutual) {
          Bj[j].TRG[0] += invR;
          Bj[j].TRG[1] += dX[0];
          Bj[j].TRG[2] += dX[1];
          Bj[j].TRG[3] += dX[2];
        }
      }
    }
    Bi[i].TRG[0] += pot;
    Bi[i].TRG[1] -= ax;
    Bi[i].TRG[2] -= ay;
    Bi[i].TRG[3] -= az;
  }
}

void Kernel::P2P(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  int i = 0;
#if 1
  for ( ; i<=n-NSIMD; i+=NSIMD) {
    simdvec zero = 0;
    ksimdvec pot = zero;
    ksimdvec ax = zero;
    ksimdvec ay = zero;
    ksimdvec az = zero;

    simdvec index = SIMD<simdvec,0,NSIMD>::setIndex(i);
    simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(B,i);
    simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(B,i);
    simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(B,i);
    simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(B,i);
    simdvec R2 = EPS2;

    simdvec xj = Xperiodic[0];
    xi -= xj;
    simdvec yj = Xperiodic[1];
    yi -= yj;
    simdvec zj = Xperiodic[2];
    zi -= zj;

    simdvec x2 = B[i+1].X[0];
    x2 -= xi;
    simdvec y2 = B[i+1].X[1];
    y2 -= yi;
    simdvec z2 = B[i+1].X[2];
    z2 -= zi;
    simdvec mj = B[i+1].SRC;

    xj = x2;
    R2 += x2 * x2;
    yj = y2;
    R2 += y2 * y2;
    zj = z2;
    R2 += z2 * z2;
    simdvec invR;

    x2 = B[i+2].X[0];
    y2 = B[i+2].X[1];
    z2 = B[i+2].X[2];
    for (int j=i+1; j<n-2; j++) {
      invR = rsqrt(R2);
      invR &= index < j;
      invR &= R2 > zero;
      R2 = EPS2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      B[j].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = B[j+1].SRC;

      xj *= invR;
      ax += xj;
      B[j].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;
      x2 = B[j+2].X[0];

      yj *= invR;
      ay += yj;
      B[j].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;
      y2 = B[j+2].X[1];

      zj *= invR;
      az += zj;
      B[j].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
      z2 = B[j+2].X[2];
    }
    if ( n-2 > i ) {
      invR = rsqrt(R2);
      invR &= index < n-2;
      invR &= R2 > zero;
      R2 = EPS2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      B[n-2].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = B[n-1].SRC;

      xj *= invR;
      ax += xj;
      B[n-2].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;

      yj *= invR;
      ay += yj;
      B[n-2].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;

      zj *= invR;
      az += zj;
      B[n-2].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
    }
    invR = rsqrt(R2);
    invR &= index < n-1;
    invR &= R2 > zero;
    mj *= invR * mi;
    pot += mj;
    B[n-1].TRG[0] += sum(mj);
    invR = invR * invR * mj;

    xj *= invR;
    ax += xj;
    B[n-1].TRG[1] -= sum(xj);
    yj *= invR;
    ay += yj;
    B[n-1].TRG[2] -= sum(yj);
    zj *= invR;
    az += zj;
    B[n-1].TRG[3] -= sum(zj);
    for (int k=0; k<NSIMD; k++) {
      B[i+k].TRG[0] += transpose(pot,k);
      B[i+k].TRG[1] += transpose(ax,k);
      B[i+k].TRG[2] += transpose(ay,k);
      B[i+k].TRG[3] += transpose(az,k);
    }
  }
#endif
  for ( ; i<n; i++) {
    kreal_t pot = 0;
    kreal_t ax = 0;
    kreal_t ay = 0;
    kreal_t az = 0;
    for (int j=i+1; j<n; j++) {
      vec3 dX = B[i].X - B[j].X;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        ax += dX[0];
        ay += dX[1];
        az += dX[2];
        B[j].TRG[0] += invR;
        B[j].TRG[1] += dX[0];
        B[j].TRG[2] += dX[1];
        B[j].TRG[3] += dX[2];
      }
    }
    B[i].TRG[0] += pot;
    B[i].TRG[1] -= ax;
    B[i].TRG[2] -= ay;
    B[i].TRG[3] -= az;
  }
}
