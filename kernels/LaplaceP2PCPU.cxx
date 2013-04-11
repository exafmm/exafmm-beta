#include "kernel.h"

#if __AVX__
const int NSIMD = 32 / sizeof(real_t);
typedef vec<NSIMD,real_t> simdvec;
inline float vecSum8s(simdvec v) {
  return v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7];
}
inline double vecSum4d(__m256d reg) {
  double mem[4] __attribute__ ((aligned(32)));
  _mm256_store_pd(mem, reg);
  return mem[0] + mem[1] + mem[2] + mem[3];
}
#elif __SSE__
const int NSIMD = 16 / sizeof(real_t);
typedef vec<NSIMD,real_t> simdvec;
inline float vecSum4s(__m128 reg) {
  float mem[4] __attribute__ ((aligned(16)));
  _mm_store_ps(mem, reg);
  return mem[0] + mem[1] + mem[2] + mem[3];
}
inline double vecSum2d(__m128d reg) {
  double mem[2] __attribute__ ((aligned(16)));
  _mm_store_pd(mem, reg);
  return mem[0] + mem[1];
}
#endif

template<typename T, int D, int N>
struct SIMD {
  static inline T setVec(B_iter B, int i) {
    T v;
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,8> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i  ].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D],
        B[i+4].X[D],B[i+5].X[D],B[i+6].X[D],B[i+7].X[D]);
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,4> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D]);
    return v;
  }
};
template<typename T, int D>
struct SIMD<T,D,2> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i].X[D],B[i+1].X[D]);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,8> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i  ].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC,
        B[i+4].SRC,B[i+5].SRC,B[i+6].SRC,B[i+7].SRC);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,4> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC);
    return v;
  }
};
template<typename T>
struct SIMD<T,3,2> {
  static inline T setVec(B_iter B, int i) {
    T v(B[i].SRC,B[i+1].SRC);
    return v;
  }
};


#if KAHAN >= KAHAN_IN_DIRECT

inline void accum(real_t & s, real_t ds, real_t & c) {
  // s += ds;
  real_t y = ds - c;
  real_t t = s + y;
  c = (t - s) - y;
  s = t;
}

inline void accumVec(vec3 & s, vec3 ds, vec3 & c) {
  // s += ds;
  vec3 y = ds - c;
  vec3 t = s + y;
  c = (t - s) - y;
  s = t;
}

void Kernel::P2PKahan(C_iter Ci, C_iter Cj, bool mutual) const {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  int i = 0;
  for ( ; i<ni; i++) {
    real_t pot = 0;
    real_t pot_c = 0;
    vec3 acc = 0;
    vec3 acc_c = 0;
    for (int j=0; j<nj; j++) {
      vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0f / R2;
        real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
	accum(pot, invR, pot_c);
	accumVec(acc, dX, acc_c);
        if (mutual) {
	  accum(Bj[j].TRG[0], invR,  Bj[j].TRGc[0]);
          accum(Bj[j].TRG[1], dX[0], Bj[j].TRGc[1]);
          accum(Bj[j].TRG[2], dX[1], Bj[j].TRGc[2]);
          accum(Bj[j].TRG[3], dX[2], Bj[j].TRGc[3]);
        }
      }
    }
    accum(Bi[i].TRG[0], pot,       Bi[i].TRGc[0]);
    accum(Bi[i].TRG[0], pot_c,     Bi[i].TRGc[0]);
    accum(Bi[i].TRG[1], -acc[0],   Bi[i].TRGc[1]);
    accum(Bi[i].TRG[1], -acc_c[0], Bi[i].TRGc[1]);
    accum(Bi[i].TRG[2], -acc[1],   Bi[i].TRGc[2]);
    accum(Bi[i].TRG[2], -acc_c[1], Bi[i].TRGc[2]);
    accum(Bi[i].TRG[3], -acc[2],   Bi[i].TRGc[3]);
    accum(Bi[i].TRG[3], -acc_c[2], Bi[i].TRGc[3]);
  }
}

#endif

#if KAHAN >= KAHAN_ALWAYS

void Kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) const {
  P2PKahan(Ci, Cj, mutual);
}

#else

void Kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) const {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  int i = 0;
  for ( ; i<=ni-NSIMD; i+=NSIMD) {
    simdvec zero = 0;
    simdvec pot = zero;
    simdvec ax = zero;
    simdvec ay = zero;
    simdvec az = zero;

    simdvec xi = SIMD<simdvec,0,NSIMD>::setVec(Bi,i);
    simdvec yi = SIMD<simdvec,1,NSIMD>::setVec(Bi,i);
    simdvec zi = SIMD<simdvec,2,NSIMD>::setVec(Bi,i);
    simdvec mi = SIMD<simdvec,3,NSIMD>::setVec(Bi,i);
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
      Bi[i+k].TRG[0] += pot[k];
      Bi[i+k].TRG[1] += ax[k];
      Bi[i+k].TRG[2] += ay[k];
      Bi[i+k].TRG[3] += az[k];
    }
  }
  for ( ; i<ni; i++) {
    real_t pot = 0;
    vec3 acc = 0;
    for (int j=0; j<nj; j++) {
      vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
        if (mutual) {
          Bj[j].TRG[0] += invR;
          Bj[j].TRG[1] += dX[0];
          Bj[j].TRG[2] += dX[1];
          Bj[j].TRG[3] += dX[2];
        }
      }
    }
    Bi[i].TRG[0] += pot;
    Bi[i].TRG[1] -= acc[0];
    Bi[i].TRG[2] -= acc[1];
    Bi[i].TRG[3] -= acc[2];
  }
}
#endif

#if KAHAN >= KAHAN_IN_DIRECT

void Kernel::P2PKahan(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  int i = 0;
  for ( ; i<n; i++) {
    real_t pot = 0;
    real_t pot_c = 0;
    vec3 acc = 0;
    vec3 acc_c = 0;
    for (int j=i+1; j<n; j++) {
      vec3 dX = B[i].X - B[j].X;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        accum(pot, invR, pot_c);
        accumVec(acc, dX, acc_c);
        accum(B[j].TRG[0], invR,  B[j].TRGc[0]);
        accum(B[j].TRG[1], dX[0], B[j].TRGc[1]);
        accum(B[j].TRG[2], dX[1], B[j].TRGc[2]);
        accum(B[j].TRG[3], dX[2], B[j].TRGc[3]);
      }
    }
    accum(B[i].TRG[0], pot,       B[i].TRGc[0]);
    accum(B[i].TRG[0], pot_c,     B[i].TRGc[0]);
    accum(B[i].TRG[1], -acc[0],   B[i].TRGc[1]);
    accum(B[i].TRG[1], -acc_c[0], B[i].TRGc[1]);
    accum(B[i].TRG[2], -acc[1],   B[i].TRGc[2]);
    accum(B[i].TRG[2], -acc_c[1], B[i].TRGc[2]);
    accum(B[i].TRG[3], -acc[2],   B[i].TRGc[3]);
    accum(B[i].TRG[3], -acc_c[2], B[i].TRGc[3]);
  }
}

#endif

#if KAHAN >= KAHAN_ALWAYS

void Kernel::P2P(C_iter Ci) const {
  P2PKahan(Ci);
}

#else
void Kernel::P2P(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  int i = 0;
#if __AVX__ 

#if !defined(REAL_TYPE) || REAL_TYPE == REAL_TYPE_FLOAT
  for ( ; i<=n-8; i+=8) {
    __m256 pot = _mm256_setzero_ps();
    __m256 ax = _mm256_setzero_ps();
    __m256 ay = _mm256_setzero_ps();
    __m256 az = _mm256_setzero_ps();

    __m256 xi = _mm256_sub_ps(_mm256_setr_ps(B[i].X[0],B[i+1].X[1],B[i+2].X[0],B[i+3].X[0],
                                             B[i+4].X[0],B[i+5].X[0],B[i+6].X[0],B[i+7].X[0]),
                              _mm256_set1_ps(Xperiodic[0]));
    __m256 yi = _mm256_sub_ps(_mm256_setr_ps(B[i].X[1],B[i+1].X[1],B[i+2].X[1],B[i+3].X[1],
                                             B[i+4].X[1],B[i+5].X[1],B[i+6].X[1],B[i+7].X[1]),
                              _mm256_set1_ps(Xperiodic[1]));
    __m256 zi = _mm256_sub_ps(_mm256_setr_ps(B[i].X[2],B[i+1].X[2],B[i+2].X[2],B[i+3].X[2],
                                             B[i+4].X[2],B[i+5].X[2],B[i+6].X[2],B[i+7].X[2]),
                              _mm256_set1_ps(Xperiodic[2]));
    __m256 mi = _mm256_setr_ps(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC,
      B[i+4].SRC,B[i+5].SRC,B[i+6].SRC,B[i+7].SRC);
    __m256 R2 = _mm256_set1_ps(EPS2);

    __m256 x2 = _mm256_set1_ps(B[i+1].X[0]);
    x2 = _mm256_sub_ps(x2, xi);
    __m256 y2 = _mm256_set1_ps(B[i+1].X[1]);
    y2 = _mm256_sub_ps(y2, yi);
    __m256 z2 = _mm256_set1_ps(B[i+1].X[2]);
    z2 = _mm256_sub_ps(z2, zi);
    __m256 mj = _mm256_set1_ps(B[i+1].SRC);

    __m256 xj = x2;
    x2 = _mm256_mul_ps(x2, x2);
    R2 = _mm256_add_ps(R2, x2);
    __m256 yj = y2;
    y2 = _mm256_mul_ps(y2, y2);
    R2 = _mm256_add_ps(R2, y2);
    __m256 zj = z2;
    z2 = _mm256_mul_ps(z2, z2);
    R2 = _mm256_add_ps(R2, z2);
    __m256 invR, mask;

    x2 = _mm256_set1_ps(B[i+2].X[0]);
    y2 = _mm256_set1_ps(B[i+2].X[1]);
    z2 = _mm256_set1_ps(B[i+2].X[2]);
    for (int j=i+1; j<n-2; j++) {
      invR = _mm256_rsqrt_ps(R2);
      mask = _mm256_cmp_ps(_mm256_setr_ps(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7),
        _mm256_set1_ps(j), _CMP_LT_OQ);
      mask = _mm256_and_ps(mask, _mm256_cmp_ps(R2, _mm256_setzero_ps(), _CMP_GT_OQ));
      invR = _mm256_and_ps(invR, mask);
      R2 = _mm256_set1_ps(EPS2);
      x2 = _mm256_sub_ps(x2, xi);
      y2 = _mm256_sub_ps(y2, yi);
      z2 = _mm256_sub_ps(z2, zi);

      mj = _mm256_mul_ps(mj, invR);
      mj = _mm256_mul_ps(mj, mi);
      pot = _mm256_add_ps(pot, mj);
      B[j].TRG[0] += vecSum8s(mj);
      invR = _mm256_mul_ps(invR, invR);
      invR = _mm256_mul_ps(invR, mj);
      mj = _mm256_set1_ps(B[j+1].SRC);

      xj = _mm256_mul_ps(xj, invR);
      ax = _mm256_add_ps(ax, xj);
      B[j].TRG[1] -= vecSum8s(xj);
      xj = x2;
      x2 = _mm256_mul_ps(x2, x2);
      R2 = _mm256_add_ps(R2, x2);
      x2 = _mm256_set1_ps(B[j+2].X[0]);

      yj = _mm256_mul_ps(yj, invR);
      ay = _mm256_add_ps(ay, yj);
      B[j].TRG[2] -= vecSum8s(yj);
      yj = y2;
      y2 = _mm256_mul_ps(y2, y2);
      R2 = _mm256_add_ps(R2, y2);
      y2 = _mm256_set1_ps(B[j+2].X[1]);

      zj = _mm256_mul_ps(zj, invR);
      az = _mm256_add_ps(az, zj);
      B[j].TRG[3] -= vecSum8s(zj);
      zj = z2;
      z2 = _mm256_mul_ps(z2, z2);
      R2 = _mm256_add_ps(R2, z2);
      z2 = _mm256_set1_ps(B[j+2].X[2]);
    }
    invR = _mm256_rsqrt_ps(R2);
    mask = _mm256_cmp_ps(_mm256_setr_ps(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7),
			 _mm256_set1_ps(n-2), _CMP_LT_OQ);
    mask = _mm256_and_ps(mask, _mm256_cmp_ps(R2, _mm256_setzero_ps(), _CMP_GT_OQ));
    invR = _mm256_and_ps(invR, mask);
    R2 = _mm256_set1_ps(EPS2);
    x2 = _mm256_sub_ps(x2, xi);
    y2 = _mm256_sub_ps(y2, yi);
    z2 = _mm256_sub_ps(z2, zi);

    mj = _mm256_mul_ps(mj, invR);
    mj = _mm256_mul_ps(mj, mi);
    pot = _mm256_add_ps(pot, mj);
    B[n-2].TRG[0] += vecSum8s(mj);
    invR = _mm256_mul_ps(invR, invR);
    invR = _mm256_mul_ps(invR, mj);
    mj = _mm256_set1_ps(B[n-1].SRC);

    xj = _mm256_mul_ps(xj, invR);
    ax = _mm256_add_ps(ax, xj);
    B[n-2].TRG[1] -= vecSum8s(xj);
    xj = x2;
    x2 = _mm256_mul_ps(x2, x2);
    R2 = _mm256_add_ps(R2, x2);

    yj = _mm256_mul_ps(yj, invR);
    ay = _mm256_add_ps(ay, yj);
    B[n-2].TRG[2] -= vecSum8s(yj);
    yj = y2;
    y2 = _mm256_mul_ps(y2, y2);
    R2 = _mm256_add_ps(R2, y2);

    zj = _mm256_mul_ps(zj, invR);
    az = _mm256_add_ps(az, zj);
    B[n-2].TRG[3] -= vecSum8s(zj);
    zj = z2;
    z2 = _mm256_mul_ps(z2, z2);
    R2 = _mm256_add_ps(R2, z2);

    invR = _mm256_rsqrt_ps(R2);
    mask = _mm256_cmp_ps(_mm256_setr_ps(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7),
			 _mm256_set1_ps(n-1), _CMP_LT_OQ);
    mask = _mm256_and_ps(mask, _mm256_cmp_ps(R2, _mm256_setzero_ps(), _CMP_GT_OQ));
    invR = _mm256_and_ps(invR, mask);
    mj = _mm256_mul_ps(mj, invR);
    mj = _mm256_mul_ps(mj, mi);
    pot = _mm256_add_ps(pot, mj);
    B[n-1].TRG[0] += vecSum8s(mj);
    invR = _mm256_mul_ps(invR, invR);
    invR = _mm256_mul_ps(invR, mj);

    xj = _mm256_mul_ps(xj, invR);
    ax = _mm256_add_ps(ax, xj);
    B[n-1].TRG[1] -= vecSum8s(xj);
    yj = _mm256_mul_ps(yj, invR);
    ay = _mm256_add_ps(ay, yj);
    B[n-1].TRG[2] -= vecSum8s(yj);
    zj = _mm256_mul_ps(zj, invR);
    az = _mm256_add_ps(az, zj);
    B[n-1].TRG[3] -= vecSum8s(zj);
    for (int k=0; k<8; k++) {
      B[i+k].TRG[0] += ((float*)&pot)[k];
      B[i+k].TRG[1] += ((float*)&ax)[k];
      B[i+k].TRG[2] += ((float*)&ay)[k];
      B[i+k].TRG[3] += ((float*)&az)[k];
    }
  }
#else  // P2P(C), AVX, double

  for ( ; i<=n-4; i+=4) {
    __m256d pot = _mm256_setzero_pd();
    __m256d ax = _mm256_setzero_pd();
    __m256d ay = _mm256_setzero_pd();
    __m256d az = _mm256_setzero_pd();

    __m256d xi = _mm256_sub_pd(_mm256_setr_pd(B[i].X[0],B[i+1].X[0],B[i+2].X[0],B[i+3].X[0]),
                               _mm256_set1_pd(Xperiodic[0]));
    __m256d yi = _mm256_sub_pd(_mm256_setr_pd(B[i].X[1],B[i+1].X[1],B[i+2].X[1],B[i+3].X[1]),
                               _mm256_set1_pd(Xperiodic[1]));
    __m256d zi = _mm256_sub_pd(_mm256_setr_pd(B[i].X[2],B[i+1].X[2],B[i+2].X[2],B[i+3].X[2]),
                               _mm256_set1_pd(Xperiodic[2]));
    __m256d mi = _mm256_setr_pd(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC);
    __m256d R2 = _mm256_set1_pd(EPS2);

    __m256d x2 = _mm256_set1_pd(B[i+1].X[0]);
    x2 = _mm256_sub_pd(x2, xi);
    __m256d y2 = _mm256_set1_pd(B[i+1].X[1]);
    y2 = _mm256_sub_pd(y2, yi);
    __m256d z2 = _mm256_set1_pd(B[i+1].X[2]);
    z2 = _mm256_sub_pd(z2, zi);
    __m256d mj = _mm256_set1_pd(B[i+1].SRC);

    __m256d xj = x2;
    x2 = _mm256_mul_pd(x2, x2);
    R2 = _mm256_add_pd(R2, x2);
    __m256d yj = y2;
    y2 = _mm256_mul_pd(y2, y2);
    R2 = _mm256_add_pd(R2, y2);
    __m256d zj = z2;
    z2 = _mm256_mul_pd(z2, z2);
    R2 = _mm256_add_pd(R2, z2);
    __m256d invR, mask;

    x2 = _mm256_set1_pd(B[i+2].X[0]);
    y2 = _mm256_set1_pd(B[i+2].X[1]);
    z2 = _mm256_set1_pd(B[i+2].X[2]);
    for (int j=i+1; j<n-2; j++) {
      invR = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_sqrt_pd(R2));
      mask = _mm256_cmp_pd(_mm256_setr_pd(i, i+1, i+2, i+3),
			   _mm256_set1_pd(j), _CMP_LT_OQ);
      mask = _mm256_and_pd(mask, _mm256_cmp_pd(R2, _mm256_setzero_pd(), _CMP_GT_OQ));
      invR = _mm256_and_pd(invR, mask);
      R2 = _mm256_set1_pd(EPS2);
      x2 = _mm256_sub_pd(x2, xi);
      y2 = _mm256_sub_pd(y2, yi);
      z2 = _mm256_sub_pd(z2, zi);

      mj = _mm256_mul_pd(mj, invR);
      mj = _mm256_mul_pd(mj, mi);
      pot = _mm256_add_pd(pot, mj);
      B[j].TRG[0] += vecSum4d(mj);
      invR = _mm256_mul_pd(invR, invR);
      invR = _mm256_mul_pd(invR, mj);
      mj = _mm256_set1_pd(B[j+1].SRC);

      xj = _mm256_mul_pd(xj, invR);
      ax = _mm256_add_pd(ax, xj);
      B[j].TRG[1] -= vecSum4d(xj);
      xj = x2;
      x2 = _mm256_mul_pd(x2, x2);
      R2 = _mm256_add_pd(R2, x2);
      x2 = _mm256_set1_pd(B[j+2].X[0]);

      yj = _mm256_mul_pd(yj, invR);
      ay = _mm256_add_pd(ay, yj);
      B[j].TRG[2] -= vecSum4d(yj);
      yj = y2;
      y2 = _mm256_mul_pd(y2, y2);
      R2 = _mm256_add_pd(R2, y2);
      y2 = _mm256_set1_pd(B[j+2].X[1]);

      zj = _mm256_mul_pd(zj, invR);
      az = _mm256_add_pd(az, zj);
      B[j].TRG[3] -= vecSum4d(zj);
      zj = z2;
      z2 = _mm256_mul_pd(z2, z2);
      R2 = _mm256_add_pd(R2, z2);
      z2 = _mm256_set1_pd(B[j+2].X[2]);
    }
    invR = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_sqrt_pd(R2));
    mask = _mm256_cmp_pd(_mm256_setr_pd(i, i+1, i+2, i+3),
			 _mm256_set1_pd(n-2), _CMP_LT_OQ);
    mask = _mm256_and_pd(mask, _mm256_cmp_pd(R2, _mm256_setzero_pd(), _CMP_GT_OQ));
    invR = _mm256_and_pd(invR, mask);
    R2 = _mm256_set1_pd(EPS2);
    x2 = _mm256_sub_pd(x2, xi);
    y2 = _mm256_sub_pd(y2, yi);
    z2 = _mm256_sub_pd(z2, zi);

    mj = _mm256_mul_pd(mj, invR);
    mj = _mm256_mul_pd(mj, mi);
    pot = _mm256_add_pd(pot, mj);
    B[n-2].TRG[0] += vecSum4d(mj);
    invR = _mm256_mul_pd(invR, invR);
    invR = _mm256_mul_pd(invR, mj);
    mj = _mm256_set1_pd(B[n-1].SRC);

    xj = _mm256_mul_pd(xj, invR);
    ax = _mm256_add_pd(ax, xj);
    B[n-2].TRG[1] -= vecSum4d(xj);
    xj = x2;
    x2 = _mm256_mul_pd(x2, x2);
    R2 = _mm256_add_pd(R2, x2);

    yj = _mm256_mul_pd(yj, invR);
    ay = _mm256_add_pd(ay, yj);
    B[n-2].TRG[2] -= vecSum4d(yj);
    yj = y2;
    y2 = _mm256_mul_pd(y2, y2);
    R2 = _mm256_add_pd(R2, y2);

    zj = _mm256_mul_pd(zj, invR);
    az = _mm256_add_pd(az, zj);
    B[n-2].TRG[3] -= vecSum4d(zj);
    zj = z2;
    z2 = _mm256_mul_pd(z2, z2);
    R2 = _mm256_add_pd(R2, z2);

    invR = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_sqrt_pd(R2));
    mask = _mm256_cmp_pd(_mm256_setr_pd(i, i+1, i+2, i+3),
			 _mm256_set1_pd(n-1), _CMP_LT_OQ);
    mask = _mm256_and_pd(mask, _mm256_cmp_pd(R2, _mm256_setzero_pd(), _CMP_GT_OQ));
    invR = _mm256_and_pd(invR, mask);
    mj = _mm256_mul_pd(mj, invR);
    mj = _mm256_mul_pd(mj, mi);
    pot = _mm256_add_pd(pot, mj);
    B[n-1].TRG[0] += vecSum4d(mj);
    invR = _mm256_mul_pd(invR, invR);
    invR = _mm256_mul_pd(invR, mj);

    xj = _mm256_mul_pd(xj, invR);
    ax = _mm256_add_pd(ax, xj);
    B[n-1].TRG[1] -= vecSum4d(xj);
    yj = _mm256_mul_pd(yj, invR);
    ay = _mm256_add_pd(ay, yj);
    B[n-1].TRG[2] -= vecSum4d(yj);
    zj = _mm256_mul_pd(zj, invR);
    az = _mm256_add_pd(az, zj);
    B[n-1].TRG[3] -= vecSum4d(zj);
    for (int k=0; k<4; k++) {
      B[i+k].TRG[0] += ((double*)&pot)[k];
      B[i+k].TRG[1] += ((double*)&ax)[k];
      B[i+k].TRG[2] += ((double*)&ay)[k];
      B[i+k].TRG[3] += ((double*)&az)[k];
    }
  }
#endif

#elif __SSE__ 
#if !defined(REAL_TYPE) || REAL_TYPE == REAL_TYPE_FLOAT
  for ( ; i<=n-4; i+=4) {
    __m128 pot = _mm_setzero_ps();
    __m128 ax = _mm_setzero_ps();
    __m128 ay = _mm_setzero_ps();
    __m128 az = _mm_setzero_ps();

    __m128 xi = _mm_sub_ps(_mm_setr_ps(B[i].X[0], B[i+1].X[0], B[i+2].X[0], B[i+3].X[0]),
                           _mm_load1_ps(&Xperiodic[0]));
    __m128 yi = _mm_sub_ps(_mm_setr_ps(B[i].X[1], B[i+1].X[1], B[i+2].X[1], B[i+3].X[1]),
                           _mm_load1_ps(&Xperiodic[1]));
    __m128 zi = _mm_sub_ps(_mm_setr_ps(B[i].X[2], B[i+1].X[2], B[i+2].X[2], B[i+3].X[2]),
                           _mm_load1_ps(&Xperiodic[2]));
    __m128 mi = _mm_setr_ps(B[i].SRC,  B[i+1].SRC,  B[i+2].SRC,  B[i+3].SRC);
    __m128 R2 = _mm_set1_ps(EPS2);

    __m128 x2 = _mm_load1_ps(&B[i+1].X[0]);
    x2 = _mm_sub_ps(x2, xi);
    __m128 y2 = _mm_load1_ps(&B[i+1].X[1]);
    y2 = _mm_sub_ps(y2, yi);
    __m128 z2 = _mm_load1_ps(&B[i+1].X[2]);
    z2 = _mm_sub_ps(z2, zi);
    __m128 mj = _mm_load1_ps(&B[i+1].SRC);

    __m128 xj = x2;
    x2 = _mm_mul_ps(x2, x2);
    R2 = _mm_add_ps(R2, x2);
    __m128 yj = y2;
    y2 = _mm_mul_ps(y2, y2);
    R2 = _mm_add_ps(R2, y2);
    __m128 zj = z2;
    z2 = _mm_mul_ps(z2, z2);
    R2 = _mm_add_ps(R2, z2);
    __m128 invR, mask;

    x2 = _mm_load_ps(&B[i+2].X[0]);
    y2 = x2;
    z2 = x2;
    for (int j=i+1; j<n-2; j++) {
      invR = _mm_rsqrt_ps(R2);
      mask = _mm_cmplt_ps(_mm_setr_ps(i, i+1, i+2, i+3), _mm_set1_ps(j));
      mask = _mm_and_ps(mask, _mm_cmpgt_ps(R2, _mm_setzero_ps()));
      invR = _mm_and_ps(invR, mask);
      R2 = _mm_set1_ps(EPS2);
      x2 = _mm_shuffle_ps(x2, x2, _MM_SHUFFLE(0,0,0,0));
      x2 = _mm_sub_ps(x2, xi);
      y2 = _mm_shuffle_ps(y2, y2, _MM_SHUFFLE(1,1,1,1));
      y2 = _mm_sub_ps(y2, yi);
      z2 = _mm_shuffle_ps(z2, z2, _MM_SHUFFLE(2,2,2,2));
      z2 = _mm_sub_ps(z2, zi);

      mj = _mm_mul_ps(mj, invR);
      mj = _mm_mul_ps(mj, mi);
      pot = _mm_add_ps(pot, mj);
      B[j].TRG[0] += vecSum4s(mj);
      invR = _mm_mul_ps(invR, invR);
      invR = _mm_mul_ps(invR, mj);
      mj = _mm_load_ps(&B[j+1].X[0]);
      mj = _mm_shuffle_ps(mj, mj, _MM_SHUFFLE(3,3,3,3));

      xj = _mm_mul_ps(xj, invR);
      ax = _mm_add_ps(ax, xj);
      B[j].TRG[1] -= vecSum4s(xj);
      xj = x2;
      x2 = _mm_mul_ps(x2, x2);
      R2 = _mm_add_ps(R2, x2);
      x2 = _mm_load_ps(&B[j+2].X[0]);

      yj = _mm_mul_ps(yj, invR);
      ay = _mm_add_ps(ay, yj);
      B[j].TRG[2] -= vecSum4s(yj);
      yj = y2;
      y2 = _mm_mul_ps(y2, y2);
      R2 = _mm_add_ps(R2, y2);
      y2 = x2;

      zj = _mm_mul_ps(zj, invR);
      az = _mm_add_ps(az, zj);
      B[j].TRG[3] -= vecSum4s(zj);
      zj = z2;
      z2 = _mm_mul_ps(z2, z2);
      R2 = _mm_add_ps(R2, z2);
      z2 = x2;
    }
    invR = _mm_rsqrt_ps(R2);
    mask = _mm_cmplt_ps(_mm_setr_ps(i, i+1, i+2, i+3), _mm_set1_ps(n-2));
    mask = _mm_and_ps(mask, _mm_cmpgt_ps(R2, _mm_setzero_ps()));
    invR = _mm_and_ps(invR, mask);
    R2 = _mm_set1_ps(EPS2);
    x2 = _mm_shuffle_ps(x2, x2, _MM_SHUFFLE(0,0,0,0));
    x2 = _mm_sub_ps(x2, xi);
    y2 = _mm_shuffle_ps(y2, y2, _MM_SHUFFLE(1,1,1,1));
    y2 = _mm_sub_ps(y2, yi);
    z2 = _mm_shuffle_ps(z2, z2, _MM_SHUFFLE(2,2,2,2));
    z2 = _mm_sub_ps(z2, zi);

    mj = _mm_mul_ps(mj, invR);
    mj = _mm_mul_ps(mj, mi);
    pot = _mm_add_ps(pot, mj);
    B[n-2].TRG[0] += vecSum4s(mj);
    invR = _mm_mul_ps(invR, invR);
    invR = _mm_mul_ps(invR, mj);
    mj = _mm_load_ps(&B[n-1].X[0]);
    mj = _mm_shuffle_ps(mj, mj, _MM_SHUFFLE(3,3,3,3));

    xj = _mm_mul_ps(xj, invR);
    ax = _mm_add_ps(ax, xj);
    B[n-2].TRG[1] -= vecSum4s(xj);
    xj = x2;
    x2 = _mm_mul_ps(x2, x2);
    R2 = _mm_add_ps(R2, x2);

    yj = _mm_mul_ps(yj, invR);
    ay = _mm_add_ps(ay, yj);
    B[n-2].TRG[2] -= vecSum4s(yj);
    yj = y2;
    y2 = _mm_mul_ps(y2, y2);
    R2 = _mm_add_ps(R2, y2);

    zj = _mm_mul_ps(zj, invR);
    az = _mm_add_ps(az, zj);
    B[n-2].TRG[3] -= vecSum4s(zj);
    zj = z2;
    z2 = _mm_mul_ps(z2, z2);
    R2 = _mm_add_ps(R2, z2);

    invR = _mm_rsqrt_ps(R2);
    mask = _mm_cmplt_ps(_mm_setr_ps(i, i+1, i+2, i+3), _mm_set1_ps(n-1));
    mask = _mm_and_ps(mask, _mm_cmpgt_ps(R2, _mm_setzero_ps()));
    invR = _mm_and_ps(invR, mask);
    mj = _mm_mul_ps(mj, invR);
    mj = _mm_mul_ps(mj, mi);
    pot = _mm_add_ps(pot, mj);
    B[n-1].TRG[0] += vecSum4s(mj);
    invR = _mm_mul_ps(invR, invR);
    invR = _mm_mul_ps(invR, mj);

    xj = _mm_mul_ps(xj, invR);
    ax = _mm_add_ps(ax, xj);
    B[n-1].TRG[1] -= vecSum4s(xj);
    yj = _mm_mul_ps(yj, invR);
    ay = _mm_add_ps(ay, yj);
    B[n-1].TRG[2] -= vecSum4s(yj);
    zj = _mm_mul_ps(zj, invR);
    az = _mm_add_ps(az, zj);
    B[n-1].TRG[3] -= vecSum4s(zj);
    for (int k=0; k<4; k++) {
      B[i+k].TRG[0] += ((float*)&pot)[k];
      B[i+k].TRG[1] += ((float*)&ax)[k];
      B[i+k].TRG[2] += ((float*)&ay)[k];
      B[i+k].TRG[3] += ((float*)&az)[k];
    }
  }

#else  // P2P(C), SSE, double

  for ( ; i<=n-2; i+=2) {
    __m128d pot = _mm_setzero_pd();
    __m128d ax = _mm_setzero_pd();
    __m128d ay = _mm_setzero_pd();
    __m128d az = _mm_setzero_pd();

    __m128d xi = _mm_sub_pd(_mm_setr_pd(B[i].X[0], B[i+1].X[0]),
                            _mm_load1_pd(&Xperiodic[0]));
    __m128d yi = _mm_sub_pd(_mm_setr_pd(B[i].X[1], B[i+1].X[1]),
                            _mm_load1_pd(&Xperiodic[1]));
    __m128d zi = _mm_sub_pd(_mm_setr_pd(B[i].X[2], B[i+1].X[2]),
                            _mm_load1_pd(&Xperiodic[2]));
    __m128d mi = _mm_setr_pd(B[i].SRC,  B[i+1].SRC);
    __m128d R2 = _mm_set1_pd(EPS2);

    __m128d x2 = _mm_load1_pd(&B[i+1].X[0]);
    x2 = _mm_sub_pd(x2, xi);
    __m128d y2 = _mm_load1_pd(&B[i+1].X[1]);
    y2 = _mm_sub_pd(y2, yi);
    __m128d z2 = _mm_load1_pd(&B[i+1].X[2]);
    z2 = _mm_sub_pd(z2, zi);
    __m128d mj = _mm_load1_pd(&B[i+1].SRC);

    __m128d xj = x2;
    x2 = _mm_mul_pd(x2, x2);
    R2 = _mm_add_pd(R2, x2);
    __m128d yj = y2;
    y2 = _mm_mul_pd(y2, y2);
    R2 = _mm_add_pd(R2, y2);
    __m128d zj = z2;
    z2 = _mm_mul_pd(z2, z2);
    R2 = _mm_add_pd(R2, z2);
    __m128d invR, mask;

    if ( n-2 > i ) {
      x2 = _mm_load_pd(&B[i+2].X[0]);
      y2 = x2;
      z2 = _mm_load_pd(&B[i+2].X[2]);
      for (int j=i+1; j<n-2; j++) {
	invR = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(R2));
	mask = _mm_cmplt_pd(_mm_setr_pd(i, i+1), _mm_set1_pd(j));
	mask = _mm_and_pd(mask, _mm_cmpgt_pd(R2, _mm_setzero_pd()));
	invR = _mm_and_pd(invR, mask);
	R2 = _mm_set1_pd(EPS2);
	x2 = _mm_shuffle_pd(x2, x2, 0x0);
	x2 = _mm_sub_pd(x2, xi);
	y2 = _mm_shuffle_pd(y2, y2, 0x3);
	y2 = _mm_sub_pd(y2, yi);
	z2 = _mm_shuffle_pd(z2, z2, 0x0);
	z2 = _mm_sub_pd(z2, zi);

	mj = _mm_mul_pd(mj, invR);
	mj = _mm_mul_pd(mj, mi);
	pot = _mm_add_pd(pot, mj);
	B[j].TRG[0] += vecSum2d(mj);
	invR = _mm_mul_pd(invR, invR);
	invR = _mm_mul_pd(invR, mj);
	mj = _mm_load_pd(&B[j+1].X[2]);
	mj = _mm_shuffle_pd(mj, mj, 0x3);

	xj = _mm_mul_pd(xj, invR);
	ax = _mm_add_pd(ax, xj);
	B[j].TRG[1] -= vecSum2d(xj);
	xj = x2;
	x2 = _mm_mul_pd(x2, x2);
	R2 = _mm_add_pd(R2, x2);
	x2 = _mm_load_pd(&B[j+2].X[0]);

	yj = _mm_mul_pd(yj, invR);
	ay = _mm_add_pd(ay, yj);
	B[j].TRG[2] -= vecSum2d(yj);
	yj = y2;
	y2 = _mm_mul_pd(y2, y2);
	R2 = _mm_add_pd(R2, y2);
	y2 = x2;

	zj = _mm_mul_pd(zj, invR);
	az = _mm_add_pd(az, zj);
	B[j].TRG[3] -= vecSum2d(zj);
	zj = z2;
	z2 = _mm_mul_pd(z2, z2);
	R2 = _mm_add_pd(R2, z2);
	z2 = _mm_load_pd(&B[j+2].X[2]);
      }
      invR = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(R2));
      mask = _mm_cmplt_pd(_mm_setr_pd(i, i+1), _mm_set1_pd(n-2));
      mask = _mm_and_pd(mask, _mm_cmpgt_pd(R2, _mm_setzero_pd()));
      invR = _mm_and_pd(invR, mask);
      R2 = _mm_set1_pd(EPS2);
      x2 = _mm_shuffle_pd(x2, x2, 0x0);
      x2 = _mm_sub_pd(x2, xi);
      y2 = _mm_shuffle_pd(y2, y2, 0x3);
      y2 = _mm_sub_pd(y2, yi);
      z2 = _mm_shuffle_pd(z2, z2, 0x0);
      z2 = _mm_sub_pd(z2, zi);

      mj = _mm_mul_pd(mj, invR);
      mj = _mm_mul_pd(mj, mi);
      pot = _mm_add_pd(pot, mj);
      B[n-2].TRG[0] += vecSum2d(mj);
      invR = _mm_mul_pd(invR, invR);
      invR = _mm_mul_pd(invR, mj);
      mj = _mm_load_pd(&B[n-1].X[2]);
      mj = _mm_shuffle_pd(mj, mj, 0x3);

      xj = _mm_mul_pd(xj, invR);
      ax = _mm_add_pd(ax, xj);
      B[n-2].TRG[1] -= vecSum2d(xj);
      xj = x2;
      x2 = _mm_mul_pd(x2, x2);
      R2 = _mm_add_pd(R2, x2);

      yj = _mm_mul_pd(yj, invR);
      ay = _mm_add_pd(ay, yj);
      B[n-2].TRG[2] -= vecSum2d(yj);
      yj = y2;
      y2 = _mm_mul_pd(y2, y2);
      R2 = _mm_add_pd(R2, y2);

      zj = _mm_mul_pd(zj, invR);
      az = _mm_add_pd(az, zj);
      B[n-2].TRG[3] -= vecSum2d(zj);
      zj = z2;
      z2 = _mm_mul_pd(z2, z2);
      R2 = _mm_add_pd(R2, z2);
    }

    invR = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(R2));
    mask = _mm_cmplt_pd(_mm_setr_pd(i, i+1), _mm_set1_pd(n-1));
    mask = _mm_and_pd(mask, _mm_cmpgt_pd(R2, _mm_setzero_pd()));
    invR = _mm_and_pd(invR, mask);
    mj = _mm_mul_pd(mj, invR);
    mj = _mm_mul_pd(mj, mi);
    pot = _mm_add_pd(pot, mj);
    B[n-1].TRG[0] += vecSum2d(mj);
    invR = _mm_mul_pd(invR, invR);
    invR = _mm_mul_pd(invR, mj);

    xj = _mm_mul_pd(xj, invR);
    ax = _mm_add_pd(ax, xj);
    B[n-1].TRG[1] -= vecSum2d(xj);
    yj = _mm_mul_pd(yj, invR);
    ay = _mm_add_pd(ay, yj);
    B[n-1].TRG[2] -= vecSum2d(yj);
    zj = _mm_mul_pd(zj, invR);
    az = _mm_add_pd(az, zj);
    B[n-1].TRG[3] -= vecSum2d(zj);
    for (int k=0; k<2; k++) {
      B[i+k].TRG[0] += ((double*)&pot)[k];
      B[i+k].TRG[1] += ((double*)&ax)[k];
      B[i+k].TRG[2] += ((double*)&ay)[k];
      B[i+k].TRG[3] += ((double*)&az)[k];
    }
  }
#endif

#endif // __SSE__

  for ( ; i<n; i++) {
    real_t pot = 0;
    vec3 acc = 0;
    for (int j=i+1; j<n; j++) {
      vec3 dX = B[i].X - B[j].X;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
        B[j].TRG[0] += invR;
        B[j].TRG[1] += dX[0];
        B[j].TRG[2] += dX[1];
        B[j].TRG[3] += dX[2];
      }
    }
    B[i].TRG[0] += pot;
    B[i].TRG[1] -= acc[0];
    B[i].TRG[2] -= acc[1];
    B[i].TRG[3] -= acc[2];
  }
}
#endif

