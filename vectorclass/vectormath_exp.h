#ifndef VECTORMATH_EXP_H
#define VECTORMATH_EXP_H  1

#include "vectormath_common.h"

static inline Vec2d vm_pow2n (Vec2d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec2d a = n + (bias + pow2_52);
  Vec2q b = reinterpret_i(a);
  Vec2q c = b << 52;
  Vec2d d = reinterpret_d(c);
  return d;
}

static inline Vec4f vm_pow2n (Vec4f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec4f a = n + (bias + pow2_23);
  Vec4i b = reinterpret_i(a);
  Vec4i c = b << 23;
  Vec4f d = reinterpret_f(c);
  return d;
}

#if __AVX__
static inline Vec4d vm_pow2n (Vec4d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec4d a = n + (bias + pow2_52);
  Vec4q b = reinterpret_i(a);
  Vec4q c = b << 52;
  Vec4d d = reinterpret_d(c);
  return d;
}

static inline Vec8f vm_pow2n (Vec8f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec8f a = n + (bias + pow2_23);
  Vec8i b = reinterpret_i(a);
  Vec8i c = b << 23;
  Vec8f d = reinterpret_f(c);
  return d;
}
#endif

#if __AVX512F__ | __MIC__
static inline Vec8d vm_pow2n (Vec8d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec8d a = n + (bias + pow2_52);
  Vec8q b = Vec8q(reinterpret_i(a));
  Vec8q c = b << 52;
  Vec8d d = Vec8d(reinterpret_d(c));
  return d;
}

static inline Vec16f vm_pow2n (Vec16f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec16f a = n + (bias + pow2_23);
  Vec16i b = Vec16i(reinterpret_i(a));
  Vec16i c = b << 23;
  Vec16f d = Vec16f(reinterpret_f(c));
  return d;
}
#endif

template<class VTYPE, class BVTYPE, int M1, int BA>
static inline VTYPE exp_f(VTYPE const & initial_x) {
  const float P0expf = 1.f/2.f;
  const float P1expf = 1.f/6.f;
  const float P2expf = 1.f/24.f;
  const float P3expf = 1.f/120.f;
  const float P4expf = 1.f/720.f;
  const float P5expf = 1.f/5040.f;
  VTYPE  x, r, x2, z, n2;
  BVTYPE inrange;
  float max_x;

  if (BA <= 1) {
    const float ln2f_hi  =  0.693359375f;
    const float ln2f_lo  = -2.12194440e-4f;
    max_x = (BA == 0) ? 87.3f : 89.0f;
    x = initial_x;
    r = round(initial_x*float(VM_LOG2E));
    x = nmul_add(r, VTYPE(ln2f_hi), x);
    x = nmul_add(r, VTYPE(ln2f_lo), x);
  } else if (BA == 2) {
    max_x = 126.f;
    r = round(initial_x);
    x = initial_x - r;
    x = x * (float)VM_LN2;
  } else if (BA == 10) {
    max_x = 37.9f;
    const float log10_2_hi = 0.301025391f;
    const float log10_2_lo = 4.60503907E-6f;
    x = initial_x;
    r = round(initial_x*float(VM_LOG2E*VM_LN10));
    x = nmul_add(r, VTYPE(log10_2_hi), x);
    x = nmul_add(r, VTYPE(log10_2_lo), x);
    x = x * (float)VM_LN10;
  } else {
    return 0.;
  }

  x2 = x * x;
  z = polynomial_5(x,P0expf,P1expf,P2expf,P3expf,P4expf,P5expf);
  z = mul_add(z, x2, x);
  if (BA == 1) r--;
  n2 = vm_pow2n(r);
  if (M1 == 0) {
    z = (z + 1.0f) * n2;
  } else {
    z = mul_add(z, n2, n2 - 1.0f);
  }
  inrange  = abs(initial_x) < max_x;
  inrange &= is_finite(initial_x);
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = select(sign_bit(initial_x), 0.f-M1, infinite_vec<VTYPE>());
    z = select(inrange, z, r);
    z = select(is_nan(initial_x), initial_x, z);
    return z;
  }
}

static inline __m128 _mm_exp_ps(__m128 const & x) {
  return __m128(exp_f<Vec4f, Vec4fb, 0, 0>(Vec4f(x)));
}
#ifdef __AVX__
static inline __m256 _mm256_exp_ps(__m256 const & x) {
  return __m256(exp_f<Vec8f, Vec8fb, 0, 0>(Vec8f(x)));
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512 _mm512_exp_ps(__m512 const & x) {
  return __m512(exp_f<Vec16f, Vec16fb, 0, 0>(Vec16f(x)));
}
#endif

template<class VTYPE, class BVTYPE, int M1, int BA>
static inline VTYPE exp_d(VTYPE const & initial_x) {
  const double p2  = 1./2.;
  const double p3  = 1./6.;
  const double p4  = 1./24.;
  const double p5  = 1./120.;
  const double p6  = 1./720.;
  const double p7  = 1./5040.;
  const double p8  = 1./40320.;
  const double p9  = 1./362880.;
  const double p10 = 1./3628800.;
  const double p11 = 1./39916800.;
  const double p12 = 1./479001600.;
  const double p13 = 1./6227020800.;
  double max_x;
  VTYPE  x, r, z, n2;
  BVTYPE inrange;

  if (BA <= 1) {
    max_x = BA == 0 ? 708.39 : 709.7;
    const double ln2d_hi = 0.693145751953125;
    const double ln2d_lo = 1.42860682030941723212E-6;
    x  = initial_x;
    r  = round(initial_x*VM_LOG2E);
    x = nmul_add(r, ln2d_hi, x);
    x = nmul_add(r, ln2d_lo, x);
  } else if (BA == 2) {
    max_x = 1022.0;
    r  = round(initial_x);
    x  = initial_x - r;
    x *= VM_LN2;
  } else if (BA == 10) {
    max_x = 307.65;
    const double log10_2_hi = 0.30102999554947019;
    const double log10_2_lo = 1.1451100899212592E-10;
    x  = initial_x;
    r  = round(initial_x*(VM_LOG2E*VM_LN10));
    x  = nmul_add(r, log10_2_hi, x);
    x  = nmul_add(r, log10_2_lo, x);
    x *= VM_LN10;
  } else  {
    return 0.;
  }
  z = polynomial_13m(x, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13);
  if (BA == 1) r--;
  n2 = vm_pow2n(r);
  if (M1 == 0) {
    z = (z + 1.0) * n2;
  } else {
    z = mul_add(z, n2, n2 - 1.0);
  }
  inrange  = abs(initial_x) < max_x;
  inrange &= is_finite(initial_x);
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = select(sign_bit(initial_x), 0.-M1, infinite_vec<VTYPE>());
    z = select(inrange, z, r);
    z = select(is_nan(initial_x), initial_x, z);
    return z;
  }
}

static inline __m128d _mm_exp_pd(__m128d const & x) {
  return __m128d(exp_d<Vec2d, Vec2db, 0, 0>(Vec2d(x)));
}
#ifdef __AVX__
static inline __m256d _mm256_exp_pd(__m256d const & x) {
  return __m256d(exp_d<Vec4d, Vec4db, 0, 0>(Vec4d(x)));
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512d _mm512_exp_pd(__m512d const & x) {
  return __m512d(exp_d<Vec8d, Vec8db, 0, 0>(Vec8d(x)));
}
#endif

#endif  // VECTORMATH_EXP_H
