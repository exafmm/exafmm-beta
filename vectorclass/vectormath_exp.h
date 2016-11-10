#ifndef VECTORMATH_EXP_H
#define VECTORMATH_EXP_H  1

#include "vectormath_trig.h"

static inline Vec2d vec_pow2n (Vec2d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec2d a = n + (bias + pow2_52);
  Vec2q b = reinterpret_i(a);
  Vec2q c = b << 52;
  Vec2d d = reinterpret_d(c);
  return d;
}
static inline Vec4f vec_pow2n (Vec4f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec4f a = n + (bias + pow2_23);
  Vec4i b = reinterpret_i(a);
  Vec4i c = b << 23;
  Vec4f d = reinterpret_f(c);
  return d;
}
template<class VTYPE>
static inline VTYPE vec_inf();
template<>
inline __m128 vec_inf() {
  return _mm_castsi128_ps(_mm_set1_epi32(0x7F800000));
}
template<>
inline __m128d vec_inf() {
  return _mm_castsi128_pd(_mm_setr_epi32(0,0x7FF00000,0,0x7FF00000));
}

#if __AVX__
static inline Vec4d vec_pow2n (Vec4d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec4d a = n + (bias + pow2_52);
  Vec4q b = reinterpret_i(a);
  Vec4q c = b << 52;
  Vec4d d = reinterpret_d(c);
  return d;
}
static inline Vec8f vec_pow2n (Vec8f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec8f a = n + (bias + pow2_23);
  Vec8i b = reinterpret_i(a);
  Vec8i c = b << 23;
  Vec8f d = reinterpret_f(c);
  return d;
}
template<>
inline __m256 vec_inf() {
  static const union {
    int i[8];
    __m256 ymm;
  } u = {{0x7F800000,0x7F800000,0x7F800000,0x7F800000,0x7F800000,0x7F800000,0x7F800000,0x7F800000}};
  return u.ymm;
}
template<>
inline __m256d vec_inf() {
  static const union {
    int i[8];
    __m256 ymm;
  } u = {{0,0x7FF00000,0,0x7FF00000,0,0x7FF00000,0,0x7FF00000}};
  return _mm256_castps_pd(u.ymm);
}
#endif

#if __AVX512F__ | __MIC__
static inline Vec8d vec_pow2n (Vec8d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  Vec8d a = n + (bias + pow2_52);
  Vec8q b = Vec8q(reinterpret_i(a));
  Vec8q c = b << 52;
  Vec8d d = Vec8d(reinterpret_d(c));
  return d;
}
static inline Vec16f vec_pow2n (Vec16f const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  Vec16f a = n + (bias + pow2_23);
  Vec16i b = Vec16i(reinterpret_i(a));
  Vec16i c = b << 23;
  Vec16f d = Vec16f(reinterpret_f(c));
  return d;
}
template<>
inline __m512 vec_inf() {
  union {
    int32_t i;
    float f;
  } u = {0x7F800000};
  return _mm512_set1_ps(u.f);
}
template<>
inline __m512d vec_inf() {
  union {
    uint64_t i;
    double f;
  } u = {0x7FF0000000000000};
  return _mm512_set1_pd(u.f);
}
#endif

template<class VTYPE, class VTYPE2, class BVTYPE, class BVTYPE2>
static inline VTYPE exp_f(VTYPE2 const & initial_x2) {
  const VTYPE zero = vec_set1_ps<VTYPE>(0.f);
  const VTYPE p2 = vec_set1_ps<VTYPE>(1.f/2.f);
  const VTYPE p3 = vec_set1_ps<VTYPE>(1.f/6.f);
  const VTYPE p4 = vec_set1_ps<VTYPE>(1.f/24.f);
  const VTYPE p5 = vec_set1_ps<VTYPE>(1.f/120.f);
  const VTYPE p6 = vec_set1_ps<VTYPE>(1.f/720.f);
  const VTYPE p7 = vec_set1_ps<VTYPE>(1.f/5040.f);
  const VTYPE log2e = vec_set1_ps<VTYPE>(1.44269504088896340736);
  const VTYPE ln2f_hi = vec_set1_ps<VTYPE>(-0.693359375f);
  const VTYPE ln2f_lo = vec_set1_ps<VTYPE>(2.12194440e-4f);
  const float max_x = 87.3f;
  VTYPE2 r, x, x2, x4, z, n2;
  VTYPE initial_x = initial_x2;
  BVTYPE2 inrange2;
  BVTYPE inrange;

  x = initial_x;
  r = vec_round(initial_x * log2e);
  x = mul_add(r, ln2f_hi, x);
  x = mul_add(r, ln2f_lo, x);

  x2 = x * x;
  x4 = x2 * x2;
  z = mul_add(mul_add(p5,x,p4), x2, mul_add(mul_add(p7,x,p6), x4, mul_add(p3,x,p2)));
  z = mul_add(z, x2, x);
  n2 = vec_pow2n(r);
  z = (z + 1.0f) * n2;
  inrange = vec_lt(vec_abs(initial_x),max_x);
  inrange = vec_and(inrange,is_finite(initial_x));
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = vec_select(sign_bit(initial_x), zero, vec_inf<VTYPE>());
    z = vec_select(inrange, z, r);
    z = vec_select(is_nan(initial_x), initial_x, z);
    return z;
  }
}

static inline __m128 _mm_exp_ps(__m128 const & x) {
  return exp_f<__m128, Vec4f, __m128i, Vec4fb>(Vec4f(x));
}
#ifdef __AVX__
static inline __m256 _mm256_exp_ps(__m256 const & x) {
  return exp_f<__m256, Vec8f, __m256i, Vec8fb>(Vec8f(x));
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512 _mm512_exp_ps(__m512 const & x) {
  return exp_f<__m512, Vec16f, __mmask16, Vec16fb>(Vec16f(x));
}
#endif

template<class VTYPE, class VTYPE2, class BVTYPE>
 inline VTYPE exp_d(VTYPE2 const & initial_x2) {
  const VTYPE zero = vec_set1_pd<VTYPE>(0.);
  const VTYPE p2  = vec_set1_pd<VTYPE>(1./2.);
  const VTYPE p3  = vec_set1_pd<VTYPE>(1./6.);
  const VTYPE p4  = vec_set1_pd<VTYPE>(1./24.);
  const VTYPE p5  = vec_set1_pd<VTYPE>(1./120.);
  const VTYPE p6  = vec_set1_pd<VTYPE>(1./720.);
  const VTYPE p7  = vec_set1_pd<VTYPE>(1./5040.);
  const VTYPE p8  = vec_set1_pd<VTYPE>(1./40320.);
  const VTYPE p9  = vec_set1_pd<VTYPE>(1./362880.);
  const VTYPE p10 = vec_set1_pd<VTYPE>(1./3628800.);
  const VTYPE p11 = vec_set1_pd<VTYPE>(1./39916800.);
  const VTYPE p12 = vec_set1_pd<VTYPE>(1./479001600.);
  const VTYPE p13 = vec_set1_pd<VTYPE>(1./6227020800.);
  const VTYPE log2e = vec_set1_pd<VTYPE>(1.44269504088896340736);
  const VTYPE ln2d_hi = vec_set1_pd<VTYPE>(-0.693145751953125);
  const VTYPE ln2d_lo = vec_set1_pd<VTYPE>(-1.42860682030941723212E-6);
  const double max_x = 708.39;
  VTYPE2  x, x2, x4, x8, r, z, n2;
  VTYPE initial_x = initial_x2;
  BVTYPE inrange;

  x = initial_x2;
  r = vec_round(initial_x2*log2e);
  x = mul_add(r, ln2d_hi, x);
  x = mul_add(r, ln2d_lo, x);

  x2 = x * x;
  x4 = x2 * x2;
  x8 = x4 * x4;
  z = mul_add(mul_add(mul_add(p13,x,p12), x4, mul_add(mul_add(p11,x,p10), x2, mul_add(p9,x,p8))), x8,
              mul_add(mul_add(mul_add(p7,x,p6), x2, mul_add(p5,x,p4)), x4, mul_add(mul_add(p3,x,p2),x2,x)));
  n2 = vec_pow2n(r);
  z = (z + 1.0) * n2;
  inrange = vec_lt(vec_abs(initial_x2),max_x);
  inrange &= is_finite(initial_x2);
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = select(sign_bit(initial_x2), zero, vec_inf<VTYPE>());
    z = select(inrange, z, r);
    z = select(is_nan(initial_x2), initial_x2, z);
    return z;
  }
}

static inline __m128d _mm_exp_pd(__m128d const & x) {
  return exp_d<__m128d, Vec2d, Vec2db>(Vec2d(x));
}
#ifdef __AVX__
static inline __m256d _mm256_exp_pd(__m256d const & x) {
  return exp_d<__m256d, Vec4d, Vec4db>(Vec4d(x));
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512d _mm512_exp_pd(__m512d const & x) {
  return exp_d<__m512d, Vec8d, Vec8db>(Vec8d(x));
}
#endif

#endif  // VECTORMATH_EXP_H
