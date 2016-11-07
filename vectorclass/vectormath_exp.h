#ifndef VECTORMATH_EXP_H
#define VECTORMATH_EXP_H  1

#include "vectormath_common.h"
//#include "vectormath_trig.h"

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
#endif

template<class VTYPE, class VTYPE2, class BVTYPE, class BVTYPE2>
static inline VTYPE exp_f(VTYPE2 const & initial_x2) {
  const float p2 = 1.f/2.f;
  const float p3 = 1.f/6.f;
  const float p4 = 1.f/24.f;
  const float p5 = 1.f/120.f;
  const float p6 = 1.f/720.f;
  const float p7 = 1.f/5040.f;
  const float log2e = 1.44269504088896340736;
  VTYPE2  x, r, x2, z, n2;
  VTYPE initial_x = initial_x2;
  BVTYPE2 inrange2;
  BVTYPE inrange;
  float max_x;

  const float ln2f_hi  =  0.693359375f;
  const float ln2f_lo  = -2.12194440e-4f;
  max_x = 87.3f;
  x = initial_x;
  r = round(initial_x * log2e);
  x = nmul_add(r, VTYPE2(ln2f_hi), x);
  x = nmul_add(r, VTYPE2(ln2f_lo), x);

  x2 = x * x;
  z = polynomial_5(x,p2,p3,p4,p5,p6,p7);
  z = mul_add(z, x2, x);
  n2 = vec_pow2n(r);
  z = (z + 1.0f) * n2;
  inrange = vec_lt(vec_abs(initial_x),max_x);
  inrange = vec_and(inrange,is_finite(initial_x));
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = select(sign_bit(initial_x2), 0.f, vec_inf<VTYPE>());
    z = select(inrange2, z, r);
    z = select(is_nan(initial_x2), initial_x2, z);
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
static inline VTYPE exp_d(VTYPE2 const & initial_x2) {
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
  const double max_x = 708.39;
  const double log2e = 1.44269504088896340736;
  const double ln2d_hi = 0.693145751953125;
  const double ln2d_lo = 1.42860682030941723212E-6;
  VTYPE2  x, r, z, n2;
  VTYPE initial_x = initial_x2;
  BVTYPE inrange;

  x  = initial_x2;
  r  = round(initial_x2*log2e);
  x = nmul_add(r, ln2d_hi, x);
  x = nmul_add(r, ln2d_lo, x);

  z = polynomial_13m(x, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13);
  n2 = vec_pow2n(r);
  z = (z + 1.0) * n2;
  inrange  = abs(initial_x2) < max_x;
  inrange &= is_finite(initial_x2);
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = select(sign_bit(initial_x2), 0., infinite_vec<VTYPE2>());
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
