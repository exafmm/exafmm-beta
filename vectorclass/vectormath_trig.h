#ifndef VECTORMATH_TRIG_H
#define VECTORMATH_TRIG_H  1
#include "vectorclass.h"

static inline __m128 vec_fmadd(__m128 const & a, __m128 const & b, __m128 const & c) {
#ifdef __FMA__
  return _mm_fmadd_ps(a, b, c);
#elif defined (__FMA4__)
  return _mm_macc_ps(a, b, c);
#else
  return a * b + c;
#endif
}
static inline __m128d vec_fmadd(__m128d const & a, __m128d const & b, __m128d const & c) {
#ifdef __FMA__
  return _mm_fmadd_pd(a, b, c);
#elif defined (__FMA4__)
  return _mm_macc_pd(a, b, c);
#else
  return a * b + c;
#endif
}
static inline __m128i vec_and(__m128i const & a, int const & b) {
  return _mm_and_si128(a,_mm_set1_epi32(b));
}
static inline __m128i vec_and(__m128i const & a, __m128i const & b) {
  return _mm_and_si128(a,b);
}
static inline __m128i vec_and_64(__m128i const & a, int64_t const & b) {
  return _mm_and_si128(a,_mm_set1_epi64x(b));
}
static inline __m128i vec_xor(__m128i const & a, __m128 const & b) {
  return _mm_xor_si128(a,_mm_castps_si128(b));
}
static inline __m128i vec_xor(__m128i const & a, __m128d const & b) {
  return _mm_xor_si128(a,_mm_castpd_si128(b));
}
static inline __m128 vec_xor(__m128 const & a, __m128i const & b) {
  return _mm_xor_ps(a,_mm_castsi128_ps(b));
}
static inline __m128d vec_xor(__m128d const & a, __m128i const & b) {
  return _mm_xor_pd(a,_mm_castsi128_pd(b));
}
static inline __m128i vec_neq(__m128i const & a, int const & b) {
  return _mm_xor_si128(_mm_cmpeq_epi32(a,_mm_set1_epi32(b)),_mm_set1_epi32(-1));
}
static inline __m128i vec_neq_64(__m128i const & a, int64_t const & b) {
  return _mm_xor_si128(_mm_cmpeq_epi64(a,_mm_set1_epi64x(b)),_mm_set1_epi64x(-1));
}
static inline __m128i vec_lt(__m128i const & a, int const & b) {
  return _mm_cmpgt_epi32(_mm_set1_epi32(b),a);
}
static inline __m128i vec_lt(__m128 const & a, float const & b) {
  return _mm_castps_si128(_mm_cmpgt_ps(_mm_set1_ps(b),a));
}
static inline __m128i vec_lt(__m128d const & a, double const & b) {
  return _mm_castpd_si128(_mm_cmpgt_pd(_mm_set1_pd(b),a));
}
static inline __m128i vec_sll(__m128i const & a, int const & b) {
  return _mm_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m128i vec_sll_64(__m128i const & a, int const & b) {
  return _mm_sll_epi64(a,_mm_cvtsi32_si128(b));
}
static inline __m128 vec_abs(__m128 const & a) {
  __m128 mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
  return _mm_and_ps(a,mask);
}
static inline __m128d vec_abs(__m128d const & a) {
  __m128d mask = _mm_castsi128_pd(_mm_setr_epi32(-1,0x7FFFFFFF,-1,0x7FFFFFFF));
  return _mm_and_pd(a,mask);
}
static inline __m128 vec_round(__m128 const & a) {
  return _mm_round_ps(a, 0);
}
static inline __m128d vec_round(__m128d const & a) {
  return _mm_round_pd(a, 0);
}
static inline __m128i vec_cvtps_epi32(__m128 const & a) {
    return _mm_cvttps_epi32(a);
}
static inline __m128i vec_cvtpd_epi32(__m128d const & a) {
  static const union {
    int i[4];
    __m128i xmm;
  } u = {{-1,-1,0,0}};
  __m128i a1 = _mm_and_si128(_mm_cvttpd_epi32(a), u.xmm);
  __m128i b1 = _mm_slli_si128(_mm_cvttpd_epi32(a), 8);
  return _mm_or_si128(a1,b1);
}
static inline __m128 vec_select(__m128i const & s, __m128 const & a, __m128 const & b) {
  return _mm_blendv_ps(b,a,_mm_castsi128_ps(s));
}
static inline __m128d vec_select(__m128i const & s, __m128d const & a, __m128d const & b) {
  return _mm_blendv_pd(b,a,_mm_castsi128_pd(s));
}
static inline __m128i is_finite(__m128 const & a) {
  __m128i t1 = _mm_castps_si128(a);
  __m128i t2 = _mm_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m128i t3 = _mm_set1_epi32(0xFF000000);
  __m128i t4 = _mm_and_si128(t2,t3);
  return _mm_xor_si128(_mm_cmpeq_epi32(t4,_mm_set1_epi32(0xFF000000)),_mm_set1_epi32(-1));
}
static inline __m128i is_finite(__m128d const & a) {
  __m128i t1 = _mm_castpd_si128(a);
  __m128i t2 = _mm_sll_epi64(t1,_mm_cvtsi32_si128(1));
  __m128i t3 = _mm_set1_epi64x(0xFFE0000000000000ll);
  __m128i t4 = _mm_and_si128(t2,t3);
  return _mm_xor_si128(_mm_cmpeq_epi64(t4,_mm_set1_epi64x(0xFFE0000000000000ll)),_mm_set1_epi64x(-1));
}
static inline __m128i is_nan(__m128 const & a) {
  __m128i t1 = _mm_castps_si128(a);
  __m128i t2 = _mm_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m128i t3 = _mm_set1_epi32(0xFF000000);
  __m128i t4 = _mm_and_si128(t2,t3);
  __m128i t5 = _mm_andnot_si128(t3,t2);
  t5 = _mm_xor_si128(_mm_cmpeq_epi32(t5,_mm_set1_epi32(0)),_mm_set1_epi32(-1));
  return _mm_and_si128(_mm_cmpeq_epi32(t4,t3),t5);
}
static inline __m128i is_nan(__m128d const & a) {
  __m128i t1 = _mm_castpd_si128(a);
  __m128i t2 = _mm_sll_epi64(t1,_mm_cvtsi32_si128(1));
  __m128i t3 = _mm_set1_epi64x(0xFFE0000000000000ll);
  __m128i t4 = _mm_and_si128(t2,t3);
  __m128i t5 = _mm_andnot_si128(t3,t2);
  t5 = _mm_xor_si128(_mm_cmpeq_epi64(t5,_mm_set1_epi64x(0)),_mm_set1_epi64x(-1));
  return _mm_and_si128(_mm_cmpeq_epi64(t4,t3),t5);
}
static inline bool horizontal_and(__m128i const & a) {
  static const union {
    int i[4];
    __m128i xmm;
  } u = {{-1,-1,-1,-1}};
  return _mm_testc_si128(a,u.xmm) != 0;
}
static inline bool horizontal_or(__m128i const & a) {
  return ! _mm_testz_si128(a,a);
}
static inline __m128i sign_bit(__m128 const & a) {
  __m128i t1 = _mm_castps_si128(a);
  return _mm_sra_epi32(t1,_mm_cvtsi32_si128(31));
}
static inline __m128i sign_bit(__m128d const & a) {
  __m128i t1 = _mm_castpd_si128(a);
  __m128i bm32 = _mm_cvtsi32_si128(31);
  __m128i sign = _mm_srai_epi32(t1,31);
  __m128i sra2 = _mm_sra_epi32(t1,bm32);
  __m128i sra3 = _mm_srli_epi64(sra2,32);
  __m128i mask = _mm_setr_epi32(0,-1,0,-1);
  return _mm_blendv_epi8(sra3,sign,mask);
}
static inline __m128 vec_pow2n (__m128 const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  __m128 a = n + (bias + pow2_23);
  __m128i b = _mm_castps_si128(a);
  __m128i c = _mm_sll_epi32(b,_mm_cvtsi32_si128(23));
  return _mm_castsi128_ps(c);
}
static inline __m128d vec_pow2n (__m128d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  __m128d a = n + (bias + pow2_52);
  __m128i b = _mm_castpd_si128(a);
  __m128i c = _mm_sll_epi64(b,_mm_cvtsi32_si128(52));
  return _mm_castsi128_pd(c);
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
template<class VTYPE>
static inline VTYPE vec_set1_ps(float const & a);
template<>
inline __m128 vec_set1_ps(float const & a) {
  return _mm_set1_ps(a);
}
template<class VTYPE>
static inline VTYPE vec_set1_pd(double const & a);
template<>
inline __m128d vec_set1_pd(double const & a) {
  return _mm_set1_pd(a);
}
static inline __m128 vec_cvtepi32_ps(__m128i const & a) {
  return _mm_cvtepi32_ps(a);
}
template<class VTYPE, class ITYPE>
static inline VTYPE vec_cvtepi32_pd(ITYPE const &);
template<>
inline __m128d vec_cvtepi32_pd<__m128d,__m128i>(__m128i const & a) {
  return _mm_cvtepi32_pd(a);
}
template<class ITYPE, class ITYPEH>
static inline ITYPE vec_cvtepi32_epi64(ITYPEH const &);
template<>
inline __m128i vec_cvtepi32_epi64<__m128i,__m128i>(__m128i const & a) {
  __m128i sign = _mm_srai_epi32(a,31);
  return _mm_unpacklo_epi32(a,sign);
}

#if __AVX__
static inline __m256 vec_fmadd(__m256 const & a, __m256 const & b, __m256 const & c) {
#ifdef __FMA__
  return _mm256_fmadd_ps(a, b, c);
#elif defined (__FMA4__)
  return _mm256_macc_ps(a, b, c);
#else
  return a * b + c;
#endif
}
static inline __m256d vec_fmadd(__m256d const & a, __m256d const & b, __m256d const & c) {
#ifdef __FMA__
  return _mm256_fmadd_pd(a, b, c);
#elif defined (__FMA4__)
  return _mm256_macc_pd(a, b, c);
#else
  return a * b + c;
#endif
}
static inline __m256i vec_and(__m256i const & a, int const & b) {
  return _mm256_and_si256(a,_mm256_set1_epi32(b));
}
static inline __m256i vec_and(__m256i const & a, __m256i const & b) {
  return _mm256_and_si256(a,b);
}
static inline __m256i vec_and_64(__m256i const & a, int64_t const & b) {
  return _mm256_and_si256(a,_mm256_set1_epi64x(b));
}
static inline __m256i vec_xor(__m256i const & a, __m256 const & b) {
  return _mm256_xor_si256(a,_mm256_castps_si256(b));
}
static inline __m256i vec_xor(__m256i const & a, __m256d const & b) {
  return _mm256_xor_si256(a,_mm256_castpd_si256(b));
}
static inline __m256 vec_xor(__m256 const & a, __m256i const & b) {
  return _mm256_xor_ps(a,_mm256_castsi256_ps(b));
}
static inline __m256d vec_xor(__m256d const & a, __m256i const & b) {
  return _mm256_xor_pd(a,_mm256_castsi256_pd(b));
}
static inline __m256i vec_neq(__m256i const & a, int const & b) {
  return _mm256_xor_si256(_mm256_cmpeq_epi32(a,_mm256_set1_epi32(b)), _mm256_set1_epi32(-1));
}
static inline __m256i vec_neq_64(__m256i const & a, int64_t const & b) {
  return _mm256_xor_si256(_mm256_cmpeq_epi64(a,_mm256_set1_epi64x(b)), _mm256_set1_epi64x(-1));
}
static inline __m256i vec_lt(__m256i const & a, int const & b) {
  return _mm256_cmpgt_epi32(_mm256_set1_epi32(b),a);
}
static inline __m256i vec_lt(__m256 const & a, float const & b) {
  return _mm256_castps_si256(_mm256_cmp_ps(a,_mm256_set1_ps(b),1));
}
static inline __m256i vec_lt(__m256d const & a, double const & b) {
  return _mm256_castpd_si256(_mm256_cmp_pd(a,_mm256_set1_pd(b),1));
}
static inline __m256i vec_sll(__m256i const & a, int const & b) {
  return _mm256_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m256i vec_sll_64(__m256i const & a, int const & b) {
  return _mm256_sll_epi64(a,_mm_cvtsi32_si128(b));
}
static inline __m256 vec_abs(__m256 const & a) {
  static const union {
    int i[8];
    __m256 ymm;
  } u = {{0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF,0x7FFFFFFF}};
  return _mm256_and_ps(a,u.ymm);
}
static inline __m256d vec_abs(__m256d const & a) {
  static const union {
    int i[8];
    __m256 ymm;
  } u = {{-1,0x7FFFFFFF,-1,0x7FFFFFFF,-1,0x7FFFFFFF,-1,0x7FFFFFFF}};
  return _mm256_and_pd(a,_mm256_castps_pd(u.ymm));
}
static inline __m256 vec_round(__m256 const & a) {
  return _mm256_round_ps(a, 0);
}
static inline __m256d vec_round(__m256d const & a) {
  return _mm256_round_pd(a, 0);
}
static inline __m256i vec_cvtps_epi32(__m256 const & a) {
  __m128i lo = _mm_cvttps_epi32(_mm256_castps256_ps128(a));
  __m128i hi = _mm_cvttps_epi32(_mm256_extractf128_ps(a,1));
  return _mm256_inserti128_si256(_mm256_castsi128_si256(lo),(hi),1);
}
static inline __m128i vec_cvtpd_epi32(__m256d const & x) {
  return _mm256_cvttpd_epi32(x);
}
static inline __m256 vec_select(__m256i const & s, __m256 const & a, __m256 const & b) {
  return _mm256_blendv_ps(b,a,_mm256_castsi256_ps(s));
}
static inline __m256d vec_select(__m256i const & s, __m256d const & a, __m256d const & b) {
  return _mm256_blendv_pd(b,a,_mm256_castsi256_pd(s));
}
static inline __m256i is_finite(__m256 const & a) {
  __m256i t1 = _mm256_castps_si256(a);
  __m256i t2 = _mm256_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m256i t3 = _mm256_set1_epi32(0xFF000000);
  __m256i t4 = _mm256_and_si256(t2,t3);
  return _mm256_xor_si256(_mm256_cmpeq_epi32(t4,_mm256_set1_epi32(0xFF000000)),_mm256_set1_epi32(-1));
}
static inline __m256i is_finite(__m256d const & a) {
  __m256i t1 = _mm256_castpd_si256(a);
  __m256i t2 = _mm256_sll_epi64(t1,_mm_cvtsi32_si128(1));
  __m256i t3 = _mm256_set1_epi64x(0xFFE0000000000000);
  __m256i t4 = _mm256_and_si256(t2,t3);
  return _mm256_xor_si256(_mm256_cmpeq_epi64(t4,_mm256_set1_epi64x(0xFFE0000000000000ll)), _mm256_set1_epi64x(-1));
}
static inline __m256i is_nan(__m256 const & a) {
  __m256i t1 = _mm256_castps_si256(a);
  __m256i t2 = _mm256_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m256i t3 = _mm256_set1_epi32(0xFF000000);
  __m256i t4 = _mm256_and_si256(t2,t3);
  __m256i t5 = _mm256_andnot_si256(t3,t2);
  t5 = _mm256_xor_si256(_mm256_cmpeq_epi32(t5,_mm256_set1_epi32(0)), _mm256_set1_epi32(-1));
  return _mm256_and_si256(_mm256_cmpeq_epi32(t4,t3),t5);
}
static inline __m256i is_nan(__m256d const & a) {
  __m256i t1 = _mm256_castpd_si256(a);
  __m256i t2 = _mm256_sll_epi64(t1,_mm_cvtsi32_si128(1));
  __m256i t3 = _mm256_set1_epi64x(0xFFE0000000000000);
  __m256i t4 = _mm256_and_si256(t2,t3);
  __m256i t5 = _mm256_andnot_si256(t3,t2);
  t5 = _mm256_xor_si256(_mm256_cmpeq_epi64(t5,_mm256_set1_epi64x(0)), _mm256_set1_epi64x(-1));
  return _mm256_and_si256(_mm256_cmpeq_epi64(t4,t3),t5);
}
static inline bool horizontal_and(__m256i const & a) {
  static const union {
    int32_t i[8];
    __m256i ymm;
  } u = {{-1,-1,-1,-1,-1,-1,-1,-1}};
  return _mm256_testc_si256(a,u.ymm) != 0;
}
static inline bool horizontal_or(__m256i const & a) {
  return ! _mm256_testz_si256(a,a);
}
static inline __m256i sign_bit(__m256 const & a) {
  __m256i t1 = _mm256_castps_si256(a);
  return _mm256_sra_epi32(t1, _mm_cvtsi32_si128(31));
}
static inline __m256i sign_bit(__m256d const & a) {
  static const union {
    int32_t i[8];
    __m256i ymm;
  } u = {{0,-1,0,-1,0,-1,0,-1}};
  __m256i t1 = _mm256_castpd_si256(a);
  __m128i bm32 = _mm_cvtsi32_si128(31);
  __m256i sign = _mm256_srai_epi32(t1,31);
  __m256i sra2 = _mm256_sra_epi32(t1,bm32);
  __m256i sra3 = _mm256_srli_epi64(sra2,32);
  return _mm256_blendv_epi8(sra3,sign,u.ymm);
}
static inline __m256 vec_pow2n (__m256 const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  __m256 a = n + (bias + pow2_23);
  __m256i b = _mm256_castps_si256(a);
  __m256i c = _mm256_sll_epi32(b,_mm_cvtsi32_si128(23));
  __m256 d = _mm256_castsi256_ps(c);
  return d;
}
static inline __m256d vec_pow2n (__m256d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  __m256d a = n + (bias + pow2_52);
  __m256i b = _mm256_castpd_si256(a);
  __m256i c = _mm256_sll_epi64(b,_mm_cvtsi32_si128(52));
  return _mm256_castsi256_pd(c);
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
inline __m256 vec_set1_ps(float const & a) {
  return _mm256_set1_ps(a);
}
template<>
inline __m256d vec_set1_pd(double const & a) {
  return _mm256_set1_pd(a);
}
static inline __m256 vec_cvtepi32_ps(__m256i const & a) {
  return _mm256_cvtepi32_ps(a);
}
template<>
inline __m256d vec_cvtepi32_pd<__m256d,__m128i>(__m128i const & a) {
  return _mm256_cvtepi32_pd(a);
}
template<>
inline __m256i vec_cvtepi32_epi64<__m256i,__m128i>(__m128i const & a) {
  __m256i a1 = _mm256_inserti128_si256(_mm256_castsi128_si256(a),a,1);
  __m256i a2 = _mm256_permute4x64_epi64(a1,16);
  __m256i sign = _mm256_srai_epi32(a2,31);
  return _mm256_unpacklo_epi32(a2,sign);
}
#endif

#if __AVX512F__ | __MIC__
static inline __m512 _mm512_castpd_ps(__m512d x) {
  union {
    __m512d a;
    __m512  b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512d _mm512_castps_pd(__m512 x) {
  union {
    __m512  a;
    __m512d b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512 _mm512_castps256_ps512(__m256 x) {
  union {
    __m256 a;
    __m512 b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512d _mm512_castpd256_pd512(__m256d x) {
  union {
    __m256d a;
    __m512d b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m256d _mm512_castpd512_pd256(__m512d x) {
  union {
    __m512d a;
    __m256d b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m256 _mm512_castps512_ps256(__m512 x) {
  union {
    __m512 a;
    __m256 b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512i _mm512_castps_si512(__m512 x) {
  union {
    __m512  a;
    __m512i b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512 _mm512_castsi512_ps(__m512i x) {
  union {
    __m512i a;
    __m512  b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512i _mm512_castpd_si512(__m512d x) {
  union {
    __m512d a;
    __m512i b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512d _mm512_castsi512_pd(__m512i x) {
  union {
    __m512i a;
    __m512d b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512i _mm512_castsi256_si512(__m256i x) {
  union {
    __m512i a;
    __m256i b;
  } u;
  u.b = x;
  return u.a;
}
static inline __m256i _mm512_castsi512_si256(__m512i x) {
  union {
    __m512i a;
    __m256i b;
  } u;
  u.a = x;
  return u.b;
}
static inline __m512 vec_fmadd(__m512 const & a, __m512 const & b, __m512 const & c) {
  return _mm512_fmadd_ps(a, b, c);
}
static inline __m512d vec_fmadd(__m512d const & a, __m512d const & b, __m512d const & c) {
  return _mm512_fmadd_pd(a, b, c);
}
static inline __m512i vec_and(__m512i const & a, int const & b) {
  return _mm512_and_epi32(a,_mm512_set1_epi32(b));
}
static inline __mmask8 vec_and(__mmask8 const & a, __mmask8 const & b) {
  return _mm512_kand(a,b);
}
static inline __mmask16 vec_and(__mmask16 const & a, __mmask16 const & b) {
  return _mm512_kand(a,b);
}
static inline __m512i vec_and_64(__m512i const & a, int64_t const & b) {
  return _mm512_and_epi64(a,_mm512_set1_epi64(b));
}
static inline __m512i vec_xor(__m512i const & a, __m512 const & b) {
  return _mm512_xor_epi32(a,_mm512_castps_si512(b));
}
static inline __m512i vec_xor(__m512i const & a, __m512d const & b) {
  return _mm512_xor_epi64(a,_mm512_castpd_si512(b));
}
static inline __m512 vec_xor(__m512 const & a, __m512i const & b) {
  return _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(a),b));
}
static inline __m512d vec_xor(__m512d const & a, __m512i const & b) {
  return _mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(a),b));
}
static inline __mmask16 vec_neq(__m512i const & a, int const & b) {
  return _mm512_cmpneq_epi32_mask(a,_mm512_set1_epi32(b));
}
static inline __mmask8 vec_neq_64(__m512i const & a, int64_t const & b) {
  return _mm512_cmpneq_epi64_mask(a,_mm512_set1_epi64(b));
}
static inline __mmask16 vec_lt(__m512i const & a, int const & b) {
  return _mm512_cmpgt_epi32_mask(_mm512_set1_epi32(b),a);
}
static inline __mmask16 vec_lt(__m512 const & a, float const & b) {
  return _mm512_cmp_ps_mask(a,_mm512_set1_ps(b),1);
}
static inline __mmask8 vec_lt(__m512d const & a, double const & b) {
  return _mm512_cmp_pd_mask(a,_mm512_set1_pd(b),1);
}
static inline __m512i vec_sll(__m512i const & a, int const & b) {
  return _mm512_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m512i vec_sll_64(__m512i const & a, int const & b) {
  return _mm512_sll_epi64(a,_mm_cvtsi32_si128(b));
}
static inline __m512 vec_abs(__m512 const & a) {
  union {
    int32_t i;
    float f;
  } u = {0x7FFFFFFF};
  union {
    __m512  a;
    __m512i b;
  } v;
  union {
    __m512  a;
    __m512i b;
  } w;
  v.a = _mm512_set1_ps(u.f);
  w.a = a;
  v.b = _mm512_and_epi32(w.b,v.b);
  return v.a;
}
static inline __m512d vec_abs(__m512d const & a) {
  union {
    __m512d a;
    __m512i b;
  } u;
  u.a = a;
  u.b = _mm512_and_epi64(u.b,_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF));
  return u.a;
}
static inline __m512 vec_round(__m512 const & a) {
  return _mm512_roundscale_ps(a, 0);
}
static inline __m512d vec_round(__m512d const & a) {
  return _mm512_roundscale_pd(a, 0);
}
static inline __m512i vec_cvtps_epi32(__m512 const & a) {
  return _mm512_cvtt_roundps_epi32(a, _MM_FROUND_NO_EXC);
}
static inline __m256i vec_cvtpd_epi32(__m512d const & a) {
  return _mm512_cvtpd_epi32(a);
}
static inline __m512 vec_select(__mmask16 const & s, __m512 const & a, __m512 const & b) {
  return _mm512_mask_mov_ps(b,s,a);
}
static inline __m512d vec_select(__mmask8 const & s, __m512d const & a, __m512d const & b) {
  return _mm512_mask_mov_pd(b,s,a);
}
static inline __mmask16 is_finite(__m512 const & a) {
  __m512i t1 = _mm512_castps_si512(a);
  __m512i t2 = _mm512_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m512i t3 = _mm512_set1_epi32(0xFF000000);
  __m512i t4 = _mm512_and_epi32(t2,t3);
  return _mm512_cmpneq_epi32_mask(t4,_mm512_set1_epi32(0xFF000000));
}
static inline __mmask8 is_finite(__m512d const & a) {
  __m512i t1 = _mm512_castpd_si512(a);
  __m512i t2 = _mm512_sll_epi64(t1,_mm_cvtsi32_si128(1));
  __m512i t3 = _mm512_set1_epi64(0xFFE0000000000000);
  __m512i t4 = _mm512_and_epi32(t2,t3);
  return _mm512_cmpneq_epi64_mask(t4,_mm512_set1_epi64(0xFFE0000000000000ll));
}
static inline __mmask16 is_nan(__m512 const & a) {
  __m512i t1 = _mm512_castps_si512(a);
  __m512i t2 = _mm512_sll_epi32(t1,_mm_cvtsi32_si128(1));
  __m512i t3 = _mm512_set1_epi32(0xFF000000);
  __m512i t4 = _mm512_and_epi32(t2,t3);
  __m512i t5 = _mm512_andnot_si512(t3,t2);
  return _mm512_kand(_mm512_cmpeq_epi32_mask(t4,t3),_mm512_cmpneq_epi32_mask(t5,_mm512_set1_epi32(0)));
}
static inline __mmask8 is_nan(__m512d const & a) {
  __m512i t1 = _mm512_castpd_si512(a);
  __m512i t2 = _mm512_sll_epi64(t1, _mm_cvtsi32_si128(1));
  __m512i t3 = _mm512_set1_epi64(0xFFE0000000000000);
  __m512i t4 = _mm512_and_epi32(t2,t3);
  __m512i t5 = _mm512_andnot_si512(t3,t2);
  return _mm512_kand(_mm512_cmpeq_epi64_mask(t4,t3),_mm512_cmpneq_epi64_mask(t5,_mm512_set1_epi64(0)));
}
static inline bool horizontal_and(__mmask16 const & a) {
  return (uint16_t)(a == 0xFFFF);
}
static inline bool horizontal_or(__mmask16 const & a) {
  return (uint16_t)(a != 0);
}
static inline __mmask16 sign_bit(__m512 const & a) {
  __m512i t1 = _mm512_castps_si512(a);
  return _mm512_cmp_ps_mask(t1,_mm512_set1_ps(0),1);
}
static inline __mmask8 sign_bit(__m512d const & a) {
  __m512i t1 = _mm512_castpd_si512(a);
  return _mm512_cmp_pd_mask(t1,_mm512_set1_pd(0),1);
}
static inline __m512 vec_pow2n (__m512 const & n) {
  const float pow2_23 =  8388608.0;
  const float bias = 127.0;
  __m512 a = n + (bias + pow2_23);
  __m512i b = _mm512_castps_si512(a);
  __m512i c = _mm512_sll_epi32(b,_mm_cvtsi32_si128(23));
  return _mm512_castsi512_ps(c);
}
static inline __m512d vec_pow2n (__m512d const & n) {
  const double pow2_52 = 4503599627370496.0;
  const double bias = 1023.0;
  __m512d a = n + (bias + pow2_52);
  __m512i b = _mm512_castpd_si512(a);
  __m512i c = _mm512_sll_epi64(b,_mm_cvtsi32_si128(52));
  return _mm512_castsi512_pd(c);
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
template<>
inline __m512 vec_set1_ps(float const & a) {
  return _mm512_set1_ps(a);
}
template<>
inline __m512d vec_set1_pd(double const & a) {
  return _mm512_set1_pd(a);
}
static inline __m512 vec_cvtepi32_ps(__m512i const & a) {
  return _mm512_cvtepi32_ps(a);
}
template<>
inline __m512d vec_cvtepi32_pd<__m512d,__m256i>(__m256i const & a) {
  return _mm512_cvtepi32_pd(a);
}
template<>
inline __m512i vec_cvtepi32_epi64<__m512i,__m256i>(__m256i const & a) {
  return _mm512_cvtepi32_epi64(_mm512_castsi512_si256(_mm512_inserti64x4(_mm512_castsi256_si512(a),a,1)));
}
#endif

template<class VTYPE, class ITYPE, class BVTYPE, int SC>
static inline VTYPE sincos_f(VTYPE * cosret, VTYPE const & xx) {
  const VTYPE zero = vec_set1_ps<VTYPE>(0.f);
  const VTYPE half = vec_set1_ps<VTYPE>(-.5f);
  const VTYPE one = vec_set1_ps<VTYPE>(1.f);
  const VTYPE pi4 = vec_set1_ps<VTYPE>(4./M_PI);
  const VTYPE d0 = vec_set1_ps<VTYPE>(0.78515625f);
  const VTYPE d1 = vec_set1_ps<VTYPE>(2.4187564849853515625E-4f);
  const VTYPE d2 = vec_set1_ps<VTYPE>(3.77489497744594108E-8f);
  const VTYPE s0 = vec_set1_ps<VTYPE>(-1.6666654611E-1f);
  const VTYPE s1 = vec_set1_ps<VTYPE>( 8.3321608736E-3f);
  const VTYPE s2 = vec_set1_ps<VTYPE>(-1.9515295891E-4f);
  const VTYPE c0 = vec_set1_ps<VTYPE>( 4.166664568298827E-2f);
  const VTYPE c1 = vec_set1_ps<VTYPE>(-1.388731625493765E-3f);
  const VTYPE c2 = vec_set1_ps<VTYPE>( 2.443315711809948E-5f);
  VTYPE xa, x, y, x2, s, c, sin1, cos1;
  ITYPE q, signsin, signcos;
  BVTYPE swap, overflow;
  xa = vec_abs(xx);
  q = vec_cvtps_epi32(xa * pi4);
  q = (q + 1) & ~1;
  y = -vec_cvtepi32_ps(q);
  x = vec_fmadd(y, d2, vec_fmadd(y, d1, vec_fmadd(y, d0, xa)));
  x2 = x * x;
  s = vec_fmadd(x2 * x2, s2, vec_fmadd(x2, s1, s0)) * (x * x2) + x;
  c = vec_fmadd(x2 * x2, c2, vec_fmadd(x2, c1, c0)) * (x2 * x2) + vec_fmadd(half, x2, one);
  swap = vec_neq(vec_and(q,2),0);
  overflow = vec_lt(q,0);
  if (horizontal_or(vec_and(overflow,is_finite(xa)))) {
    s = vec_select(overflow, zero, s);
    c = vec_select(overflow, one, c);
  }
  if (SC & 5) {
    sin1 = vec_select(swap, c, s);
    signsin = vec_and(vec_xor(vec_sll(q,29),xx),1 << 31);
    sin1 = vec_xor(sin1,signsin);
  }
  if (SC & 6) {
    cos1 = vec_select(swap, s, c);
    signcos = vec_and(vec_sll(q+2,29),1 << 31);
    cos1 = vec_xor(cos1,signcos);
  }
  if      (SC == 1) return sin1;
  else if (SC == 2) return cos1;
  else if (SC == 3) {
    *cosret = cos1;
    return sin1;
  }
}

static inline __m128 _mm_sin_ps(__m128 const & x) {
  return sincos_f<__m128, __m128i, __m128i, 1>(0, x);
}
static inline __m128 _mm_cos_ps(__m128 const & x) {
  return sincos_f<__m128, __m128i, __m128i, 2>(0, x);
}
static inline __m128 _mm_sincos_ps(__m128 * cosret, __m128 const & x) {
  return sincos_f<__m128, __m128i, __m128i, 3>(cosret, x);
}
static inline float _mm_reduce_add_ps(__m128 const & in) {
  union {
    __m128 temp;
    float out[4];
  };
  temp = _mm_hadd_ps(in,in);
  temp = _mm_hadd_ps(temp,temp);
  return out[0];
}

#ifdef __AVX__
static inline __m256 _mm256_sin_ps(__m256 const & x) {
  return sincos_f<__m256, __m256i, __m256i, 1>(0, x);
}
static inline __m256 _mm256_cos_ps(__m256 const & x) {
  return sincos_f<__m256, __m256i, __m256i, 2>(0, x);
}
static inline __m256 _mm256_sincos_ps(__m256 * cosret, __m256 const & x) {
  return sincos_f<__m256, __m256i, __m256i, 3>(cosret, x);
}
static inline float _mm256_reduce_add_ps(__m256 const &in) {
  union {
    __m256 temp;
    float out[8];
  };
  temp = _mm256_permute2f128_ps(in,in,1);
  temp = _mm256_add_ps(temp,in);
  temp = _mm256_hadd_ps(temp,temp);
  temp = _mm256_hadd_ps(temp,temp);
  return out[0];
}
#endif

#if __AVX512F__ | __MIC__
static inline __m512 _mm512_sin_ps(__m512 const & x) {
  return sincos_f<__m512, __m512i, __mmask16, 1>(0, x);
}
static inline __m512 _mm512_cos_ps(__m512 const & x) {
  return sincos_f<__m512, __m512i, __mmask16, 2>(0, x);
}
static inline __m512 _mm512_sincos_ps(__m512 * cosret, __m512 const & x) {
  return sincos_f<__m512, __m512i, __mmask16, 3>(cosret, x);
}
static inline float _mm512_reduce_add_ps(__m512 const & in) {
  __m256 temp = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(in),1));
  temp = _mm256_add_ps(temp,_mm512_castps512_ps256(in));
  return _mm256_reduce_add_ps(temp);
}
#endif

template<class VTYPE, class ITYPE, class ITYPEH, class BVTYPE, int SC>
static inline VTYPE sincos_d(VTYPE * cosret, VTYPE const & xx) {
  const VTYPE zero = vec_set1_pd<VTYPE>(0.);
  const VTYPE half = vec_set1_pd<VTYPE>(-.5);
  const VTYPE one = vec_set1_pd<VTYPE>(1.);
  const VTYPE pi4 = vec_set1_pd<VTYPE>(4./M_PI);
  const VTYPE d0 = vec_set1_pd<VTYPE>( 7.853981554508209228515625E-1);
  const VTYPE d1 = vec_set1_pd<VTYPE>( 7.94662735614792836714E-9);
  const VTYPE d2 = vec_set1_pd<VTYPE>( 3.06161699786838294307E-17);
  const VTYPE s0 = vec_set1_pd<VTYPE>(-1.66666666666666307295E-1);
  const VTYPE s1 = vec_set1_pd<VTYPE>( 8.33333333332211858878E-3);
  const VTYPE s2 = vec_set1_pd<VTYPE>(-1.98412698295895385996E-4);
  const VTYPE s3 = vec_set1_pd<VTYPE>( 2.75573136213857245213E-6);
  const VTYPE s4 = vec_set1_pd<VTYPE>(-2.50507477628578072866E-8);
  const VTYPE s5 = vec_set1_pd<VTYPE>( 1.58962301576546568060E-10);
  const VTYPE c0 = vec_set1_pd<VTYPE>( 4.16666666666665929218E-2);
  const VTYPE c1 = vec_set1_pd<VTYPE>(-1.38888888888730564116E-3);
  const VTYPE c2 = vec_set1_pd<VTYPE>( 2.48015872888517045348E-5);
  const VTYPE c3 = vec_set1_pd<VTYPE>(-2.75573141792967388112E-7);
  const VTYPE c4 = vec_set1_pd<VTYPE>( 2.08757008419747316778E-9);
  const VTYPE c5 = vec_set1_pd<VTYPE>(-1.13585365213876817300E-11);
  VTYPE xa, x, y, x2, x4, x8, s, c, sin1, cos1;
  ITYPE qq, signsin, signcos;
  ITYPEH q;
  BVTYPE swap, overflow;
  xa = vec_abs(xx);
  q = vec_cvtpd_epi32(xa * pi4);
  q = (q + 1) & ~1;
  y = -vec_cvtepi32_pd<VTYPE>(q);
  x = vec_fmadd(y, d2, vec_fmadd(y, d1, vec_fmadd(y, d0, xa)));
  x2 = x * x;
  x4 = x2 * x2;
  x8 = x4 * x4;
  s = vec_fmadd(vec_fmadd(s3,x2,s2), x4, vec_fmadd(vec_fmadd(s5,x2,s4), x8, vec_fmadd(s1,x2,s0)));
  c = vec_fmadd(vec_fmadd(c3,x2,c2), x4, vec_fmadd(vec_fmadd(c5,x2,c4), x8, vec_fmadd(c1,x2,c0)));
  s = vec_fmadd(x * x2, s, x);
  c = vec_fmadd(x2 * x2, c, vec_fmadd(x2, half, one));
  qq = vec_cvtepi32_epi64<ITYPE>(q);
  swap = vec_neq_64(vec_and_64(qq,2),int64_t(0));
  if (horizontal_or(vec_lt(q,0))) {
    overflow = vec_and(vec_lt(y,0),is_finite(xa));
    s = vec_select(overflow, zero, s);
    c = vec_select(overflow, one, c);
  }
  if (SC & 1) {
    sin1 = vec_select(swap, c, s);
    signsin = vec_and_64(vec_xor(vec_sll_64(qq,61),xx),1ULL << 63);
    sin1 = vec_xor(sin1,signsin);
  }
  if (SC & 2) {
    cos1 = vec_select(swap, s, c);
    signcos = vec_and_64(vec_sll_64(qq+2,61),1ULL << 63);
    cos1 = vec_xor(cos1,signcos);
  }
  if (SC == 3) {
    *cosret = cos1;
  }
  if (SC & 1) return sin1; else return cos1;
}

static inline __m128d _mm_sin_pd(__m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, __m128i, 1>(0, x);
}
static inline __m128d _mm_cos_pd(__m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, __m128i, 2>(0, x);
}
static inline __m128d _mm_sincos_pd(__m128d * cosret, __m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, __m128i, 3>(cosret, x);
}
static inline double _mm_reduce_add_pd(__m128d const & in) {
  union {
    __m128d temp;
    double out[2];
  };
  temp = _mm_hadd_pd(in,in);
  return out[0];
}

#ifdef __AVX__
static inline __m256d _mm256_sin_pd(__m256d const & x) {
  return sincos_d<__m256d, __m256i, __m128i, __m256i, 1>(0, x);
}
static inline __m256d _mm256_cos_pd(__m256d const & x) {
  return sincos_d<__m256d, __m256i, __m128i, __m256i, 2>(0, x);
}
static inline __m256d _mm256_sincos_pd(__m256d * cosret, __m256d const & x) {
  return sincos_d<__m256d, __m256i, __m128i, __m256i, 3>(cosret, x);
}
static inline double _mm256_reduce_add_pd(__m256d const & in) {
  union {
    __m256d temp;
    double out[4];
  };
  temp = _mm256_permute2f128_pd(in,in,1);
  temp = _mm256_add_pd(temp,in);
  temp = _mm256_hadd_pd(temp,temp);
  return out[0];
}
#endif

#if __AVX512F__ | __MIC__
static inline __m512d _mm512_sin_pd(__m512d const & x) {
  return sincos_d<__m512d, __m512i, __m256i, __mmask8, 1>(0, x);
}
static inline __m512d _mm512_cos_pd(__m512d const & x) {
  return sincos_d<__m512d, __m512i, __m256i, __mmask8, 2>(0, x);
}
static inline __m512d _mm512_sincos_pd(__m512d * cosret, __m512d const & x) {
  return sincos_d<__m512d, __m512i, __m256i, __mmask8, 3>(cosret, x);
}
static inline double _mm512_reduce_add_pd(__m512d const & in) {
  __m256d temp = _mm512_extractf64x4_pd(in,1);
  temp = _mm256_add_pd(temp,_mm512_castpd512_pd256(in));
  return _mm256_reduce_add_pd(temp);
}
#endif

template<class VTYPE, class BVTYPE>
static inline VTYPE exp_f(VTYPE const & initial_x) {
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
  VTYPE r, x, x2, x4, z, n2;
  BVTYPE inrange;
  x = initial_x;
  r = vec_round(initial_x * log2e);
  x = vec_fmadd(r, ln2f_hi, x);
  x = vec_fmadd(r, ln2f_lo, x);
  x2 = x * x;
  x4 = x2 * x2;
  z = vec_fmadd(vec_fmadd(p5,x,p4), x2, vec_fmadd(vec_fmadd(p7,x,p6), x4, vec_fmadd(p3,x,p2)));
  z = vec_fmadd(z, x2, x);
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
  return exp_f<__m128, __m128i>(x);
}
#ifdef __AVX__
static inline __m256 _mm256_exp_ps(__m256 const & x) {
  return exp_f<__m256, __m256i>(x);
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512 _mm512_exp_ps(__m512 const & x) {
  return exp_f<__m512, __mmask16>(x);
}
#endif

template<class VTYPE, class BVTYPE>
inline VTYPE exp_d(VTYPE const & initial_x) {
  const VTYPE zero = vec_set1_pd<VTYPE>(0.);
  const VTYPE p2 = vec_set1_pd<VTYPE>(1./2.);
  const VTYPE p3 = vec_set1_pd<VTYPE>(1./6.);
  const VTYPE p4 = vec_set1_pd<VTYPE>(1./24.);
  const VTYPE p5 = vec_set1_pd<VTYPE>(1./120.);
  const VTYPE p6 = vec_set1_pd<VTYPE>(1./720.);
  const VTYPE p7 = vec_set1_pd<VTYPE>(1./5040.);
  const VTYPE p8 = vec_set1_pd<VTYPE>(1./40320.);
  const VTYPE p9 = vec_set1_pd<VTYPE>(1./362880.);
  const VTYPE p10 = vec_set1_pd<VTYPE>(1./3628800.);
  const VTYPE p11 = vec_set1_pd<VTYPE>(1./39916800.);
  const VTYPE p12 = vec_set1_pd<VTYPE>(1./479001600.);
  const VTYPE p13 = vec_set1_pd<VTYPE>(1./6227020800.);
  const VTYPE log2e = vec_set1_pd<VTYPE>(1.44269504088896340736);
  const VTYPE ln2d_hi = vec_set1_pd<VTYPE>(-0.693145751953125);
  const VTYPE ln2d_lo = vec_set1_pd<VTYPE>(-1.42860682030941723212E-6);
  const double max_x = 708.39;
  VTYPE x, x2, x4, x8, r, z, n2;
  BVTYPE inrange;

  x = initial_x;
  r = vec_round(initial_x*log2e);
  x = vec_fmadd(r, ln2d_hi, x);
  x = vec_fmadd(r, ln2d_lo, x);
  x2 = x * x;
  x4 = x2 * x2;
  x8 = x4 * x4;
  z = vec_fmadd(vec_fmadd(vec_fmadd(p13,x,p12), x4, vec_fmadd(vec_fmadd(p11,x,p10), x2,
      vec_fmadd(p9,x,p8))), x8, vec_fmadd(vec_fmadd(vec_fmadd(p7,x,p6), x2,
      vec_fmadd(p5,x,p4)), x4, vec_fmadd(vec_fmadd(p3,x,p2),x2,x)));
  n2 = vec_pow2n(r);
  z = (z + 1.0) * n2;
  inrange = vec_lt(vec_abs(initial_x),max_x);
  inrange &= is_finite(initial_x);
  if (horizontal_and(inrange)) {
    return z;
  } else {
    r = vec_select(sign_bit(initial_x), zero, vec_inf<VTYPE>());
    z = vec_select(inrange, z, r);
    z = vec_select(is_nan(initial_x), initial_x, z);
    return z;
  }
}

static inline __m128d _mm_exp_pd(__m128d const & x) {
  return exp_d<__m128d, __m128i>(x);
}
#ifdef __AVX__
static inline __m256d _mm256_exp_pd(__m256d const & x) {
  return exp_d<__m256d, __m256i>(x);
}
#endif
#if __AVX512F__ | __MIC__
static inline __m512d _mm512_exp_pd(__m512d const & x) {
  return exp_d<__m512d, __mmask8>(x);
}
#endif

#endif
