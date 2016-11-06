#ifndef VECTORMATH_TRIG_H
#define VECTORMATH_TRIG_H  1
#include "vectormath_common.h"

static inline __m128i vec_and(__m128i const & a, int const & b) {
  return _mm_and_si128(a,_mm_set1_epi32(b));
}
static inline __m128i vec_and(__m128i const & a, __m128i const & b) {
  return _mm_and_si128(a,b);
}
static inline __m128i vec_and_64(__m128i const & a, int64_t const & b) {
  return _mm_and_si128(a,_mm_set1_epi64x(b));
}
static inline __m128 vec_xor(__m128 const & a, __m128 const & b) {
  return _mm_xor_ps(a,b);
}
static inline __m128i vec_xor(__m128i const & a, __m128i const & b) {
  return _mm_xor_si128(a,b);
}
static inline __m128i vec_xor_64(__m128i const & a, __m128d const & b) {
  return _mm_xor_si128(a,_mm_castpd_si128(b));
}
static inline __m128d vec_xor(__m128d const & a, __m128d const & b) {
  return _mm_xor_pd(a,b);
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
static inline __m128i vec_lt(__m128d const & a, double const & b) {
  return _mm_castpd_si128(_mm_cmplt_pd(_mm_set1_pd(b),a));
}
static inline __m128i vec_sll(__m128i const & a, int const & b) {
  return _mm_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m128i vec_sll_64(__m128i const & a, int const & b) {
  return _mm_sll_epi64(a,_mm_cvtsi32_si128(b));
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
static inline __m128i truncate_to_int(__m128 const & a) {
    return _mm_cvttps_epi32(a);
}
static inline __m128 select(__m128i const & s, __m128 const & a, __m128 const & b) {
  return _mm_blendv_ps(b,a,_mm_castsi128_ps(s));
}
static inline __m128d selectd(__m128i const & s, __m128d const & a, __m128d const & b) {
  return _mm_blendv_pd(b,a,_mm_castsi128_pd(s));
}
static inline __m128i is_finite(__m128 const & a) {
  return vec_neq(vec_and(vec_sll(_mm_castps_si128(a),1),0xFF000000),0xFF000000);
}
static inline __m128i is_finite(__m128d const & a) {
  return vec_neq_64(vec_and_64(vec_sll_64(_mm_castpd_si128(a),1),0xFFE0000000000000ll),0xFFE0000000000000ll);
}
static inline bool horizontal_or (__m128i const & a) {
  return ! _mm_testz_si128(a,a);
}

static inline __m128i vm_truncate_low_to_int(__m128d const & x) {
  static const union {
    int     i[4];
    __m128i xmm;
  } u = {{-1,-1,0,0}};
  __m128i a1 = _mm_and_si128(_mm_cvttpd_epi32(x), u.xmm);
  __m128i b1 = _mm_slli_si128(_mm_cvttpd_epi32(x), 8);
  return _mm_or_si128(a1,b1);
}
template<class VTYPE, class ITYPE>
static inline VTYPE vm_half_int_vector_to_double(ITYPE const & x);
template<>
inline __m128d vm_half_int_vector_to_double<__m128d,__m128i>(__m128i const & x) {
  return _mm_cvtepi32_pd(x);
}
template<class ITYPE, class ITYPEH>
static inline ITYPE vm_half_int_vector_to_full(ITYPEH const & x);
template<>
inline __m128i vm_half_int_vector_to_full<__m128i,__m128i>(__m128i const & x) {
  __m128i sign = _mm_srai_epi32(x,31);
  return _mm_unpacklo_epi32(x,sign);
}

#if __AVX__
static inline __m256i vec_and(__m256i const & a, int const & b) {
  return _mm256_and_si256(a,_mm256_set1_epi32(b));
}
static inline __m256i vec_and(__m256i const & a, __m256i const & b) {
  return _mm256_and_si256(a,b);
}
static inline __m256i vec_and_64(__m256i const & a, int64_t const & b) {
  return _mm256_and_si256(a,_mm256_set1_epi64x(b));
}
static inline __m256 vec_xor(__m256 const & a, __m256 const & b) {
  return _mm256_xor_ps(a,b);
}
static inline __m256i vec_xor(__m256i const & a, __m256i const & b) {
  return _mm256_xor_si256(a,b);
}
static inline __m256i vec_xor_64(__m256i const & a, __m256d const & b) {
  return _mm256_xor_si256(a,_mm256_castpd_si256(b));
}
static inline __m256d vec_xor(__m256d const & a, __m256d const & b) {
  return _mm256_xor_pd(a,b);
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
static inline __m256i vec_lt(__m256d const & a, double const & b) {
  return _mm256_castpd_si256(_mm256_cmp_pd(_mm256_set1_pd(b),a,1));
}
static inline __m256i vec_sll(__m256i const & a, int const & b) {
  return _mm256_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m256i vec_sll_64(__m256i const & a, int const & b) {
  return _mm256_sll_epi64(a,_mm_cvtsi32_si128(b));
}
template<>
inline __m256 vec_set1_ps(float const & a) {
  return _mm256_set1_ps(a);
}
template<>
inline __m256d vec_set1_pd(double const & a) {
  return _mm256_set1_pd(a);
}
static inline __m256i truncate_to_int(__m256 const & a) {
  __m128i lo = _mm_cvttps_epi32(_mm256_castps256_ps128(a));
  __m128i hi = _mm_cvttps_epi32(_mm256_extractf128_ps(a,1));
  return _mm256_inserti128_si256(_mm256_castsi128_si256(lo),(hi),1);
}
static inline __m256 select(__m256i const & s, __m256 const & a, __m256 const & b) {
  return _mm256_blendv_ps(b,a,_mm256_castsi256_ps(s));
}
static inline __m256d selectd(__m256i const & s, __m256d const & a, __m256d const & b) {
  return _mm256_blendv_pd(b,a,_mm256_castsi256_pd(s));
}
static inline __m256i is_finite(__m256 const & a) {
  return vec_neq(vec_and(vec_sll(_mm256_castps_si256(a),1),0xFF000000),0xFF000000);
}
static inline __m256i is_finite(__m256d const & a) {
  return vec_neq_64(vec_and_64(vec_sll_64(_mm256_castpd_si256(a),1),0xFFE0000000000000ll),0xFFE0000000000000ll);
}
static inline bool horizontal_or (__m256i const & a) {
  return ! _mm256_testz_si256(a,a);
}

static inline __m128i vm_truncate_low_to_int(__m256d const & x) {
  return _mm256_cvttpd_epi32(x);
}
template<>
inline __m256d vm_half_int_vector_to_double<__m256d,__m128i>(__m128i const & x) {
  return _mm256_cvtepi32_pd(x);
}
template<>
inline __m256i vm_half_int_vector_to_full<__m256i,__m128i>(__m128i const & x) {
  __m256i a = _mm256_inserti128_si256(_mm256_castsi128_si256(x),x,1);
  __m256i a2 = permute4q<0,-256,1,-256>(Vec4q(a));
  __m256i sign = _mm256_srai_epi32(a2,31);
  return _mm256_unpacklo_epi32(a2,sign);
}
#endif

#if __AVX512F__ | __MIC__
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
static inline __m512 vec_xor(__m512 const & a, __m512 const & b) {
  return _mm512_castsi512_ps(Vec16i(_mm512_castps_si512(a))^Vec16i(_mm512_castps_si512(b)));
}
static inline __m512i vec_xor(__m512i const & a, __m512i const & b) {
  return _mm512_xor_epi32(a,b);
}
static inline __m512i vec_xor_64(__m512i const & a, __m512d const & b) {
  return _mm512_xor_epi64(a,_mm512_castpd_si512(b));
}
static inline __m512d vec_xor(__m512d const & a, __m512d const & b) {
  return _mm512_castsi512_pd(Vec8q(_mm512_castpd_si512(a))^Vec8q(_mm512_castpd_si512(b)));
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
static inline __mmask8 vec_lt(__m512d const & a, double const & b) {
  return _mm512_cmp_pd_mask(a,_mm512_set1_pd(b),1);
}
static inline __m512i vec_sll(__m512i const & a, int const & b) {
  return _mm512_sll_epi32(a,_mm_cvtsi32_si128(b));
}
static inline __m512i vec_sll_64(__m512i const & a, int const & b) {
  return _mm512_sll_epi64(a,_mm_cvtsi32_si128(b));
}
template<>
inline __m512 vec_set1_ps(float const & a) {
  return _mm512_set1_ps(a);
}
template<>
inline __m512d vec_set1_pd(double const & a) {
  return _mm512_set1_pd(a);
}
static inline __m512i truncate_to_int(__m512 const & a) {
    return _mm512_cvtt_roundps_epi32(a, _MM_FROUND_NO_EXC);
}
static inline __m512 select(__mmask16 const & s, __m512 const & a, __m512 const & b) {
  return _mm512_mask_mov_ps(b,s,a);
}
static inline __m512d selectd(__mmask8 const & s, __m512d const & a, __m512d const & b) {
  return _mm512_mask_mov_pd(b,s,a);
}
static inline __mmask16 is_finite(__m512 const & a) {
  return vec_neq(vec_and(vec_sll(_mm512_castps_si512(a),1),0xFF000000),0xFF000000);
}
static inline __mmask8 is_finite(__m512d const & a) {
  return vec_neq_64(vec_and_64(vec_sll_64(_mm512_castpd_si512(a),1),0xFFE0000000000000ll),0xFFE0000000000000ll);
}
static inline bool horizontal_or (__mmask16 const & a) {
  return (uint16_t)(a != 0);
}

static inline __m256i vm_truncate_low_to_int(__m512d const & x) {
  return _mm512_cvtpd_epi32(x);
}
template<>
inline __m512d vm_half_int_vector_to_double<__m512d,__m256i>(__m256i const & x) {
  return _mm512_cvtepi32_pd(x);
}
template<>
inline __m512i vm_half_int_vector_to_full<__m512i,__m256i>(__m256i const & x) {
  return _mm512_cvtepi32_epi64(_mm512_castsi512_si256(_mm512_inserti64x4(_mm512_castsi256_si512(x),x,1)));
}
#endif

template<class VTYPE, class ITYPE, class BVTYPE, int SC>
static inline VTYPE sincos_f(VTYPE * cosret, VTYPE const & xx) {
    const float ONEOPIO4f = (float)(4./VM_PI);

    const float DP1F = 0.78515625f;
    const float DP2F = 2.4187564849853515625E-4f;
    const float DP3F = 3.77489497744594108E-8f;

    const float P0sinf = -1.6666654611E-1f;
    const float P1sinf =  8.3321608736E-3f;
    const float P2sinf = -1.9515295891E-4f;

    const float P0cosf =  4.166664568298827E-2f;
    const float P1cosf = -1.388731625493765E-3f;
    const float P2cosf =  2.443315711809948E-5f;

    VTYPE  xa, x, y, x2, s, c, sin1, cos1;
    ITYPE  q, signsin, signcos;
    BVTYPE swap, overflow;

    xa = abs(xx);

    // Find quadrant
    //      0 -   pi/4 => 0
    //   pi/4 - 3*pi/4 => 2
    // 3*pi/4 - 5*pi/4 => 4
    // 5*pi/4 - 7*pi/4 => 6
    // 7*pi/4 - 8*pi/4 => 8
    q = truncate_to_int(xa * ONEOPIO4f);
    q = (q + 1) & ~1;

    y = to_float(q);  // quadrant, as float

    // Reduce by extended precision modular arithmetic
    x = nmul_add(y, DP3F, nmul_add(y, DP2F, nmul_add(y, DP1F, xa))); // x = ((xa - y * DP1F) - y * DP2F) - y * DP3F;

    // A two-step reduction saves time at the cost of precision for very big x:
    //x = (xa - y * DP1F) - y * (DP2F+DP3F);

    // Taylor expansion of sin and cos, valid for -pi/4 <= x <= pi/4
    x2 = x * x;
    s = polynomial_2(x2, P0sinf, P1sinf, P2sinf) * (x*x2)  + x;
    c = polynomial_2(x2, P0cosf, P1cosf, P2cosf) * (x2*x2) + nmul_add(0.5f, x2, 1.0f);

    // correct for quadrant
    swap = vec_neq(vec_and(q,2),0);

    // check for overflow
    overflow = vec_lt(q,0);
    if (horizontal_or(vec_and(overflow,is_finite(xa)))) {
      s = select(overflow, vec_set1_ps<VTYPE>(0.f), s);
      c = select(overflow, vec_set1_ps<VTYPE>(1.f), c);
    }

    if (SC & 5) {  // calculate sin
      sin1 = select(swap, c, s);
      signsin = vec_and(vec_xor(vec_sll(q,29),reinterpret_i(xx)),1 << 31);
      sin1 = vec_xor(sin1,reinterpret_f(signsin));
    }
    if (SC & 6) {  // calculate cos
      cos1 = select(swap, s, c);
      signcos = vec_and(vec_sll(q+2,29),1 << 31);
      cos1 = vec_xor(cos1,reinterpret_f(signcos));
    }
    if      (SC == 1) return sin1;
    else if (SC == 2) return cos1;
    else if (SC == 3) {  // calculate both. cos returned through pointer
      *cosret = cos1;
      return sin1;
    }
}

// instantiations of sincos_f template:

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
    // define constants
  const double ONEOPIO4 = 4./VM_PI;

  const double P0sin =-1.66666666666666307295E-1;
  const double P1sin = 8.33333333332211858878E-3;
  const double P2sin =-1.98412698295895385996E-4;
  const double P3sin = 2.75573136213857245213E-6;
  const double P4sin =-2.50507477628578072866E-8;
  const double P5sin = 1.58962301576546568060E-10;

  const double P0cos = 4.16666666666665929218E-2;
  const double P1cos =-1.38888888888730564116E-3;
  const double P2cos = 2.48015872888517045348E-5;
  const double P3cos =-2.75573141792967388112E-7;
  const double P4cos = 2.08757008419747316778E-9;
  const double P5cos =-1.13585365213876817300E-11;

  const double DP1 = 7.853981554508209228515625E-1;
  const double DP2 = 7.94662735614792836714E-9;
  const double DP3 = 3.06161699786838294307E-17;

  VTYPE  xa, x, y, x2, s, c, sin1, cos1;
  ITYPE qq, signsin, signcos;
  ITYPEH q;
  BVTYPE swap, overflow;

  xa = abs(xx);

  // Find quadrant
  //      0 -   pi/4 => 0
  //   pi/4 - 3*pi/4 => 2
  // 3*pi/4 - 5*pi/4 => 4
  // 5*pi/4 - 7*pi/4 => 6
  // 7*pi/4 - 8*pi/4 => 8

  // truncate to integer (magic number conversion is not faster here)
  q = vm_truncate_low_to_int(xa * ONEOPIO4);
  q = (q + 1) & ~1;

  y = vm_half_int_vector_to_double<VTYPE>(q);  // quadrant, as double

  // Reduce by extended precision modular arithmetic
  x = nmul_add(y, DP3, nmul_add(y, DP2, nmul_add(y, DP1, xa)));    // x = ((xa - y * DP1) - y * DP2) - y * DP3;

  // Expansion of sin and cos, valid for -pi/4 <= x <= pi/4
  x2 = x * x;
  s = polynomial_5(x2, P0sin, P1sin, P2sin, P3sin, P4sin, P5sin);
  c = polynomial_5(x2, P0cos, P1cos, P2cos, P3cos, P4cos, P5cos);
  s = mul_add(x * x2, s, x);                                       // s = x + (x * x2) * s;
  c = mul_add(x2 * x2, c, nmul_add(x2, 0.5, 1.0));                 // c = 1.0 - x2 * 0.5 + (x2 * x2) * c;

  // correct for quadrant
  qq = vm_half_int_vector_to_full<ITYPE,ITYPEH>(q);
  swap = vec_neq_64(vec_and_64(qq,2),int64_t(0));

  // check for overflow
  if (horizontal_or(q < 0)) {
    overflow = vec_and(vec_lt(y,0),is_finite(xa));
    s = selectd(overflow,vec_set1_pd<VTYPE>(0.),s);
    c = selectd(overflow,vec_set1_pd<VTYPE>(1.),c);
  }
  if (SC & 1) {  // calculate sin
    sin1 = selectd(swap, c, s);
    signsin = vec_and_64(vec_xor_64(vec_sll_64(qq,61),xx),1ULL << 63);
    sin1 = vec_xor(sin1,reinterpret_d(signsin));
  }
  if (SC & 2) {  // calculate cos
    cos1 = selectd(swap, s, c);
    signcos = vec_and_64(vec_sll_64(qq+2,61),1ULL << 63);
    cos1 = vec_xor(cos1,reinterpret_d(signcos));
  }
  if (SC == 3) {  // calculate both. cos returned through pointer
    *cosret = cos1;
  }
  if (SC & 1) return sin1; else return cos1;
}

// instantiations of sincos_d template:

static inline __m128d _mm_sin_pd(__m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, Vec2db, 1>(0, x);
}
static inline __m128d _mm_cos_pd(__m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, Vec2db, 2>(0, x);
}
static inline __m128d _mm_sincos_pd(__m128d * cosret, __m128d const & x) {
  return sincos_d<__m128d, __m128i, __m128i, Vec2db, 3>(cosret, x);
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
  return sincos_d<__m256d, __m256i, __m128i, Vec4db, 1>(0, x);
}
static inline __m256d _mm256_cos_pd(__m256d const & x) {
  return sincos_d<__m256d, __m256i, __m128i, Vec4db, 2>(0, x);
}
static inline __m256d _mm256_sincos_pd(__m256d * cosret, __m256d const & x) {
  return sincos_d<__m256d, __m256i, __m128i, Vec4db, 3>(cosret, x);
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
  return sincos_d<__m512d, __m512i, __m256i, Vec8db, 1>(0, x);
}
static inline __m512d _mm512_cos_pd(__m512d const & x) {
  return sincos_d<__m512d, __m512i, __m256i, Vec8db, 2>(0, x);
}
static inline __m512d _mm512_sincos_pd(__m512d * cosret, __m512d const & x) {
  return sincos_d<__m512d, __m512i, __m256i, Vec8db, 3>(cosret, x);
}
static inline double _mm512_reduce_add_pd(__m512d const & in) {
  __m256d temp = _mm512_extractf64x4_pd(in,1);
  temp = _mm256_add_pd(temp,_mm512_castpd512_pd256(in));
  return _mm256_reduce_add_pd(temp);
}
#endif

#endif
