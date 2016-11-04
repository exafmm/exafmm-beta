#ifndef VECTORMATH_TRIG_H
#define VECTORMATH_TRIG_H  1
#include "vectormath_common.h"

static inline __m128i vec_and(__m128i const & a, __m128i const & b) {
  return _mm_and_si128(a, b);
}
static inline __m128 vec_xor(__m128 const & a, __m128 const & b) {
  return _mm_xor_ps(a, b);
}
static inline Vec4fb vec_neq(__m128i const & a, __m128i const & b) {
  return Vec4fb(_mm_xor_si128(_mm_cmpeq_epi32(a, b), _mm_set1_epi32(-1)));
}
static inline Vec4fb vec_lt(__m128i const & a, __m128i const & b) {
  return Vec4fb(_mm_cmpgt_epi32(b, a));
}
template<class ITYPE>
static inline ITYPE vec_set_i32(int const & a);
template<>
inline __m128i vec_set_i32(int const & a) {
  return _mm_set1_epi32(a);
}
static inline Vec4i vm_truncate_low_to_int(Vec2d const & x) {
    return truncate_to_int(x,x);
}
template<class VTYPE, class ITYPE>
static inline VTYPE vm_half_int_vector_to_double(ITYPE const & x);

template<>
inline Vec2d vm_half_int_vector_to_double<Vec2d, Vec4i>(Vec4i const & x) {
    return to_double_low(x);
}
template<class ITYPE, class ITYPEH>
static inline ITYPE vm_half_int_vector_to_full(ITYPEH const & x);

template<>
inline Vec2q vm_half_int_vector_to_full<Vec2q,Vec4i>(Vec4i const & x) {
    return extend_low(x);
}

#if __AVX__
static inline __m256i vec_and(__m256i const & a, __m256i const & b) {
  return _mm256_and_si256(a, b);
}
static inline __m256 vec_xor(__m256 const & a, __m256 const & b) {
  return _mm256_xor_ps(a, b);
}
static inline Vec8fb vec_neq(__m256i const & a, __m256i const & b) {
  return Vec8fb(_mm256_xor_si256(_mm256_cmpeq_epi32(a, b), _mm256_set1_epi32(-1)));
}
static inline Vec8fb vec_lt(__m256i const & a, __m256i const & b) {
  return Vec8fb(_mm256_cmpgt_epi32(b, a));
}
template<>
inline __m256i vec_set_i32(int const & a) {
  return _mm256_set1_epi32(a);
}
static inline Vec4i vm_truncate_low_to_int(Vec4d const & x) {
    return truncate_to_int(x);
}
template<>
inline Vec4d vm_half_int_vector_to_double<Vec4d, Vec4i>(Vec4i const & x) {
    return to_double(x);
}
template<>
inline Vec4q vm_half_int_vector_to_full<Vec4q,Vec4i>(Vec4i const & x) {
    return extend_low(Vec8i(x,x));
}
#endif

#if __AVX512F__ | __MIC__
static inline __m512i vec_and(__m512i const & a, __m512i const & b) {
  return _mm512_and_epi32(a, b);
}
static inline __m512 vec_xor(__m512 const & a, __m512 const & b) {
  return _mm512_castsi512_ps(Vec16i(_mm512_castps_si512(a)) ^ Vec16i(_mm512_castps_si512(b)));
}
static inline Vec16fb vec_neq(__m512i const & a, __m512i const & b) {
  return Vec16fb(_mm512_cmpneq_epi32_mask(a,b));
}
static inline Vec16fb vec_lt(__m512i const & a, __m512i const & b) {
  return Vec16fb(_mm512_cmpgt_epi32_mask(b, a));
}
template<>
inline __m512i vec_set_i32(int const & a) {
  return _mm512_set1_epi32(a);
}
static inline Vec8i vm_truncate_low_to_int(Vec8d const & x) {
    return truncate_to_int(x);
}
template<>
inline Vec8d vm_half_int_vector_to_double<Vec8d, Vec8i>(Vec8i const & x) {
    return to_double(x);
}
template<>
inline Vec8q vm_half_int_vector_to_full<Vec8q,Vec8i>(Vec8i const & x) {
    return extend_low(Vec16i(x,x));
}
#endif

template<class VTYPE, class ITYPE, class ITYPE2, class BVTYPE, int SC>
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

    VTYPE  xa, x, y, x2, s, c, sin1, cos1;  // data vectors
    ITYPE  q, signsin, signcos;             // integer vectors
    ITYPE2 qq;
    BVTYPE swap, overflow;                  // boolean vectors

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
    swap = vec_neq(vec_and(q,vec_set_i32<ITYPE>(2)),vec_set_i32<ITYPE>(0));
    qq = q;

    // check for overflow
    overflow = vec_lt(qq,vec_set_i32<ITYPE>(0));
    if (horizontal_or(overflow & is_finite(xa))) {
        s = select(overflow, 0.f, s);
        c = select(overflow, 1.f, c);
    }

    if (SC & 5) {  // calculate sin
        sin1 = select(swap, c, s);
        signsin = ((qq << 29) ^ reinterpret_i(xx)) & vec_set_i32<ITYPE>(1 << 31);
        sin1 = vec_xor(sin1,reinterpret_f(signsin));
    }
    if (SC & 6) {  // calculate cos
        cos1 = select(swap, s, c);
        signcos = ((qq + 2) << 29) & (1 << 31);
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
  return sincos_f<__m128, __m128i, Vec4i, Vec4fb, 1>(0, x);
}
static inline __m128 _mm_cos_ps(__m128 const & x) {
  return sincos_f<__m128, __m128i, Vec4i, Vec4fb, 2>(0, x);
}
static inline __m128 _mm_sincos_ps(__m128 * cosret, __m128 const & x) {
  return sincos_f<__m128, __m128i, Vec4i, Vec4fb, 3>(cosret, x);
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
  return sincos_f<__m256, __m256i, Vec8i, Vec8fb, 1>(0, x);
}
static inline __m256 _mm256_cos_ps(__m256 const & x) {
  return sincos_f<__m256, __m256i, Vec8i, Vec8fb, 2>(0, x);
}
static inline __m256 _mm256_sincos_ps(__m256 * cosret, __m256 const & x) {
  return sincos_f<__m256, __m256i, Vec8i, Vec8fb, 3>(cosret, x);
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
  return sincos_f<__m512, __m512i, Vec16i, Vec16fb, 1>(0, x);
}
static inline __m512 _mm512_cos_ps(__m512 const & x) {
  return sincos_f<__m512, __m512i, Vec16i, Vec16fb, 2>(0, x);
}
static inline __m512 _mm512_sincos_ps(__m512 * cosret, __m512 const & x) {
  return sincos_f<__m512, __m512i, Vec16i, Vec16fb, 3>(cosret, x);
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
    /*
    const double DP1sc = 7.85398125648498535156E-1;
    const double DP2sc = 3.77489470793079817668E-8;
    const double DP3sc = 2.69515142907905952645E-15;
    */
    VTYPE  xa, x, y, x2, s, c, sin1, cos1;       // data vectors
    ITYPEH q;                                    // integer vectors, 32 bit
    ITYPE  qq, signsin, signcos;                 // integer vectors, 64 bit
    BVTYPE swap, overflow;                       // boolean vectors

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
    swap = BVTYPE((qq & 2) != 0);

    // check for overflow
    if (horizontal_or(q < 0)) {
        overflow = (y < 0) & is_finite(xa);
        s = select(overflow, 0., s);
        c = select(overflow, 1., c);
    }
    if (SC & 1) {  // calculate sin
        sin1 = select(swap, c, s);
        signsin = ((qq << 61) ^ ITYPE(reinterpret_i(xx))) & ITYPE(1ULL << 63);
        sin1 ^= reinterpret_d(signsin);
    }
    if (SC & 2) {  // calculate cos
        cos1 = select(swap, s, c);
        signcos = ((qq + 2) << 61) & (1ULL << 63);
        cos1 ^= reinterpret_d(signcos);
    }
    if (SC == 3) {  // calculate both. cos returned through pointer
        *cosret = cos1;
    }
    if (SC & 1) return sin1; else return cos1;
}

// instantiations of sincos_d template:

static inline __m128d _mm_sin_pd(__m128d const & x) {
  return sincos_d<Vec2d, Vec2q, Vec4i, Vec2db, 1>(0, x);
}
static inline __m128d _mm_cos_pd(__m128d const & x) {
  return sincos_d<Vec2d, Vec2q, Vec4i, Vec2db, 2>(0, x);
}
static inline __m128d _mm_sincos_pd(__m128d * cosret, __m128d const & x) {
  return sincos_d<Vec2d, Vec2q, Vec4i, Vec2db, 3>((Vec2d*)cosret, x);
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
  return sincos_d<Vec4d, Vec4q, Vec4i, Vec4db, 1>(0, x);
}
static inline __m256d _mm256_cos_pd(__m256d const & x) {
  return sincos_d<Vec4d, Vec4q, Vec4i, Vec4db, 2>(0, x);
}
static inline __m256d _mm256_sincos_pd(__m256d * cosret, __m256d const & x) {
  return sincos_d<Vec4d, Vec4q, Vec4i, Vec4db, 3>((Vec4d*)cosret, x);
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
  return sincos_d<Vec8d, Vec8q, Vec8i, Vec8db, 1>(0, x);
}
static inline __m512d _mm512_cos_pd(__m512d const & x) {
  return sincos_d<Vec8d, Vec8q, Vec8i, Vec8db, 2>(0, x);
}
static inline __m512d _mm512_sincos_pd(__m512d * cosret, __m512d const & x) {
  return sincos_d<Vec8d, Vec8q, Vec8i, Vec8db, 3>((Vec8d*)cosret, x);
}
static inline double _mm512_reduce_add_pd(__m512d const & in) {
  __m256d temp = _mm512_extractf64x4_pd(in,1);
  temp = _mm256_add_pd(temp,_mm512_castpd512_pd256(in));
  return _mm256_reduce_add_pd(temp);
}
#endif

#endif
