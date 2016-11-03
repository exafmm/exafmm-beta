/****************************  vectormath_exp.h   ******************************
* Author:        Agner Fog
* Date created:  2014-04-18
* Last modified: 2015-02-10
* Version:       1.16
* Project:       vector classes
* Description:
* Header file containing inline vector functions of logarithms, exponential
* and power functions:
* exp         exponential function
* exmp1       exponential function minus 1
* log         natural logarithm
* log1p       natural logarithm of 1+x
* cbrt        cube root
* pow         raise vector elements to power
* pow_ratio   raise vector elements to rational power
*
* Theory, methods and inspiration based partially on these sources:
* > Moshier, Stephen Lloyd Baluk: Methods and programs for mathematical functions.
*   Ellis Horwood, 1989.
* > VDT library developed on CERN by Danilo Piparo, Thomas Hauth and
*   Vincenzo Innocente, 2012, https://svnweb.cern.ch/trac/vdt
* > Cephes math library by Stephen L. Moshier 1992,
*   http://www.netlib.org/cephes/
*
* For detailed instructions, see vectormath_common.h and VectorClass.pdf
*
* (c) Copyright 2014 GNU General Public License http://www.gnu.org/licenses
******************************************************************************/

#ifndef VECTORMATH_EXP_H
#define VECTORMATH_EXP_H  1

#include "vectormath_common.h"


/******************************************************************************
*                 Exponential functions
******************************************************************************/

// Helper functions, used internally:

// This function calculates pow(2,n) where n must be an integer. Does not check for overflow or underflow
static inline Vec2d vm_pow2n (Vec2d const & n) {
    const double pow2_52 = 4503599627370496.0;   // 2^52
    const double bias = 1023.0;                  // bias in exponent
    Vec2d a = n + (bias + pow2_52);              // put n + bias in least significant bits
    Vec2q b = reinterpret_i(a);                  // bit-cast to integer
    Vec2q c = b << 52;                           // shift left 52 places to get into exponent field
    Vec2d d = reinterpret_d(c);                  // bit-cast back to double
    return d;
}

static inline Vec4f vm_pow2n (Vec4f const & n) {
    const float pow2_23 =  8388608.0;            // 2^23
    const float bias = 127.0;                    // bias in exponent
    Vec4f a = n + (bias + pow2_23);              // put n + bias in least significant bits
    Vec4i b = reinterpret_i(a);                  // bit-cast to integer
    Vec4i c = b << 23;                           // shift left 23 places to get into exponent field
    Vec4f d = reinterpret_f(c);                  // bit-cast back to float
    return d;
}

#if MAX_VECTOR_SIZE >= 256

static inline Vec4d vm_pow2n (Vec4d const & n) {
    const double pow2_52 = 4503599627370496.0;   // 2^52
    const double bias = 1023.0;                  // bias in exponent
    Vec4d a = n + (bias + pow2_52);              // put n + bias in least significant bits
    Vec4q b = reinterpret_i(a);                  // bit-cast to integer
    Vec4q c = b << 52;                           // shift left 52 places to get value into exponent field
    Vec4d d = reinterpret_d(c);                  // bit-cast back to double
    return d;
}

static inline Vec8f vm_pow2n (Vec8f const & n) {
    const float pow2_23 =  8388608.0;            // 2^23
    const float bias = 127.0;                    // bias in exponent
    Vec8f a = n + (bias + pow2_23);              // put n + bias in least significant bits
    Vec8i b = reinterpret_i(a);                  // bit-cast to integer
    Vec8i c = b << 23;                           // shift left 23 places to get into exponent field
    Vec8f d = reinterpret_f(c);                  // bit-cast back to float
    return d;
}

#endif // MAX_VECTOR_SIZE >= 256

#if MAX_VECTOR_SIZE >= 512

static inline Vec8d vm_pow2n (Vec8d const & n) {
    const double pow2_52 = 4503599627370496.0;   // 2^52
    const double bias = 1023.0;                  // bias in exponent
    Vec8d a = n + (bias + pow2_52);              // put n + bias in least significant bits
    Vec8q b = Vec8q(reinterpret_i(a));           // bit-cast to integer
    Vec8q c = b << 52;                           // shift left 52 places to get value into exponent field
    Vec8d d = Vec8d(reinterpret_d(c));           // bit-cast back to double
    return d;
}

static inline Vec16f vm_pow2n (Vec16f const & n) {
    const float pow2_23 =  8388608.0;            // 2^23
    const float bias = 127.0;                    // bias in exponent
    Vec16f a = n + (bias + pow2_23);             // put n + bias in least significant bits
    Vec16i b = Vec16i(reinterpret_i(a));         // bit-cast to integer
    Vec16i c = b << 23;                          // shift left 23 places to get into exponent field
    Vec16f d = Vec16f(reinterpret_f(c));         // bit-cast back to float
    return d;
}

#endif // MAX_VECTOR_SIZE >= 512


// Template for exp function, double precision
// The limit of abs(x) is defined by max_x below
// This function does not produce denormals
// Template parameters:
// VTYPE:  double vector type
// BVTYPE: boolean vector type
// M1: 0 for exp, 1 for expm1
// BA: 0 for exp, 1 for 0.5*exp, 2 for pow(2,x), 10 for pow(10,x)

#if 1  // choose method

// Taylor expansion
template<class VTYPE, class BVTYPE, int M1, int BA>
static inline VTYPE exp_d(VTYPE const & initial_x) {

    // Taylor coefficients, 1/n!
    // Not using minimax approximation because we prioritize precision close to x = 0
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

    // maximum abs(x), value depends on BA, defined below
    // The lower limit of x is slightly more restrictive than the upper limit.
    // We are specifying the lower limit, except for BA = 1 because it is not used for negative x
    double max_x;

    // data vectors
    VTYPE  x, r, z, n2;
    BVTYPE inrange;                              // boolean vector

    if (BA <= 1) { // exp(x)
        max_x = BA == 0 ? 708.39 : 709.7; // lower limit for 0.5*exp(x) is -707.6, but we are using 0.5*exp(x) only for positive x in hyperbolic functions
        const double ln2d_hi = 0.693145751953125;
        const double ln2d_lo = 1.42860682030941723212E-6;
        x  = initial_x;
        r  = round(initial_x*VM_LOG2E);
        // subtraction in two steps for higher precision
        x = nmul_add(r, ln2d_hi, x);             //  x -= r * ln2d_hi;
        x = nmul_add(r, ln2d_lo, x);             //  x -= r * ln2d_lo;
    }
    else if (BA == 2) { // pow(2,x)
        max_x = 1022.0;
        r  = round(initial_x);
        x  = initial_x - r;
        x *= VM_LN2;
    }
    else if (BA == 10) { // pow(10,x)
        max_x = 307.65;
        const double log10_2_hi = 0.30102999554947019;     // log10(2) in two parts
        const double log10_2_lo = 1.1451100899212592E-10;
        x  = initial_x;
        r  = round(initial_x*(VM_LOG2E*VM_LN10));
        // subtraction in two steps for higher precision
        x  = nmul_add(r, log10_2_hi, x);         //  x -= r * log10_2_hi;
        x  = nmul_add(r, log10_2_lo, x);         //  x -= r * log10_2_lo;
        x *= VM_LN10;
    }
    else  {  // undefined value of BA
        return 0.;
    }

    z = polynomial_13m(x, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13);

    if (BA == 1) r--;  // 0.5 * exp(x)

    // multiply by power of 2
    n2 = vm_pow2n(r);

    if (M1 == 0) {
        // exp
        z = (z + 1.0) * n2;
    }
    else {
        // expm1
        z = mul_add(z, n2, n2 - 1.0);            // z = z * n2 + (n2 - 1.0);
    }

    // check for overflow
    inrange  = abs(initial_x) < max_x;
    // check for INF and NAN
    inrange &= is_finite(initial_x);

    if (horizontal_and(inrange)) {
        // fast normal path
        return z;
    }
    else {
        // overflow, underflow and NAN
        r = select(sign_bit(initial_x), 0.-M1, infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        z = select(inrange, z, r);                                     // +/- underflow
        z = select(is_nan(initial_x), initial_x, z);                   // NAN goes through
        return z;
    }
}

#else

// Pade expansion uses less code and fewer registers, but is slower
template<class VTYPE, class BVTYPE, int M1, int BA>
static inline VTYPE exp_d(VTYPE const & initial_x) {

    // define constants
    const double ln2p1   = 0.693145751953125;
    const double ln2p2   = 1.42860682030941723212E-6;
    const double log2e   = VM_LOG2E;
    const double max_exp = 708.39;
    // coefficients of pade polynomials
    const double P0exp = 9.99999999999999999910E-1;
    const double P1exp = 3.02994407707441961300E-2;
    const double P2exp = 1.26177193074810590878E-4;
    const double Q0exp = 2.00000000000000000009E0;
    const double Q1exp = 2.27265548208155028766E-1;
    const double Q2exp = 2.52448340349684104192E-3;
    const double Q3exp = 3.00198505138664455042E-6;

    VTYPE  x, r, xx, px, qx, y, n2;              // data vectors
    BVTYPE inrange;                              // boolean vector

    x = initial_x;
    r = round(initial_x*log2e);

    // subtraction in one step would gives loss of precision
    x -= r * ln2p1;
    x -= r * ln2p2;

    xx = x * x;

    // px = x * P(x^2).
    px = polynomial_2(xx, P0exp, P1exp, P2exp) * x;

    // Evaluate Q(x^2).
    qx = polynomial_3(xx, Q0exp, Q1exp, Q2exp, Q3exp);

    // e^x = 1 + 2*P(x^2)/( Q(x^2) - P(x^2) )
    y = (2.0 * px) / (qx - px);

    // Get 2^n in double.
    // n  = round_to_int64_limited(r);
    // n2 = exp2(n);
    n2 = vm_pow2n(r);  // this is faster

    if (M1 == 0) {
        // exp
        y = (y + 1.0) * n2;
    }
    else {
        // expm1
        y = y * n2 + (n2 - 1.0);
    }

    // overflow
    inrange  = abs(initial_x) < max_exp;
    // check for INF and NAN
    inrange &= is_finite(initial_x);

    if (horizontal_and(inrange)) {
        // fast normal path
        return y;
    }
    else {
        // overflow, underflow and NAN
        r = select(sign_bit(initial_x), 0.-M1, infinite_vec<VTYPE>()); // value in case of overflow or INF
        y = select(inrange, y, r);                                     // +/- overflow
        y = select(is_nan(initial_x), initial_x, y);                   // NAN goes through
        return y;
    }
}
#endif

// instances of exp_d template
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

// Template for exp function, single precision
// The limit of abs(x) is defined by max_x below
// This function does not produce denormals
// Template parameters:
// VTYPE:  float vector type
// BVTYPE: boolean vector type
// M1: 0 for exp, 1 for expm1
// BA: 0 for exp, 1 for 0.5*exp, 2 for pow(2,x), 10 for pow(10,x)

template<class VTYPE, class BVTYPE, int M1, int BA>
static inline VTYPE exp_f(VTYPE const & initial_x) {

    // Taylor coefficients
    const float P0expf   =  1.f/2.f;
    const float P1expf   =  1.f/6.f;
    const float P2expf   =  1.f/24.f;
    const float P3expf   =  1.f/120.f;
    const float P4expf   =  1.f/720.f;
    const float P5expf   =  1.f/5040.f;

    VTYPE  x, r, x2, z, n2;                      // data vectors
    BVTYPE inrange;                              // boolean vector

    // maximum abs(x), value depends on BA, defined below
    // The lower limit of x is slightly more restrictive than the upper limit.
    // We are specifying the lower limit, except for BA = 1 because it is not used for negative x
    float max_x;

    if (BA <= 1) { // exp(x)
        const float ln2f_hi  =  0.693359375f;
        const float ln2f_lo  = -2.12194440e-4f;
        max_x = (BA == 0) ? 87.3f : 89.0f;

        x = initial_x;
        r = round(initial_x*float(VM_LOG2E));
        x = nmul_add(r, VTYPE(ln2f_hi), x);      //  x -= r * ln2f_hi;
        x = nmul_add(r, VTYPE(ln2f_lo), x);      //  x -= r * ln2f_lo;
    }
    else if (BA == 2) {                          // pow(2,x)
        max_x = 126.f;
        r = round(initial_x);
        x = initial_x - r;
        x = x * (float)VM_LN2;
    }
    else if (BA == 10) {                         // pow(10,x)
        max_x = 37.9f;
        const float log10_2_hi = 0.301025391f;   // log10(2) in two parts
        const float log10_2_lo = 4.60503907E-6f;
        x = initial_x;
        r = round(initial_x*float(VM_LOG2E*VM_LN10));
        x = nmul_add(r, VTYPE(log10_2_hi), x);   //  x -= r * log10_2_hi;
        x = nmul_add(r, VTYPE(log10_2_lo), x);   //  x -= r * log10_2_lo;
        x = x * (float)VM_LN10;
    }
    else  {  // undefined value of BA
        return 0.;
    }

    x2 = x * x;
    z = polynomial_5(x,P0expf,P1expf,P2expf,P3expf,P4expf,P5expf);
    z = mul_add(z, x2, x);                       // z *= x2;  z += x;

    if (BA == 1) r--;                            // 0.5 * exp(x)

    // multiply by power of 2
    n2 = vm_pow2n(r);

    if (M1 == 0) {
        // exp
        z = (z + 1.0f) * n2;
    }
    else {
        // expm1
        z = mul_add(z, n2, n2 - 1.0f);           //  z = z * n2 + (n2 - 1.0f);
    }

    // check for overflow
    inrange  = abs(initial_x) < max_x;
    // check for INF and NAN
    inrange &= is_finite(initial_x);

    if (horizontal_and(inrange)) {
        // fast normal path
        return z;
    }
    else {
        // overflow, underflow and NAN
        r = select(sign_bit(initial_x), 0.f-M1, infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        z = select(inrange, z, r);                                      // +/- underflow
        z = select(is_nan(initial_x), initial_x, z);                    // NAN goes through
        return z;
    }
}

// instances of exp_f template
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


#endif  // VECTORMATH_EXP_H
