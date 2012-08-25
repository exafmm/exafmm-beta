#ifndef macros_h
#define macros_h

// Intel compiler warning disable
#ifdef __INTEL_COMPILER
#pragma warning(disable:193 383 444 981 1572 2259)
#endif

// Override assertion
#if ASSERT
#include <cassert>
#else
#define assert(x)
#endif

// SIMD instruction
#if __AVX__
#include <immintrin.h>
#elif __SSE4_2__
#include <nmmintrin.h>
#elif __SSE4_1__
#include <smmintrin.h>
#elif __SSSE3__
#include <tmmintrin.h>
#elif __SSE3__
#include <pmmintrin.h>
#elif __SSE2__
#include <emmintrin.h>
#elif __SSE__
#include <xmmintrin.h>
#endif

#endif
