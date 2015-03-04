#ifndef macros_h
#define macros_h

// Disable a few Intel compiler warnings
#ifdef __INTEL_COMPILER
#pragma warning(disable:161 193 383 444 981 1572 2259)
#endif

// Override assertion
#if ASSERT
#include <cassert>
#else
#define assert(x)
#endif

// Detect SIMD Byte length of architecture
#if __MIC__
const int SIMD_BYTES = 64;                                      //!< SIMD byte length of MIC
#elif __AVX__ | __bgq__
const int SIMD_BYTES = 32;                                      //!< SIMD byte length of AVX and BG/Q
#elif __SSE__ | __bgp__ | __sparc_v9__
const int SIMD_BYTES = 16;                                      //!< SIMD byte length of SSE and BG/P
#else
#error no SIMD
#endif

// Bluegene/Q and K computer don't have single precision arithmetic
#if __bgp__ | __bgq__ | __sparc_v9__
#ifndef FP64
#error Please use FP64 for BG/P, BG/Q, and K computer
#endif
#endif

#endif
