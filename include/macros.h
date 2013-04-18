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

// Detect SIMD Byte length of architecture
#if __MIC__
const int SIMD_BYTES = 64;                                      //!< SIMD byte length of MIC
#elif __AVX__
const int SIMD_BYTES = 32;                                      //!< SIMD byte length of AVX
#elif __SSE__
const int SIMD_BYTES = 16;                                      //!< SIMD byte length of SSE
#endif

#endif
