#ifndef macros_h
#define macros_h

// Disable a few Intel compiler warnings
#ifdef __INTEL_COMPILER
#pragma warning(disable:161 193 383 444 981 1572 2259)
#endif

// Override assertion
#if EXAFMM_ASSERT
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
#ifdef EXAFMM_SINGLE
#error Please use double precision for BG/P, BG/Q, FX10, FX100
#endif
#endif

// Check for equation and basis
#ifndef EXAFMM_EXPANSION
#error EXAFMM_EXPANSION undefined
#endif
#if defined EXAFMM_CARTESIAN || EXAFMM_SPHERICAL
#else
#error Please define EXAFMM_CARTESIAN or EXAFMM_SPHERICAL
#endif
#if defined EXAFMM_LAPLACE || EXAFMM_HELMHOLTZ || EXAFMM_STOKES
#else
#error Please define EXAFMM_LAPLACE or EXAFMM_HELMHOLTZ or EXAFMM_STOKES
#endif

// Check for mismatching equation and basis
#if EXAFMM_HELMHOLTZ
#if EXAFMM_CARTESIAN
#error Use Spherical for Helmholtz
#endif
#endif

#endif
