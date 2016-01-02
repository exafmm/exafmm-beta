#ifndef macros_h
#define macros_h

// Detect SIMD Byte length of architecture
#if __MIC__
const int SIMD_BYTES = 64;                                      //!< SIMD byte length of MIC
#elif __AVX__ | __bgq__
const int SIMD_BYTES = 32;                                      //!< SIMD byte length of AVX and BG/Q
#elif __SSE__ | __sparc_v9__ | _SX
const int SIMD_BYTES = 16;                                      //!< SIMD byte length of SSE, FX, SX
#else
#error no SIMD
#endif

#if _SX
#define __attribute__(x)
#endif

// Use Agner's vectormath for x86 SIMD
#if __MIC__ | __AVX__ | __SSE__
#define EXAFMM_USE_VECTORCLASS 1
#endif

// Bluegene/Q and K computer don't have single precision arithmetic
#if __bgq__ | __sparc_v9__
#ifdef EXAFMM_SINGLE
#error Please use double precision for BG/Q, FX10, FX100
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
#if defined EXAFMM_LAPLACE || EXAFMM_HELMHOLTZ || EXAFMM_STOKES || EXAFMM_BIOTSAVART
#else
#error Please define EXAFMM_LAPLACE or EXAFMM_HELMHOLTZ or EXAFMM_STOKES or EXAFMM_BIOTSAVART
#endif

// Check for mismatching equation and basis
#if EXAFMM_HELMHOLTZ
#if EXAFMM_CARTESIAN
#error Use Spherical for Helmholtz
#endif
#endif

#endif
