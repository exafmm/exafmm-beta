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

#endif
