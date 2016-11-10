/****************************  vectorclass.h   ********************************
* Author:        Agner Fog
* Date created:  2012-05-30
* Last modified: 2015-08-08
* Version:       1.18
* Project:       vector classes
* Description:
* Header file defining vector classes as interface to intrinsic functions
* in x86 microprocessors with SSE2 and later instruction sets up to AVX512.
*
* Instructions:
* Use Gnu, Clang, Intel or Microsoft C++ compiler. Compile for the desired
* instruction set, which must be at least SSE2. Specify the supported
* instruction set by a command line define, e.g. __SSE4_1__ if the
* compiler does not automatically do so.
*
* Each vector object is represented internally in the CPU as a vector
* register with 128, 256 or 512 bits.
*
* This header file includes the appropriate header files depending on the
* supported instruction set
*
* For detailed instructions, see VectorClass.pdf
*
* (c) Copyright 2012 - 2015 GNU General Public License www.gnu.org/licenses
******************************************************************************/
#ifndef VECTORCLASS_H
#define VECTORCLASS_H  116

// Maximum vector size, bits. Allowed values are 128, 256, 512
#ifndef MAX_VECTOR_SIZE
#define MAX_VECTOR_SIZE 512
#endif

#include <x86intrin.h>
#include <stdint.h>
#include <stdlib.h>

#define __FMA__  1
#define INSTRSET 9
#define GCC_VERSION  ((__GNUC__) * 10000 + (__GNUC_MINOR__) * 100 + (__GNUC_PATCHLEVEL__))

template <int32_t  n> class Const_int_t  {};
template <uint32_t n> class Const_uint_t {};
#define const_int(n)  (Const_int_t <n>())
#define const_uint(n) (Const_uint_t<n>())

template <bool> class Static_error_check {
public:  Static_error_check(){};
};
template <> class Static_error_check<false> {
private: Static_error_check(){};
};

//#include "vectori128.h"
//#include "vectorf128.h"
#if __AVX__
#include "vectori256.h"
#include "vectorf256.h"
#endif
#if __AVX512F__ | __MIC__
#include "vectori512.h"
#include "vectorf512.h"
#endif

#endif  // VECTORCLASS_H
