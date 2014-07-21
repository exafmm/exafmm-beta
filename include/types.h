#ifndef types_h
#define types_h
#include "align.h"
#include <complex>
#include "kahan.h"
#include "macros.h"
#include <stdint.h>
#include <vector>
#include "vec.h"

// Basic type definitions
#if FP64
typedef double               real_t;                            //!< Floating point type is double precision
const real_t EPS = 1e-16;                                       //!< Double precision epsilon
#else
typedef float                real_t;                            //!< Floating point type is single precision
const real_t EPS = 1e-8;                                        //!< Single precision epsilon
#endif
typedef std::complex<real_t> complex_t;                         //!< Complex type
typedef vec<3,real_t>        vec3;                              //!< Vector of 3 real_t types

// SIMD vector types for MIC, AVX, and SSE
const int NSIMD = SIMD_BYTES / sizeof(real_t);                  //!< SIMD vector length (SIMD_BYTES defined in macros.h)
typedef vec<NSIMD,real_t> simdvec;                              //!< SIMD vector type

// Kahan summation types (Achieves quasi-double precision using single precision types)
#if KAHAN
typedef kahan<real_t>  kreal_t;                                 //!< Floating point type with Kahan summation
typedef vec<4,kreal_t> kvec4;                                   //!< Vector of 4 floats with Kahan summaiton
typedef kahan<simdvec> ksimdvec;                                //!< SIMD vector type with Kahan summation
#else
typedef real_t         kreal_t;                                 //!< Floating point type
typedef vec<4,real_t>  kvec4;                                   //!< Vector of 4 floating point types
typedef simdvec        ksimdvec;                                //!< SIMD vector type
#endif

// Multipole/local expansion coefficients
const int P = EXPANSION;                                        //!< Order of expansions
#if Cartesian
const int NTERM = P*(P+1)*(P+2)/6;                              //!< Number of Cartesian mutlipole/local terms
typedef vec<NTERM,real_t> vecP;                                 //!< Multipole/local coefficient type for Cartesian
#elif Spherical
const int NTERM = P*(P+1)/2;                                    //!< Number of Spherical multipole/local terms
typedef vec<NTERM,complex_t> vecP;                              //!< Multipole/local coefficient type for spherical
#endif

//! Center and radius of bounding box
struct Box {
  vec3   X;                                                     //!< Box center
  real_t R;                                                     //!< Box radius
};

//! Min & max bounds of bounding box
struct Bounds {
  vec3 Xmin;                                                    //!< Minimum value of coordinates
  vec3 Xmax;                                                    //!< Maximum value of coordinates
};

//! Structure of aligned source for SIMD
struct Source {
  vec3   X;                                                     //!< Position
  real_t SRC;                                                   //!< Scalar source values
} __attribute__ ((aligned (16)));

//! Structure of bodies
struct Body : public Source {
  int      IBODY;                                               //!< Initial body numbering for sorting back
  int      IRANK;                                               //!< Initial rank numbering for partitioning back
  uint64_t ICELL;                                               //!< Cell index   
  real_t   WEIGHT;                                              //!< Weight for partitioning
  kvec4    TRG;                                                 //!< Scalar+vector3 target values
};
typedef AlignedAllocator<Body,SIMD_BYTES> BodyAllocator;        //!< Body alignment allocator
typedef std::vector<Body,BodyAllocator>   Bodies;               //!< Vector of bodies
//typedef std::vector<Body>                 Bodies;               //!< Vector of bodies
typedef Bodies::iterator                  B_iter;               //!< Iterator of body vector

//! Structure of cells
struct Cell {
  int       IPARENT;                                            //!< Index of parent cell
  int       ICHILD;                                             //!< Index of first child cell
  int       NCHILD;                                             //!< Number of child cells
  int       IBODY;                                              //!< Index of first body
  int       NBODY;                                              //!< Number of descendant bodies
  B_iter    BODY;                                               //!< Iterator of first body
  uint64_t  ICELL;                                              //!< Cell index
  real_t    WEIGHT;                                             //!< Weight for partitioning
  vec3      X;                                                  //!< Cell center
  real_t    R;                                                  //!< Cell radius
  vecP      M;                                                  //!< Multipole coefficients
  vecP      L;                                                  //!< Local coefficients
};
typedef std::vector<Cell> Cells;                                //!< Vector of cells
typedef Cells::iterator   C_iter;                               //!< Iterator of cell vector

#endif
