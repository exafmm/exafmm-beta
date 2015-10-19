#ifndef types_h
#define types_h
#include "align.h"
#include <complex>
#include "kahan.h"
#include "macros.h"
#include <stdint.h>
#include <vector>
#include "vec.h"

namespace exafmm {
  // Basic type definitions
#if EXAFMM_SINGLE
  typedef float                real_t;                          //!< Floating point type is single precision
  const real_t EPS = 1e-8;                                      //!< Single precision epsilon
#else
  typedef double               real_t;                          //!< Floating point type is double precision
  const real_t EPS = 1e-16;                                     //!< Double precision epsilon
#endif
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  typedef vec<3,int>           ivec3;                           //!< Vector of 3 int types
  typedef vec<3,real_t>        vec3;                            //!< Vector of 3 real_t types
  typedef vec<3,complex_t>     cvec3;                           //!< Vector of 3 complex_t types
  typedef vec<4,complex_t>     cvec4;                           //!< Vector of 4 complex_t types

  // SIMD vector types for MIC, AVX, and SSE
  const int NSIMD = SIMD_BYTES / sizeof(real_t);                //!< SIMD vector length (SIMD_BYTES defined in macros.h)
  typedef vec<NSIMD,real_t>     simdvec;                        //!< SIMD vector type
  typedef std::complex<simdvec> csimdvec;                       //!< Complex SIMD vector type

  // Kahan summation types (Achieves quasi-double precision using single precision types)
#if EXAFMM_USE_KAHAN
  typedef kahan<real_t>   kreal_t;                              //!< Floating point type with Kahan summation
  typedef vec<4,kreal_t>  kvec4;                                //!< Vector of 4 floats with Kahan summaiton
  typedef kahan<simdvec>  ksimdvec;                             //!< SIMD vector type with Kahan summation
  typedef kahan<csimdvec> kcsimdvec;                            //!< Complex SIMD vector type with Kahan summation
#else
  typedef real_t          kreal_t;                              //!< Floating point type
  typedef vec<4,real_t>   kvec4;                                //!< Vector of 4 floating point types
  typedef simdvec         ksimdvec;                             //!< SIMD vector type
  typedef csimdvec        kcsimdvec;                            //!< Complex SIMD vector type
#endif

  // Multipole/local expansion coefficients
  const int P = EXAFMM_EXPANSION;                               //!< Order of expansions
#if EXAFMM_CARTESIAN
  const int NTERM = P*(P+1)*(P+2)/6;                            //!< Number mutlipole/local terms
  typedef vec<NTERM,real_t> vecP;                               //!< Multipole/local coefficient type
#elif EXAFMM_SPHERICAL
#if EXAFMM_LAPLACE
  const int NTERM = P*(P+1)/2;                                  //!< Number of multipole/local terms
#elif EXAFMM_HELMHOLTZ
  const int NTERM = P*P;                                        //!< Number of multipole/local terms
#endif
  typedef vec<NTERM,complex_t> vecP;                            //!< Multipole/local coefficient type
#endif

  //! Center and radius of bounding box
  struct Box {
    vec3   X;                                                   //!< Box center
    real_t R;                                                   //!< Box radius
  };

  //! Min & max bounds of bounding box
  struct Bounds {
    vec3 Xmin;                                                  //!< Minimum value of coordinates
    vec3 Xmax;                                                  //!< Maximum value of coordinates
  };

  //! Structure of aligned source for SIMD
  struct Source {
    vec3   X;                                                   //!< Position
#if EXAFMM_LAPLACE
    real_t SRC;                                                 //!< Scalar source values
#elif EXAFMM_HELMHOLTZ
    complex_t SRC;                                              //!< Scalar source values
#endif
  } __attribute__ ((aligned (16)));

  //! Structure of bodies
  struct Body : public Source {
    int      IBODY;                                             //!< Initial body numbering for sorting back
    int      IRANK;                                             //!< Initial rank numbering for partitioning back
    uint64_t ICELL;                                             //!< Cell index   
    real_t   WEIGHT;                                            //!< Weight for partitioning
#if EXAFMM_LAPLACE
    kvec4    TRG;                                               //!< Scalar+vector3 target values
#elif EXAFMM_HELMHOLTZ
    cvec4    TRG;                                               //!< Scalar+vector3 target values
#endif
  };
  typedef AlignedAllocator<Body,SIMD_BYTES> BodyAllocator;      //!< Body alignment allocator
  typedef std::vector<Body,BodyAllocator>   Bodies;             //!< Vector of bodies
  typedef Bodies::iterator                  B_iter;             //!< Iterator of body vector

  //! Structure of cells
  struct Cell {
    int       IPARENT;                                          //!< Index of parent cell
    int       ICHILD;                                           //!< Index of first child cell
    int       NCHILD;                                           //!< Number of child cells
    int       IBODY;                                            //!< Index of first body
    int       NBODY;                                            //!< Number of descendant bodies
#if EXAFMM_COUNT_LIST
    int       numP2P;                                           //!< Size of P2P interaction list per cell
    int       numM2L;                                           //!< Size of M2L interaction list per cell
#endif
    B_iter    BODY;                                             //!< Iterator of first body
    uint64_t  ICELL;                                            //!< Cell index
    real_t    WEIGHT;                                           //!< Weight for partitioning
    real_t    SCALE;                                            //!< Scale for Helmholtz kernel
    vec3      X;                                                //!< Cell center
    real_t    R;                                                //!< Cell radius
    vecP      M;                                                //!< Multipole coefficients
    vecP      L;                                                //!< Local coefficients
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells
  typedef Cells::iterator   C_iter;                             //!< Iterator of cell vector
}
#endif
