#ifndef types_h
#define types_h
#include <assert.h>                                             // Some compilers don't have cassert
#include <complex>
#include "kahan.h"
#include "macros.h"
#include <stdint.h>
#include <vector>
#include "vec.h"

namespace exafmm {
  // Basic type definitions
#if EXAFMM_SINGLE
  typedef float real_t;                                         //!< Floating point type is single precision
  const real_t EPS = 1e-8;                                      //!< Single precision epsilon
#else
  typedef double real_t;                                        //!< Floating point type is double precision
  const real_t EPS = 1e-16;                                     //!< Double precision epsilon
#endif
  typedef std::complex<real_t> complex_t;                       //!< Complex type
  const complex_t I(0.,1.);                                     //!< Imaginary unit

  typedef vec<3,int> ivec3;                                     //!< Vector of 3 int types
  typedef vec<3,real_t> vec3;                                   //!< Vector of 3 real_t types
  typedef vec<4,real_t> vec4;                                   //!< Vector of 4 real_t types
  typedef vec<3,complex_t> cvec3;                               //!< Vector of 3 complex_t types

  // SIMD vector types for AVX512, AVX, and SSE
  const int NSIMD = SIMD_BYTES / sizeof(real_t);                //!< SIMD vector length (SIMD_BYTES defined in macros.h)
  typedef vec<NSIMD,real_t> simdvec;                            //!< SIMD vector type

  // Kahan summation types (Achieves quasi-double precision using single precision types)
#if EXAFMM_USE_KAHAN
  typedef kahan<real_t> kreal_t;                                //!< Real type with Kahan summation
  typedef kahan<complex_t> kcomplex_t;                          //!< Complex type with Kahan summation
  typedef kahan<simdvec> ksimdvec;                              //!< SIMD vector type with Kahan summation
#else
  typedef real_t kreal_t;                                       //!< Real type (dummy Kahan)
  typedef complex_t kcomplex_t;                                 //!< Complex type (dummy Kahan)
  typedef simdvec ksimdvec;                                     //!< SIMD vector type (dummy Kahan)
#endif
  typedef vec<4,kreal_t> kvec4;                                 //!< Vector of 4 real types with Kahan summaiton
  typedef vec<4,kcomplex_t> kcvec4;                             //!< Vector of 4 complex types with Kahan summaiton

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

  //! Equations supported
  enum Equation {
    Empty,                                                      //!< Empty kernel
    Laplace,                                                    //!< Laplace kernel
    Helmholtz,                                                  //!< Helmholtz kernel
    BiotSavart                                                  //!< Biot-Savart kernel
  };

  //! Structure of aligned source for SIMD
  template<Equation equation=Empty>
  struct Source {                                               //!< Base components of source structure
    vec3      X;                                                //!< Position
  };
  template<>
  struct Source<Laplace> : public Source<> {                    //!< Special components for Laplace
    real_t    SRC;                                              //!< Scalar real values
  };
  template<>
  struct Source<Helmholtz> : public Source<> {                  //!< Special components for Helmholtz
    complex_t SRC;                                              //!< Scalar complex values
  };
  template<>
  struct Source<BiotSavart> : public Source<> {                 //!< Special components for Biot-Savart
    vec4      SRC;                                              //!< Vector real values
  };

  //! Structure of bodies
  template<Equation equation=Empty>
  struct Body {                                                 //!< Base components of body structure
    int     IBODY;                                              //!< Initial body numbering for sorting back
    int     IRANK;                                              //!< Initial rank numbering for partitioning back
    int64_t ICELL;                                              //!< Cell index   
    real_t  WEIGHT;                                             //!< Weight for partitioning
  };
  template<>
  struct Body<Laplace> : public Source<Laplace>, Body<> {       //!< Special components for Laplace
    kvec4   TRG;                                                //!< Scalar+vector3 real values
  };
  template<>
  struct Body<Helmholtz> : public Source<Helmholtz>, Body<> {   //!< Special components for Helmholtz
    kcvec4  TRG;                                                //!< Scalar+vector3 complex values
  };
  template<>
  struct Body<BiotSavart> : public Source<BiotSavart>, Body<> { //!< Special components for Biot-Savart
    kvec4   TRG;                                                //!< Scalar+vector3 real values
  };

  // Multipole/local expansion coefficients
  const int P = EXAFMM_EXPANSION;                               //!< Order of expansions
  const int NTERM_LC = P*(P+1)*(P+2)/6;                         //!< # of terms for Lapalace Cartesian 
  const int NTERM_LS = P*(P+1)/2;                               //!< # of terms for Laplace Spherical
  const int NTERM_HS = P*P;                                     //!< # of terms for Helmholtz Spherical
  const int NTERM_BS = 3*P*(P+1)/2;                             //!< # of terms for Biot-Savart Spherical
  typedef vec<NTERM_LC,real_t> vecLC;                           //!< Coef vector for Laplace Cartesian
  typedef vec<NTERM_LS,complex_t> vecLS;                        //!< Coef vector for Laplace Spherical
  typedef vec<NTERM_HS,complex_t> vecHS;                        //!< Coef vector for Helmholtz Spherical
  typedef vec<NTERM_BS,complex_t> vecBS;                        //!< Coef vector for Biot-Savart Spherical

  //! Structure of cells
  template<Equation equation=Empty>
  struct Cell {                                                 //!< Base components of cell structure
    int      IPARENT;                                           //!< Index of parent cell
    int      ICHILD;                                            //!< Index of first child cell
    int      NCHILD;                                            //!< Number of child cells
    int      IBODY;                                             //!< Index of first body
    int      NBODY;                                             //!< Number of descendant bodies
#if EXAFMM_COUNT_LIST
    int      numP2P;                                            //!< Size of P2P interaction list per cell
    int      numM2L;                                            //!< Size of M2L interaction list per cell
#endif
    uint64_t ICELL;                                             //!< Cell index
    real_t   WEIGHT;                                            //!< Weight for partitioning
    real_t   SCALE;                                             //!< Scale for Helmholtz kernel
    vec3     X;                                                 //!< Cell center
    real_t   R;                                                 //!< Cell radius
  };
  template<>
  struct Cell<Laplace> : public Cell<> {                        //!< Special components for Laplace
    typedef std::vector<Body<Laplace> >::iterator B_iter;       //!< Iterator type for body vector
    B_iter BODY;                                                //!< Iterator of first body
#if EXAFMM_CARTESIAN
    vecLC  M, L;                                                //!< Multipole/local coefficients
#elif EXAFMM_SPHERICAL
    vecLS  M, L;                                                //!< Multipole/local coefficients
#endif
  };
  template<>
  struct Cell<Helmholtz> : public Cell<> {                      //!< Special components for Helmholtz
    typedef std::vector<Body<Helmholtz> >::iterator B_iter;     //!< Iterator type for body vector
    B_iter BODY;                                                //!< Iterator of first body
    vecHS  M, L;                                                //!< Multipole/local coefficients
  };
  template<>
  struct Cell<BiotSavart> : public Cell<> {                     //!< Special components for Biot-Savart
    typedef std::vector<Body<BiotSavart> >::iterator B_iter;    //!< Iterator type for body vector
    B_iter BODY;                                                //!< Iterator of first body
    vecBS  M, L;                                                //!< Multipole/local coefficients
  };

  struct Kernel {
    static vec3 Xperiodic;                                      //!< Periodic coordinate offset
    static real_t eps2;                                         //!< Epslion squared
    static complex_t wavek;                                     //!< Helmholtz wave number
  };
}
#endif
