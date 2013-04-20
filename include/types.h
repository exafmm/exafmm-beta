#ifndef types_h
#define types_h
#include "align.h"
#include <complex>
#include "kahan.h"
#include "macros.h"
#include <vector>
#include "vec.h"

// Basic type definitions
#if FP64
typedef double               real_t;                            //!< Floating point type is double precision
#else
typedef float                real_t;                            //!< Floating point type is single precision
#endif
typedef std::complex<real_t> complex_t;                         //!< Complex type
typedef vec<3,real_t>        vec3;                              //!< Vector of 3 floating point types
typedef vec<3,float>         fvec3;                             //!< Vector of 3 single precision types

// SIMD vector types for MIC, AVX, and SSE
const int NSIMD = SIMD_BYTES / sizeof(real_t);                  //!< SIMD vector length
typedef vec<NSIMD,real_t> simdvec;                              //!< SIMD vector type

// Kahan summation types
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
const real_t EPS2 = 0.0;                                        //!< Softening parameter (squared)
#if COMkernel
const int MTERM = P*(P+1)*(P+2)/6-3;                            //!< Number of Cartesian mutlipole terms
#else
const int MTERM = P*(P+1)*(P+2)/6;                              //!< Number of Cartesian mutlipole terms
#endif
const int LTERM = (P+1)*(P+2)*(P+3)/6;                          //!< Number of Cartesian local terms
const int NTERM = P*(P+1)/2;                                    //!< Number of Spherical multipole/local terms
#if Cartesian
typedef vec<MTERM,real_t> vecM;                                 //!< Multipole coefficient type for Cartesian
typedef vec<LTERM,real_t> vecL;                                 //!< Local coefficient type for Cartesian
#elif Spherical
typedef vec<NTERM,complex_t> vecM;                              //!< Multipole coefficient type for spherical
typedef vec<NTERM,complex_t> vecL;                              //!< Local coefficient type for spherical
#endif

//! Structure for defining bounding box
struct Box {
  vec3   X;                                                     //!< Box center
  real_t R;                                                     //!< Box radius
};

//! Structure of aligned source for SIMD
struct Source {
  vec3   X;                                                     //!< Position
  real_t SRC;                                                   //!< Scalar source values
} __attribute__ ((aligned (16)));

//! Structure of bodies
struct Body : public Source {
  int    IBODY;                                                 //!< Initial body numbering for sorting back
  int    IPROC;                                                 //!< Initial process numbering for partitioning back
  int    ICELL;                                                 //!< Cell index
  kvec4  TRG;                                                   //!< Scalar+vector3 target values
};
typedef AlignedAllocator<Body,SIMD_BYTES>         BodyAllocator;//!< Body alignment allocator
typedef std::vector<Body,BodyAllocator>           Bodies;       //!< Vector of bodies
typedef std::vector<Body,BodyAllocator>::iterator B_iter;       //!< Iterator of body vector

//! Structure of cells
struct Cell {
  int       NCHILD;                                             //!< Number of child cells
  int       NCBODY;                                             //!< Number of child bodies
  int       NDBODY;                                             //!< Number of descendant bodies
  int       PARENT;                                             //!< Index of parent cell
  int       CHILD;                                              //!< Index of child cells
  long long ICELL;                                              //!< Cell index
  B_iter    BODY;                                               //!< Iterator of first body
  vec3      X;                                                  //!< Cell center
  real_t    R;                                                  //!< Cell radius
  real_t    RMAX;                                               //!< Max cell radius
  real_t    RCRIT;                                              //!< Critical cell radius
  vecM      M;                                                  //!< Multipole coefficients
  vecL      L;                                                  //!< Local coefficients
};
typedef std::vector<Cell>           Cells;                      //!< Vector of cells
typedef std::vector<Cell>::iterator C_iter;                     //!< Iterator of cell vector

#endif
