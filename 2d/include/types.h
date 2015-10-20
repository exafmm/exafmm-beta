#ifndef types_h
#define types_h
#include "align.h"
#include <complex>
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
typedef vec<2,real_t>        vec2;                              //!< Vector of 3 floating point types
typedef vec<2,float>         fvec2;                             //!< Force float (Used only for communication)

// Multipole/local expansion coefficients
const int P = EXPANSION;                                        //!< Order of expansions
typedef vec<P,complex_t> vecP;                                  //!< Multipole/local coefficient type

//! Structures for defining bounding box
struct Box {
  vec2   X;                                                     //!< Box center
  real_t R;                                                     //!< Box radius
};
struct Bounds {
  vec2 Xmin;                                                    //!< Minimum value of coordinates
  vec2 Xmax;                                                    //!< Maximum value of coordinates
};

//! Structure of aligned source for SIMD
struct Source {
  vec2   X;                                                     //!< Position
  real_t SRC;                                                   //!< Scalar source values
} __attribute__ ((aligned (16)));

//! Structure of bodies
struct Body : public Source {
  int    IBODY;                                                 //!< Initial body numbering for sorting back
  int    IPROC;                                                 //!< Initial process numbering for partitioning back
  int    ICELL;                                                 //!< Cell index
  real_t TRG;                                                   //!< Scalar+vector3 target values
};
typedef AlignedAllocator<Body,SIMD_BYTES> BodyAllocator;        //!< Body alignment allocator
//typedef std::vector<Body,BodyAllocator>   Bodies;               //!< Vector of bodies
typedef std::vector<Body>                 Bodies;               //!< Vector of bodies
typedef Bodies::iterator                  B_iter;               //!< Iterator of body vector

//! Structure of cells
struct Cell {
  int       NCHILD;                                             //!< Number of child cells
  int       NCBODY;                                             //!< Number of child bodies
  int       NDBODY;                                             //!< Number of descendant bodies
  int       PARENT;                                             //!< Index of parent cell
  int       CHILD;                                              //!< Index of child cells
  long long ICELL;                                              //!< Cell index
  B_iter    BODY;                                               //!< Iterator of first body
  vec2      X;                                                  //!< Cell center
  real_t    R;                                                  //!< Cell radius
  real_t    RCRIT;                                              //!< Critical cell radius
  vecP      M;                                                  //!< Multipole coefficients
  vecP      L;                                                  //!< Local coefficients
};
typedef std::vector<Cell>           Cells;                      //!< Vector of cells
typedef std::vector<Cell>::iterator C_iter;                     //!< Iterator of cell vector
#endif
