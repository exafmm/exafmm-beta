#ifndef types_h
#define types_h

#ifdef __INTEL_COMPILER
#pragma warning(disable:193 383 444 981 1572 2259)
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <pthread.h>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"

#if SIMDIZATION
typedef enum {
  simdize_none,
  simdize_sse,
  simdize_avx
} simdize_option;

#include <immintrin.h>
#endif

#if CMDLINE_ARGS
#include "exafmm_config.h"
#endif

#if SSE
#include <xmmintrin.h>
#endif

#if 0
#if _OPENMP
#include <omp.h>
#else
int omp_get_thread_num() { return 0; }
#endif
#endif

#if 0
#if MTHREADS
#include <task_group.h>
#endif

#else

//#include <common.cilkh>
#include "task_parallel.h"

#endif

#if PAPI
#include <papi.h>
#endif

typedef float       real_t;                                     //!< Real number type on CPU
typedef vec<3,real_t> vec3;                                     //!< 3-D vector type

#ifndef KERNEL
int MPIRANK    = 0;                                             //!< MPI comm rank
int MPISIZE    = 1;                                             //!< MPI comm size
int IMAGES     = 0;                                             //!< Number of periodic image sublevels
real_t THETA   = .5;                                            //!< Multipole acceptance criteria
vec3 Xperiodic = .0;                                            //!< Coordinate offset of periodic image
#if PAPI
int PAPIEVENT  = PAPI_NULL;                                     //!< PAPI event handle
#endif
#else
extern int MPIRANK;                                             //!< MPI comm rank
extern int MPISIZE;                                             //!< MPI comm size
extern int IMAGES;                                              //!< Number of periodic image sublevels
extern real_t THETA;                                            //!< Multipole acceptance criteria
extern vec3 Xperiodic;                                          //!< Coordinate offset of periodic image
#if PAPI
extern int PAPIEVENT;                                           //!< PAPI event handle
#endif
#endif

const int    P      = 3;                                        //!< Order of expansions
#if CMDLINE_ARGS
int    NCRIT  = 10;
#else
const int    NCRIT  = 10;                                       //!< Number of bodies per cell
#endif
const real_t EPS    = 1e-6;                                     //!< Single precision epsilon
const real_t EPS2   = .0;                                       //!< Softening parameter (squared)

#if COMkernel
const int MTERM = P*(P+1)*(P+2)/6-3;                            //!< Number of Cartesian mutlipole terms
#else
const int MTERM = P*(P+1)*(P+2)/6;                              //!< Number of Cartesian mutlipole terms
#endif
const int LTERM = (P+1)*(P+2)*(P+3)/6;                          //!< Number of Cartesian local terms

typedef vec<MTERM,real_t>  vecM;                                //!< Multipole coefficient type for Cartesian
typedef vec<LTERM,real_t>  vecL;                                //!< Local coefficient type for Cartesian

//! Structure for pthread based trace
struct Trace {
  pthread_t thread;                                             //!< pthread id
  double    begin;                                              //!< Begin timer of trace
  double    end;                                                //!< End timer of trace
  int       color;                                              //!< Color of trace
};
typedef std::map<pthread_t,double>             ThreadTrace;     //!< Map of pthread id to traced value
typedef std::map<pthread_t,int>                ThreadMap;       //!< Map of pthread id to thread id
typedef std::queue<Trace>                      Traces;          //!< Queue of traces
typedef std::map<std::string,double>           Timer;           //!< Map of timer event name to timed value
typedef std::map<std::string,double>::iterator TI_iter;         //!< Iterator for timer event name map

enum Equation {                                                 //!< Equation type enumeration
  Laplace,                                                      //!< Laplace equation
  Yukawa,                                                       //!< Yukawa equation
  Helmholtz,                                                    //!< Helmholtz equation
  Stokes,                                                       //!< Stokes equation
  VanDerWaals                                                   //!< Van der Walls equation
};

//! Structure of source bodies (stuff to send)
struct Body {
  int  IBODY;                                                   //!< Initial body numbering for sorting back
  int  IPROC;                                                   //!< Initial process numbering for partitioning back
  int  ICELL;                                                   //!< Cell index
  vec3 X;                                                       //!< Position
  real_t SRC;                                                   //!< Scalar source values
  vec<4,real_t> TRG;                                            //!< Scalar+vector target values
};
typedef std::vector<Body>            Bodies;                    //!< Vector of bodies
typedef std::vector<Body>::iterator  B_iter;                    //!< Iterator for body vector

//! Linked list of leafs (only used in fast/topdown.h)
struct Leaf {
  int   I;                                                      //!< Unique index for every leaf
  vec3  X;                                                      //!< Coordinate of leaf
  Leaf *NEXT;                                                   //!< Pointer to next leaf
};
typedef std::vector<Leaf>           Leafs;                      //!< Vector of leafs
typedef std::vector<Leaf>::iterator L_iter;                     //!< Iterator for leaf vector

//! Structure of nodes (only used in fast/topdown.h)
struct Node {
  bool  NOCHILD;                                                //!< Flag for twig nodes
  int   LEVEL;                                                  //!< Level in the tree structure
  int   NLEAF;                                                  //!< Number of descendant leafs
  int   CHILD[8];                                               //!< Index of child node
  vec3  X;                                                      //!< Coordinate at center
  Leaf *LEAF;                                                   //!< Pointer to first leaf
};
typedef std::vector<Node>           Nodes;                      //!< Vector of nodes
typedef std::vector<Node>::iterator N_iter;                     //!< Iterator for node vector

//! Structure of cells
struct Cell {
  int       NCHILD;                                             //!< Number of child cells
  int       NCLEAF;                                             //!< Number of child leafs
  int       NDLEAF;                                             //!< Number of descendant leafs
  int       PARENT;                                             //!< Iterator offset of parent cell
  int       CHILD;                                              //!< Iterator offset of child cells
  long long ICELL;                                              //!< Cell index
  B_iter    LEAF;                                               //!< Iterator of first leaf
  vec3      X;                                                  //!< Cell center
  real_t    R;                                                  //!< Cell radius
  real_t    RMAX;                                               //!< Max cell radius
  real_t    RCRIT;                                              //!< Critical cell radius
  vecM      M;                                                  //!< Multipole coefficients
  vecL      L;                                                  //!< Local coefficients
};
typedef std::vector<Cell>           Cells;                      //!< Vector of cells
typedef std::vector<Cell>::iterator C_iter;                     //!< Iterator for cell vector
typedef std::queue<C_iter>          CellQueue;                  //!< Queue of cell iterators
typedef std::pair<C_iter,C_iter>    Pair;                       //!< Pair of interacting cells
typedef std::deque<Pair>            PairQueue;                  //!< Queue of interacting cell pairs

//! Structure for Ewald summation
struct Ewald {
  vec3   K;                                                     //!< 3-D wave number vector
  real_t REAL;                                                  //!< real part of wave
  real_t IMAG;                                                  //!< imaginary part of wave
};
typedef std::vector<Ewald>           Ewalds;                    //!< Vector of Ewald summation types
typedef std::vector<Ewald>::iterator E_iter;                    //!< Iterator for Ewald summation types

#endif
