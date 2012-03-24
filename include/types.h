/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef types_h
#define types_h

#ifdef __INTEL_COMPILER
#pragma warning(disable:193 383 444 981 1572 2259)
#endif

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <omp.h>
#include <queue>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"                                                //!< My vector type with operator overloading
#if PAPI
#include <papi.h>
#endif
#if QUARK
#include "quark.h"
#endif

typedef unsigned           bigint;                              //!< Big integer type
typedef float              real;                                //!< Real number type on CPU
typedef float              gpureal;                             //!< Real number type on GPU
typedef std::complex<real> complex;                             //!< Complex number type
typedef vec<3,real>        vect;                                //!< 3-D vector type


#ifndef KERNEL
int MPIRANK    = 0;                                             //!< MPI comm rank
int MPISIZE    = 1;                                             //!< MPI comm size
int DEVICE     = 0;                                             //!< GPU device ID
int IMAGES     = 0;                                             //!< Number of periodic image sublevels
real THETA     = .5;                                            //!< Multipole acceptance criteria
vect Xperiodic = 0;                                             //!< Coordinate offset of periodic image
#if PAPI
int PAPIEVENT  = PAPI_NULL;                                     //!< PAPI event handle
#endif
#else
extern int MPIRANK;                                             //!< MPI comm rank
extern int MPISIZE;                                             //!< MPI comm size
extern int DEVICE;                                              //!< GPU device ID
extern int IMAGES;                                              //!< Number of periodic image sublevels
extern real THETA;                                              //!< Multipole acceptance criteria
extern vect Xperiodic;                                          //!< Coordinate offset of periodic image
#if PAPI
extern int PAPIEVENT;                                           //!< PAPI event handle
#endif
#endif

const int  P        = 3;                                        //!< Order of expansions
const int  NCRIT    = 10;                                       //!< Number of bodies per cell
const int  MAXBODY  = 200000;                                   //!< Maximum number of bodies per GPU kernel
const int  MAXCELL  = 10000000;                                 //!< Maximum number of bodies/coefs in cell per GPU kernel
const real CLET     = 2;                                        //!< LET opening critetia
const real EPS      = 1e-6;                                     //!< Single precision epsilon
const real EPS2     = 0;                                        //!< Softening parameter (squared)
const real R2MIN    = 0.25;                                     //!< Minimum value for L-J R^2
const real R2MAX    = 64;                                       //!< Maximum value for L-J R^2
const int  GPUS     = 3;                                        //!< Number of GPUs per node
const int  THREADS  = 64;                                       //!< Number of threads per thread-block
const int  PTHREADS = 4;                                        //!< Number of pthreads in quark

const int MTERM = P*(P+1)*(P+2)/6;                              //!< Number of Cartesian mutlipole terms
const int LTERM = (P+1)*(P+2)*(P+3)/6;                          //!< Number of Cartesian local terms
const int NTERM = P*(P+1)/2;                                    //!< Number of Spherical multipole/local terms

#if SPHERICAL
#if STOKES
typedef vec<4*NTERM,complex>                     Mset;            //!< Multipole coefficient type for spherical
typedef vec<4*NTERM,complex>                     Lset;            //!< Local coefficient type for spherical
#else
typedef vec<NTERM,complex>                     Mset;            //!< Multipole coefficient type for spherical
typedef vec<NTERM,complex>                     Lset;            //!< Local coefficient type for spherical
#endif
#else
typedef vec<MTERM,real>                        Mset;            //!< Multipole coefficient type for Cartesian
typedef vec<LTERM,real>                        Lset;            //!< Local coefficient type for Cartesian
#endif
typedef std::vector<bigint>                    Bigints;         //!< Vector of big integer types

//! Structure for pthread based trace
struct Trace {
  pthread_t thread;
  double    begin;
  double    end;
  int       color;
};
typedef std::map<pthread_t,double>             ThreadTrace;     //!< Map of pthread id to traced value
typedef std::map<pthread_t,int>                ThreadMap;       //!< Map of pthread id to thread id
typedef std::queue<Trace>                      Traces;          //!< Queue of traces
typedef std::map<std::string,double>           Timer;           //!< Map of timer event name to timed value
typedef std::map<std::string,double>::iterator TI_iter;         //!< Iterator for timer event name map

enum Equation {                                                 //!< Equation type enumeration
  Laplace,                                                      //!< Laplace potential + force
  VanDerWaals,                                                  //!< Van der Walls potential + force
  Stokes                                                        //!< Stokes kernels
};

//! Structure of source bodies (stuff to send)
struct JBody {
  int         IBODY;                                            //!< Initial body numbering for sorting back
  int         IPROC;                                            //!< Initial process numbering for sending back
  bigint      ICELL;                                            //!< Cell index
  vect        X;                                                //!< Position
#if STOKES
  vect        FORCE;
#endif  
  real        SRC;                                              //!< Scalar source values
  JBody() :
#if STOKES
  FORCE(0),
#endif
  IBODY(0), IPROC(0), ICELL(0), X(0), SRC(0) {}       //!< Constructor
};
typedef std::vector<JBody>             JBodies;                 //!< Vector of source bodies
typedef std::vector<JBody>::iterator   JB_iter;                 //!< Iterator for source body vector

//! Structure of bodies
struct Body : public JBody {
  vec<4,real> TRG;                                              //!< Scalar+vector target values
  bool operator<(const Body &rhs) const {                       //!< Overload operator for comparing body index
    return this->IBODY < rhs.IBODY;                             //!< Comparison function for body index
  }                                                             //!< End operator overload
  Body() : TRG(0) {}                                            //!< Constructor
};
typedef std::vector<Body>              Bodies;                  //!< Vector of bodies
typedef std::vector<Body>::iterator    B_iter;                  //!< Iterator for body vector

//! Linked list of leafs (only used in fast/topdown.h)
struct Leaf {
  int I;                                                        //!< Unique index for every leaf
  vect X;                                                       //!< Coordinate of leaf
  Leaf *NEXT;                                                   //!< Pointer to next leaf
  Leaf() : I(0), X(0), NEXT(NULL) {}                            //!< Constructor
  ~Leaf() {}                                                    //!< Destructor
//! Copy constructor
  Leaf(const Leaf &leaf) : I(0), X(0), NEXT(NULL) {
    I    = leaf.I;                                              // Copy I
    X    = leaf.X;                                              // Copy X
    NEXT = leaf.NEXT;                                           // Copy NEXT
  }
//! Overload assignment
  Leaf &operator=(const Leaf &leaf) {
    I    = leaf.I;                                              // Copy I
    X    = leaf.X;                                              // Copy X
    NEXT = leaf.NEXT;                                           // Copy NEXT
    return *this;
  }
};
typedef std::vector<Leaf>              Leafs;                   //!< Vector of leafs
typedef std::vector<Leaf>::iterator    L_iter;                  //!< Iterator for leaf vector

//! Structure of nodes (only used in fast/topdown.h)
struct Node {
  bool NOCHILD;                                                 //!< Flag for twig nodes
  int  LEVEL;                                                   //!< Level in the tree structure
  int  NLEAF;                                                   //!< Number of descendant leafs
  vec<8,int> CHILD;                                             //!< Index of child node
  vect X;                                                       //!< Coordinate at center
  Leaf *LEAF;                                                   //!< Pointer to first leaf
  Node() : NOCHILD(true), LEVEL(0), NLEAF(0), CHILD(-1), X(0), LEAF(NULL) {}//!< Constructor
  ~Node() {}                                                    //!< Destructor
//! Copy constructor
  Node(const Node &node) : NOCHILD(true), LEVEL(0), NLEAF(0), CHILD(-1), X(0), LEAF(NULL) {
    NOCHILD = node.NOCHILD;                                     // Copy NOCHILD
    LEVEL   = node.LEVEL;                                       // Copy LEVEL
    NLEAF   = node.NLEAF;                                       // Copy NLEAF
    CHILD   = node.CHILD;                                       // Copy CHILD
    X       = node.X;                                           // Copy X
    LEAF    = node.LEAF;                                        // Copy LEAF
  }
//! Overload assignment
  Node &operator=(const Node &node) {
    NOCHILD = node.NOCHILD;                                     // Copy NOCHILD
    LEVEL   = node.LEVEL;                                       // Copy LEVEL
    NLEAF   = node.NLEAF;                                       // Copy NLEAF
    CHILD   = node.CHILD;                                       // Copy CHILD
    X       = node.X;                                           // Copy X
    LEAF    = node.LEAF;                                        // Copy LEAF
    return *this;
  }
};
typedef std::vector<Node>              Nodes;                   //!< Vector of nodes
typedef std::vector<Node>::iterator    N_iter;                  //!< Iterator for node vector

//! Structure of source cells (stuff to send)
struct JCell {
  bigint ICELL;                                                 //!< Cell index
  Mset   M;                                                     //!< Multipole coefficients
};
typedef std::vector<JCell>             JCells;                  //!< Vector of source cells
typedef std::vector<JCell>::iterator   JC_iter;                 //!< Iterator for source cell vector

//! Structure of cells
struct Cell {
  bigint ICELL;                                                 //!< Cell index
  int    NCHILD;                                                //!< Number of child cells
  int    NCLEAF;                                                //!< Number of child leafs
  int    NDLEAF;                                                //!< Number of descendant leafs
  int    PARENT;                                                //!< Iterator offset of parent cell
  int    CHILD;                                                 //!< Iterator offset of child cells
  int    ILEAF;                                                 //!< Iterator offset of first leaf
  B_iter LEAF;                                                  //!< Iterator of first leaf
  vect   X;                                                     //!< Cell center
  real   R;                                                     //!< Cell radius
  real   RMAX;                                                  //!< Max cell radius
  real   RCRIT;                                                 //!< Critical cell radius
  Mset   M;                                                     //!< Multipole coefficients
  Lset   L;                                                     //!< Local coefficients
  Cell() : ICELL(0), NCHILD(0), NCLEAF(0), NDLEAF(0), PARENT(0), CHILD(0),
           LEAF(), X(0), R(0), RMAX(0), RCRIT(0), M(0), L(0) {} //!< Constructor
};
typedef std::vector<Cell>              Cells;                   //!< Vector of cells
typedef std::vector<Cell>::iterator    C_iter;                  //!< Iterator for cell vector
typedef std::queue<C_iter>             CellQueue;               //!< Queue of cell iterators
typedef std::stack<C_iter>             CellStack;               //!< Stack of cell iterators
typedef std::pair<C_iter,C_iter>       Pair;                    //!< Pair of interacting cells
typedef std::queue<Pair>               PairQueue;               //!< Queue of interacting cell pairs
typedef std::stack<Pair>               PairStack;               //!< Stack of interacting cell pairs
typedef std::list<C_iter>              List;                    //!< Interaction list
typedef std::list<C_iter>::iterator    LC_iter;                 //!< Iterator for interaction list
typedef std::vector<List>              Lists;                   //!< Vector of interaction lists
typedef std::map<C_iter,int>           Map;                     //!< Map of interaction lists
typedef std::map<C_iter,int>::iterator MC_iter;                 //!< Iterator for interation list map
typedef std::vector<Map>               Maps;                    //!< Vector of map of interaction lists

//! Structure for Ewald summation
struct Ewald {
  vect K;                                                       //!< 3-D wave number vector
  real REAL;                                                    //!< Real part of wave
  real IMAG;                                                    //!< Imaginary part of wave
  Ewald() : K(0), REAL(0), IMAG(0) {}                           //!< Constructor
};
typedef std::vector<Ewald>             Ewalds;                  //!< Vector of Ewald summation types
typedef std::vector<Ewald>::iterator   E_iter;                  //!< Iterator for Ewald summation types

#endif
