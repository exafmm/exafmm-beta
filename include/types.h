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
#include "quark.h"
#include <queue>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"                                                //!< My vector type with operator overloading

typedef unsigned           bigint;                              //!< Big integer type
typedef float              real;                                //!< Real number type on CPU
typedef float              gpureal;                             //!< Real number type on GPU
typedef std::complex<real> complex;                             //!< Complex number type
typedef vec<3,real>        vect;                                //!< 3-D vector type


#ifndef KERNEL
int MPIRANK = 0;                                                //!< MPI comm rank
int MPISIZE = 1;                                                //!< MPI comm size
int DEVICE  = 0;                                                //!< GPU device ID
int IMAGES  = 0;                                                //!< Number of periodic image sublevels
real THETA  = .5;                                               //!< Box opening criteria
vect Xperiodic = 0;                                             //!< Coordinate offset of periodic image
#else
extern int MPIRANK;                                             //!< MPI comm rank
extern int MPISIZE;                                             //!< MPI comm size
extern int DEVICE;                                              //!< GPU device ID
extern int IMAGES;                                              //!< Number of periodic image sublevels
extern real THETA;                                              //!< Box opening criteria
extern vect Xperiodic;                                          //!< Coordinate offset of periodic image
#endif

const int  P        = 10;                                       //!< Order of expansions
const int  NCRIT    = 100;                                      //!< Number of bodies per cell
const int  MAXBODY  = 200000;                                   //!< Maximum number of bodies per GPU kernel
const int  MAXCELL  = 10000000;                                 //!< Maximum number of bodies/coefs in cell per GPU kernel
const real CLET     = 2;                                        //!< LET opening critetia
const real EPS      = 1e-6;                                     //!< Single precision epsilon
const real EPS2     = 0;                                        //!< Softening parameter (squared)
const int  GPUS     = 3;                                        //!< Number of GPUs per node
const int  THREADS  = 64;                                       //!< Number of threads per thread-block
const int  PTHREADS = 4;                                        //!< Number of pthreads in quark

const int MTERM = P*(P+1)*(P+2)/6-3;                            //!< Number of Cartesian mutlipole terms
const int LTERM = (P+1)*(P+2)*(P+3)/6;                          //!< Number of Cartesian local terms
const int NTERM = P*(P+1)/2;                                    //!< Number of Spherical multipole/local terms

#if Cartesian
typedef vec<MTERM,real>                        Mset;            //!< Multipole coefficient type for Cartesian
typedef vec<LTERM,real>                        Lset;            //!< Local coefficient type for Cartesian
#elif Spherical
typedef vec<NTERM,complex>                     Mset;            //!< Multipole coefficient type for spherical
typedef vec<NTERM,complex>                     Lset;            //!< Local coefficient type for spherical
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
  CoulombVdW                                                    //!< Coulomb + Van der Walls force
};

//! Structure of source bodies (stuff to send)
struct JBody {
  int         IBODY;                                            //!< Initial body numbering for sorting back
  int         IPROC;                                            //!< Initial process numbering for partitioning back
  unsigned    ICELL;                                            //!< Cell index
  vect        X;                                                //!< Position
  real        SRC;                                              //!< Scalar source values
};
//! Structure of bodies
struct Body : JBody {
  vec<4,real> TRG;                                              //!< Scalar+vector target values
  bool operator<(const Body &rhs) const {                       //!< Overload operator for comparing body index
    return this->IBODY < rhs.IBODY;                             //!< Comparison function for body index
  }
};
typedef std::vector<Body>              Bodies;                  //!< Vector of bodies
typedef std::vector<Body>::iterator    B_iter;                  //!< Iterator for body vector
typedef std::vector<JBody>             JBodies;                 //!< Vector of source bodies
typedef std::vector<JBody>::iterator   JB_iter;                 //!< Iterator for source body vector

struct Leaf {
  int I;
  vect X;
  Leaf *NEXT;
};
typedef std::vector<Leaf>           Leafs;                      // Vector of leafs
typedef std::vector<Leaf>::iterator L_iter;                     // Iterator for leaf vector

struct Node {
  bool NOCHILD;
  int  LEVEL;
  int  NLEAF;
  int  CHILD[8];
  vect X;
  Leaf *LEAF;
};
typedef std::vector<Node>           Nodes;
typedef std::vector<Node>::iterator N_iter;

//! Structure of source cells (stuff to send)
struct JCell {
  unsigned ICELL;                                               //!< Cell index
  Mset   M;                                                     //!< Multipole coefficients
};
//! Structure of cells
struct Cell {
  unsigned ICELL;
  int      NCHILD;                                              //!< Number of child cells
  int      NCLEAF;                                              //!< Number of child leafs
  int      NDLEAF;                                              //!< Number of descendant leafs
  int      PARENT;                                              //!< Iterator offset of parent cell
  int      CHILD;                                               //!< Iterator offset of child cells
  B_iter   LEAF;                                                //!< Iterator of first leaf
  vect     X;                                                   //!< Cell center
  real     R;                                                   //!< Cell radius
  real     RCRIT;                                               //!< Critical cell radius
  Mset     M;                                                   //!< Multipole coefficients
  Lset     L;                                                   //!< Local coefficients
};
typedef std::vector<Cell>              Cells;                   //!< Vector of cells
typedef std::vector<Cell>::iterator    C_iter;                  //!< Iterator for cell vector
typedef std::vector<JCell>             JCells;                  //!< Vector of source cells
typedef std::vector<JCell>::iterator   JC_iter;                 //!< Iterator for source cell vector

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

#endif
