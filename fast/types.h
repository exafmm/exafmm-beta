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
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <omp.h>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"

typedef long                 bigint;                            // Big integer type
typedef float                real;                              // Real number type on CPU
typedef float                gpureal;                           // Real number type on GPU
typedef std::complex<double> complex;                           // Complex number type

#ifndef KERNEL
int MPIRANK = 0;                                                // MPI comm rank
int MPISIZE = 1;                                                // MPI comm size
int DEVICE  = 0;                                                // GPU device ID
int IMAGES;                                                     // Number of periodic image sublevels
real THETA;                                                     // Box opening criteria
#else
extern int MPIRANK;                                             // MPI comm rank
extern int MPISIZE;                                             // MPI comm size
extern int DEVICE;                                              // GPU device ID
extern int IMAGES;                                              // Number of periodic image sublevels
extern real THETA;                                              // Box opening criteria
#endif

const int  P       = 3;                                         // Order of expansions
const int  NCRIT   = 8;                                         // Number of bodies per cell
const int  MAXBODY = 200000;                                    // Maximum number of bodies per GPU kernel
const int  MAXCELL = 10000000;                                  // Maximum number of bodies/coefs in cell per GPU kernel
const real CLET    = 2;                                         // LET opening critetia
const real EPS2    = 1e-6;                                      // Softening parameter
const int  GPUS    = 3;                                         // Number of GPUs per node
const int  THREADS = 64;                                        // Number of threads per thread-block

const int MTERM = P*(P+1)*(P+2)/6-3;
const int LTERM = (P+1)*(P+2)*(P+3)/6;
const int NTERM = P*(P+1)/2;

typedef vec<3 ,real> vect;
#if Cartesian
typedef vec<MTERM,real> Mset;
typedef vec<LTERM,real> Lset;
#elif Spherical
typedef vec<3*NTERM,complex> Mset;
typedef vec<3*NTERM,complex> Lset;
#endif

typedef std::vector<bigint>                    Bigints;         // Vector of big integer types
typedef std::map<std::string,double>           Event;           // Map of event name to logged value
typedef std::map<std::string,double>::iterator E_iter;          // Iterator for event name map

struct JBody {                                                  // Source properties of a body (stuff to send)
  int         IBODY;                                            // Initial body numbering for sorting back
  int         IPROC;                                            // Initial process numbering for partitioning back
  unsigned    ICELL;                                            // Cell index
  vect        X;                                                // Position
  vec<1,real> SRC;                                              // Source values
};
struct Body : JBody {                                           // All properties of a body
  vec<4,real> TRG;                                              // Target values
  bool operator<(const Body &rhs) const {                       // Overload operator for comparing body index
    return this->IBODY < rhs.IBODY;                             // Comparison function for body index
  }
};
typedef std::vector<Body>           Bodies;                     // Vector of bodies
typedef std::vector<Body>::iterator B_iter;                     // Iterator for body vector

struct Leaf {
  int I;
  vect X;
  Leaf *NEXT;
};
typedef std::vector<Leaf>           Leafs;                      // Vector of leafs
typedef std::vector<Leaf>::iterator L_iter;                     // Iterator for leaf vector

struct Node {
  bool NOCHILD;
  int LEVEL;
  int NLEAF;
  int CHILD[8];
  vect X;
  Leaf *LEAF;
};
typedef std::vector<Node>           Nodes;
typedef std::vector<Node>::iterator N_iter;

struct Cell {
  unsigned ICELL;
  int NCHILD;
  int NCLEAF;
  int NDLEAF;
  int PARENT;
  int CHILD;
  B_iter LEAF;
  vect X;
  real R;
  real RCRIT;
  Mset M;
  Lset L;
};
typedef std::vector<Cell>           Cells;
typedef std::vector<Cell>::iterator C_iter;
#endif
