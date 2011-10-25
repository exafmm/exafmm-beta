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

const int  P       = 4;                                         // Order of expansions
const int  NCRIT   = 8;                                         // Number of bodies per cell
const int  MAXBODY = 200000;                                    // Maximum number of bodies per GPU kernel
const int  MAXCELL = 10000000;                                  // Maximum number of bodies/coefs in cell per GPU kernel
const real CLET    = 2;                                         // LET opening critetia
const real EPS2    = 1e-6;                                      // Softening parameter
const int  GPUS    = 3;                                         // Number of GPUs per node
const int  THREADS = 64;                                        // Number of threads per thread-block

const int MCOEF = P*(P+1)*(P+2)/6-3;
const int LCOEF = (P+1)*(P+2)*(P+3)/6;
const int NCOEF = P*(P+1)*(P+2)/6;
typedef vec<3 ,real> vect;
#if CART
typedef vec<NCOEF,real> Mset;
typedef vec<NCOEF,real> Lset;
#elif SPHE
typedef vec<MCOEF,complex> Mset;
typedef vec<NCOEF,complex> Lset;
#else
typedef vec<MCOEF,real> Mset;
typedef vec<LCOEF,real> Lset;
#endif

typedef std::vector<bigint>                    Bigints;         // Vector of big integer types
typedef std::map<std::string,double>           Event;           // Map of event name to logged value
typedef std::map<std::string,double>::iterator E_iter;          // Iterator for event name map

struct JBody {                                                  // Source properties of a body (stuff to send)
  int         IBODY;                                            // Initial body numbering for sorting back
  int         IPROC;                                            // Initial process numbering for partitioning back
  bigint      ICELL;                                            // Cell index
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
  int  LEVEL;
  int  NLEAF;
  int  CHILD[8];
  vect X;
  Leaf *LEAF;
};
typedef std::vector<Node>           Nodes;
typedef std::vector<Node>::iterator N_iter;

struct Cell {
  unsigned ICELL;
  unsigned NCHILD;
  unsigned NCLEAF;
  unsigned NDLEAF;
  unsigned PARENT;
  unsigned CHILD;
  B_iter   LEAF;
  vect X;
  real R;
  real RCRIT;
  Mset M;
  Lset L;
};
typedef std::vector<Cell>           Cells;
typedef std::vector<Cell>::iterator C_iter;
#endif
