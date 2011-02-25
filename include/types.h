#ifndef types_h
#define types_h
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <stack>
#include <utility>
#include <vector>
#include "vec.h"
#ifndef KERNEL
#include "gettime.h"
int MPIRANK = 0;                                                // MPI rank (for debugging serial class in MPI run)
int MPISIZE = 1;                                                // MPI size (for debugging serial class in MPI run)
#else
extern int MPIRANK;
extern int MPISIZE;
#endif
typedef long                 bigint;                            // Big integer type
typedef float                real;                              // Real number type
typedef std::complex<double> complex;                           // Complex number type

const int  P     = 3;                                           // Order of expansions
//const int  NCOEF = P*(P+1)*(P+2)/6;                             // Number of coefficients for Taylor expansion
const int  NCOEF   = P*(P+1)/2;                                 // Number of coefficients for spherical harmonics
const int  NCRIT   = 100;                                       // Number of bodies per cell
const real THETA   = 0.5;                                       // Box opening criteria
const real EPS2    = 1e-4;                                      // Softening parameter
const int  GPUS    = 4;                                         // Number of GPUs per node
const int  THREADS = 256;                                       // Number of threads per thread-block

typedef vec<3,real>                    vect;                    // 3-D vector type
//typedef vec<NCOEF,real>                coef;                    // Multipole coefficient type for Taylor expansion
typedef vec<NCOEF,complex>             coef;                    // Multipole coefficient type for spherical harmonics
typedef std::vector<bigint>            Bigints;                 // Vector of big integer types

struct JBody {                                                  // Source properties of a body (stuff to send)
  bigint I;                                                     // Cell index
  vect   pos;                                                   // Position
  real   scal;                                                  // Mass/charge
};
struct Body : JBody {                                           // All properties of a body
  vect acc;                                                     // Acceleration
  real pot;                                                     // Potential
};
typedef std::vector<Body>              Bodies;                  // Vector of bodies
typedef std::vector<Body>::iterator    B_iter;                  // Iterator for body vector
typedef std::vector<JBody>             JBodies;                 // Vector of source bodies
typedef std::vector<JBody>::iterator   JB_iter;                 // Iterator for source body vector

struct JCell {                                                  // Source properties of a cell (stuff to send)
  bigint I;                                                     // Cell index
  coef   M;                                                     // Multipole coefficients
};
struct Cell : JCell {                                           // All properties of a cell
  int    NCHILD;                                                // Number of child cells
  int    NLEAF;                                                 // Number of leafs
  int    PARENT;                                                // Iterator offset of parent cell
  int    CHILD[8];                                              // Iterator offset of child cells
  B_iter LEAF;                                                  // Iterator of first leaf
  vect   X;                                                     // Cell center
  real   R;                                                     // Cell radius
  coef   L;                                                     // Local coefficients
};
typedef std::vector<Cell>              Cells;                   // Vector of cells
typedef std::vector<Cell>::iterator    C_iter;                  // Iterator for cell vector
typedef std::vector<JCell>             JCells;                  // Vector of source cells
typedef std::vector<JCell>::iterator   JC_iter;                 // Iterator for source cell vector

typedef std::pair<C_iter,C_iter>       Pair;                    // Pair of interacting cells
typedef std::stack<Pair>               Pairs;                   // Stack of interacting cell pairs
typedef std::list<C_iter>              List;                    // Interaction list
typedef std::list<C_iter>::iterator    L_iter;                  // Iterator for interaction list vector
typedef std::vector<List>              Lists;                   // Vector of interaction lists
typedef std::map<C_iter,int>           Map;                     // Map of interaction lists
typedef std::map<C_iter,int>::iterator M_iter;                  // Iterator for interation list map

#endif
