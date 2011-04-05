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
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"
#ifndef KERNEL
int MPIRANK = 0;                                                // MPI rank (for debugging serial class in MPI run)
int MPISIZE = 1;                                                // MPI size (for debugging serial class in MPI run)
#else
extern int MPIRANK;
extern int MPISIZE;
#endif

typedef long                 bigint;                            // Big integer type
typedef float                real;                              // Real number type
typedef std::complex<double> complex;                           // Complex number type

const int  P       = 7;                                         // Order of expansions
const int  NCRIT   = 100;                                       // Number of bodies per cell
const real THETA   = 1/sqrtf(3);                                // Box opening criteria
const real CLET    = 3;                                         // LET opening critetia
const real EPS2    = 1e-4;                                      // Softening parameter
const int  IMAGES  = 0;                                         // Number of periodic image sublevels
const int  GPUS    = 3;                                         // Number of GPUs per node
const int  THREADS = 64;                                        // Number of threads per thread-block

#if Cartesian
const int  NTERM   = P*(P+1)*(P+2)/6;                           // Number of terms for cartesian expansion
#elif Spherical
const int  NTERM   = P*(P+1)/2;                                 // Number of terms for spherical harmonics
#endif
const int  NCOEF   = 3 * NTERM;                                 // 3-D vector of coefficients

typedef vec<3,real>                            vect;            // 3-D vector type
#if Cartesian
typedef vec<NCOEF,real>                        coef;            // Multipole coefficient type for Cartesian
#elif Spherical
typedef vec<NCOEF,complex>                     coef;            // Multipole coefficient type for spherical
#endif
typedef std::vector<bigint>                    Bigints;         // Vector of big integer types
typedef std::map<std::string,double>           Event;           // Map of event name to logged value
typedef std::map<std::string,double>::iterator E_iter;          // Iterator for event name map

struct JBody {                                                  // Source properties of a body (stuff to send)
  int         IBODY;                                            // Initial body numbering for sorting back
  int         IPROC;                                            // Initial process numbering for partitioning back
  bigint      ICELL;                                            // Cell index
  vect        X;                                                // Position
  vec<4,real> SRC;                                              // Source values
};
struct Body : JBody {                                           // All properties of a body
  vec<4,real> TRG;                                              // Target values
  bool operator<(const Body &rhs) const {                       // Overload operator for comparing body index
    return this->IBODY < rhs.IBODY;                             // Comparison function for body index
  }
};
typedef std::vector<Body>              Bodies;                  // Vector of bodies
typedef std::vector<Body>::iterator    B_iter;                  // Iterator for body vector
typedef std::vector<JBody>             JBodies;                 // Vector of source bodies
typedef std::vector<JBody>::iterator   JB_iter;                 // Iterator for source body vector

struct JCell {                                                  // Source properties of a cell (stuff to send)
  bigint ICELL;                                                 // Cell index
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
typedef std::vector<Map>               Maps;                    // Vector of map of interaction lists

#endif
