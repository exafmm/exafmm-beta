#ifndef types_h
#define types_h
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <vector>
#include "gettime.h"
#include "vec.h"

typedef long         bigint;                                    // Big integer type
typedef float        real;                                      // Real number type
typedef vec<3,real>  vect;                                      // 3-D vector type
typedef vec<10,real> coef;                                      // Multipole coefficient type
typedef std::vector<bigint>           Bigints;                  // Vector of big integer types
typedef std::vector<bigint>::iterator BI_iter;                  // Vector of big integer types

int  const NCRIT(100);                                          // Number of bodies per cell
real const THETA(0.5);                                          // Box opening criteria
real const EPS(0.000001);                                       // Single precision epsilon
real const EPS2(0.0001);                                        // Softening parameter

struct jbody {                                                  // Source properties of a body
  vect pos;                                                     // Position
  real scal;                                                    // Mass/charge
};
struct body : jbody {                                           // All properties of a body
  vect acc;                                                     // Acceleration
  real pot;                                                     // Potential
};
typedef std::vector<body>             Bodies;                   // Vector of bodies
typedef std::vector<body>::iterator   B_iter;                   // Iterator for body vector

struct jcell {                                                  // Source properties of a cell
  bigint I;                                                     // Cell index
  coef   M;                                                     // Multipole coefficients
};
struct cell : jcell {                                           // All properties of a cell
  typedef std::vector<cell>::iterator C_iter;                   // Iterator for cell vector
  int    NLEAF;                                                 // Number of leafs
  int    NCHILD;                                                // Number of child cells
  B_iter LEAF;                                                  // Pointer to first leaf
  C_iter PARENT;                                                // Pointer to parent cell
  C_iter CHILD[8];                                              // Pointer to child cells
  vect   X;                                                     // Cell center
  real   R;                                                     // Cell radius
  coef   L;                                                     // Local coefficients
};
typedef std::vector<cell>             Cells;                    // Vector of cells
typedef std::vector<cell>::iterator   C_iter;                   // Iterator for cell vector

struct pair {                                                   // Structure for pair of interacting cells
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator
  pair(C_iter ci, C_iter cj) : CI(ci), CJ(cj) {}                // Constructor
};
typedef std::stack<pair>              Pairs;                    // Stack of interacting cells

#endif
