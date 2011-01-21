#ifndef types_h
#define types_h
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <vector>
#include "gettime.h"
#include "vec.h"

typedef long long    bigint;                                    // Big integer type
typedef float        real;                                      // Real number type
typedef vec<3,real>  vect;                                      // 3-D vector type
typedef vec<10,real> coef;                                      // Multipole coefficient type

const int  NCRIT = 100;                                         // Number of bodies per cell
const real THETA = 0.5;                                         // Box opening criteria
const real  EPS2 = 0.0001;                                      // Softening parameter

struct body {                                                   // Structure for body
  vect pos;                                                     // Position
  real scal;                                                    // Mass/Charge
  vect acc;                                                     // Acceleration
  real pot;                                                     // Potential
};
typedef std::vector<body>           Bodies;                     // Vector of bodies
typedef std::vector<body>::iterator B_iter;                     // Iterator for body vector

#define C_iter std::vector<cell>::iterator                      // Temporary macro for cell iterator
struct cell {                                                   // Structure for cell
  int    NLEAF;                                                 // Number of leafs
  int    NCHILD;                                                // Number of child cells
  bigint I;                                                     // Morton index
  B_iter LEAF;                                                  // Pointer to first leaf
  C_iter PARENT;                                                // Pointer to parent cell
  C_iter CHILD[8];                                              // Pointer to child cells
  vect   X;                                                     // Cell center
  real   R;                                                     // Cell radius
  coef   M;                                                     // Pointer to multipole coefficients
  coef   L;                                                     // Pointer to local coefficients
};
#undef C_iter
typedef std::vector<cell>           Cells;                      // Vector of cells
typedef std::vector<cell>::iterator C_iter;                     // Iterator for cell vector

struct pair {                                                   // Structure for pair of interacting cells
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator
  pair(C_iter ci, C_iter cj) : CI(ci), CJ(cj) {}                // Constructor
};
typedef std::stack<pair>            Stack;                      // Stack of interacting cells

#endif
