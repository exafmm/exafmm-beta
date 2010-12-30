#ifndef types_h
#define types_h
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "gettime.h"
#include "vec.h"

typedef long long    bigint;
typedef float        real;
typedef vec<3,real>  vect;
typedef vec<20,real> coef;

const int  NCRIT = 100;
const real   EPS = 0.01;

struct body
{
  vect pos;
  real scal;
  vect acc;
  real pot;
};
typedef std::vector<body>           Bodies;
typedef std::vector<body>::iterator B_iter;

#define C_iter std::vector<cell>::iterator
struct cell {
  int    NLEAF;
  int    NCHILD;
  bigint I;
  B_iter LEAF;
  C_iter PARENT;
  C_iter CHILD[8];
  coef   *M,*L;
};
#undef C_iter
typedef std::vector<cell>           Cells;
typedef std::vector<cell>::iterator C_iter;

#endif
