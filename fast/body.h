#ifndef body_h
#define body_h
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <sys/time.h>
#include <vector>
#include <vec.h>

typedef long                 bigint;                            // Big integer type
typedef float                real;                              // Real number type on CPU
typedef float                gpureal;                           // Real number type on GPU
typedef std::complex<double> complex;                           // Complex number type

const int P = 3;
const int NCOEF = (P+1)*(P+2)*(P+3)/6;
const int NCRIT = 8;
const float THETA = 0.6;
const float EPS = 0.01;
const float EQ = EPS * EPS;
typedef vec<3 ,real> vect;
typedef vec<7 ,real> Mset;
typedef vec<20,real> Lset;
const real zero = 0.;

double get_time(void) {
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

struct JBody {                                                  // Source properties of a body (stuff to send)
  int         IBODY;                                            // Initial body numbering for sorting back
//  int         IPROC;                                            // Initial process numbering for partitioning back
//  bigint      ICELL;                                            // Cell index
  vect        X;                                                // Position
  vec<1,real> SRC;                                              // Source values
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

struct Cell {
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
#endif
