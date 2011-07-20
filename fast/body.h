#ifndef body_h
#define body_h
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
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

struct Cell {
  unsigned LEVEL;
  unsigned NCLEAF;
  unsigned NCCELL;
  unsigned NDLEAF;
  B_iter   FCLEAF;
  Cell    *FCCELL;
  Cell    *PARENT;
  vect X;
  real RCRIT;
  Mset M;
  Lset L;
};

template<typename A, typename B>
struct Pair {
  A first;                                       // first object
  B second;                                      // second object
  Pair() : first(0), second(0) {}                // constructor
  void set(A a, B b) {                           // set pair
    first=a;  second=b;
  }
};

template<typename A, typename B>
class Stack {
private:
  typedef Pair<A,B> pair;                        // type of stack objects
  pair    *P0, *P;                               // first & current element
public:
  explicit
  Stack (unsigned const&M)                       // constructor
    : P0 ( new pair [M] ),                       //   allocate memory
      P  ( P0 - 1 ) {}                           //   set pter to activ
  ~Stack () { delete[] P0; }                     // destructor: deallocate
  bool empty() const { return P<P0; }            // if stack empty?
  pair pop() { return *(P--); }                  // give last:   pop
  void push(A a, B b) {                          // add element: push
    ++P;
    P->set(a,b);
  }
};

#endif
