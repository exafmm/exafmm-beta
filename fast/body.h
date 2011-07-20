#ifndef body_h
#define body_h
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <vec.h>

typedef float real;
const int P = 3;
const int NCOEF = (P+1)*(P+2)*(P+3)/6;
const int NCRIT = 8;
const float THETA = 0.6;
const float EPS = 0.01;
const float EQ = EPS * EPS;
typedef vec<3 ,real> vect;                      ///< a vector of 3 reals
typedef vec<7 ,real> Mset;                      ///< set of multipoles
typedef vec<20,real> Lset;                      ///< set of coefficients
const real zero = 0.;                          ///< real: zero

double get_time(void) {
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

struct body
{
  vect pos;
  real scal;
  vect acc;
  real pot;
};

class bodies
{
  unsigned I;
  unsigned const N;
  body *P;
public:
  bodies(unsigned n) : N(n) {
    P = new body [N];
  }
  ~bodies() {
    delete[] P;
  }
  vect &pos() const {
    return P[I].pos;
  }
  real &scal() const {
    return P[I].scal;
  }
  vect &acc() const {
    return P[I].acc;
  }
  real &pot() const {
    return P[I].pot;
  }
  vect &pos(unsigned i) const {
    return P[i].pos;
  }
  real &scal(unsigned i) const {
    return P[i].scal;
  }
  vect &acc(unsigned i) const {
    return P[i].acc;
  }
  real &pot(unsigned i) const {
    return P[i].pot;
  }
  unsigned begin() {
    I = 0;
    return I;
  }
  unsigned end() const {
    return N;
  }
  unsigned size() const {
    return N;
  }
  unsigned index() const {
    return I;
  }
  body const &operator[](unsigned const i) const {
    return P[i];
  }
  unsigned const &operator=(unsigned i) {
    return I = i;
  }
  bodies const &operator++() {
    ++I;
    return *this;
  }
  bool operator!=(unsigned i) const {
    return I != i;
  }
  operator unsigned () {return I;}
  friend std::ostream &operator<<(std::ostream &s, bodies const &B) {
    s<<B.I;
    return s;
  }
};

struct Leaf {
  unsigned I;
  vect X;
  real Q;
  vec<4,real> TRG;
};

struct Cell {
  unsigned LEVEL;
  unsigned NCLEAF;
  unsigned NCCELL;
  unsigned NDLEAF;
  Leaf    *FCLEAF;
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
