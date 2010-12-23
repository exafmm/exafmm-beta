#ifndef body_h
#define body_h
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "gettime.h"
#include "vec.h"

typedef unsigned long long u64;
typedef float real;
const int NCRIT = 6;
typedef vec<3,real> vect;
const real EPS = 0.01;

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
  unsigned const NMAX;
  body *P;
public:
  bodies(unsigned const n, unsigned const nmax) : N(n),NMAX(nmax) {
    P = new body [NMAX];
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
  body &operator[](unsigned const i) const {
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
  friend std::ostream &operator<<(std::ostream &s, bodies const &P) {
    s<<P.I;
    return s;
  }
};
#endif
