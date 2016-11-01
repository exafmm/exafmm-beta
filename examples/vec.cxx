#define EXAFMM_VEC_VERBOSE 1
#include <assert.h>
#include <cmath>
#include <iostream>
#include "vec.h"

using namespace exafmm;

template<int S, typename T>
struct SIMD {
  static const T EPS = 1e-5;
  static inline vec<S,T> sequence() {
    vec<S,T> v(0);
    return v;
  }
  static inline void print() {
    std::cout << "Is not a SIMD type" << std::endl;
  }
};

#if defined __AVX512F__ || defined __MIC__
template<>
struct SIMD<16,float> {
  static const float EPS = 1e-5;
  static inline vec<16,float> sequence() {
    vec<16,float> v(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    return v;
  }
  static inline void print() {
    std::cout << "AVX512 float  : pass" << std::endl;
  }
};

template<>
struct SIMD<8,double> {
  static const double EPS = 1e-13;
  static inline vec<8,double> sequence() {
    vec<8,double> v(0,1,2,3,4,5,6,7);
    return v;
  }
  static inline void print() {
    std::cout << "AVX512 double : pass" << std::endl;
  }
};
#endif

#ifdef __AVX__
template<>
struct SIMD<8,float> {
  static const float EPS = 1e-5;
  static inline vec<8,float> sequence() {
    vec<8,float> v(0,1,2,3,4,5,6,7);
    return v;
  }
  static inline void print() {
    std::cout << "AVX float     : pass" << std::endl;
  }
};

template<>
struct SIMD<4,double> {
  static const double EPS = 1e-13;
  static inline vec<4,double> sequence() {
    vec<4,double> v(0,1,2,3);
    return v;
  }
  static inline void print() {
    std::cout << "AVX double    : pass" << std::endl;
  }
};
#endif

#ifdef __SSE__
template<>
struct SIMD<4,float> {
  static const float EPS = 1e-5;
  static inline vec<4,float> sequence() {
    vec<4,float> v(0,1,2,3);
    return v;
  }
  static inline void print() {
    std::cout << "SSE float     : pass" << std::endl;
  }
};

template<>
struct SIMD<2,double> {
  static const double EPS = 1e-13;
  static inline vec<2,double> sequence() {
    vec<2,double> v(0,1);
    return v;
  }
  static inline void print() {
    std::cout << "SSE double    : pass" << std::endl;
  }
};
#endif

template<int S, typename T>
void test() {
  T EPS = SIMD<S,T>::EPS;
  vec<S,T> a(1);
  vec<S,T> b(a+1);
  vec<S,T> c = SIMD<S,T>::sequence();
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-1)<EPS);
    assert(std::abs(b[i]-2)<EPS);
    assert(std::abs(c[i]-i)<EPS);
  }
  a = 2;
  b = c;
  c += 1;
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-2)<EPS);
    assert(std::abs(b[i]-i)<EPS);
    assert(std::abs(c[i]-i-1)<EPS);
  }
  a -= 1;
  b *= 2;
  c /= 10;
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-1)<EPS);
    assert(std::abs(b[i]-2*i)<EPS);
    assert(std::abs(c[i]-(i+1.)/10)<EPS);
  }
  b &= c > (a / 2);
  a = b + c;
  c = a - b;
  for (int i=0; i<S; i++) {
    if (i < 5) {
      assert(std::abs(a[i]-(i+1.)/10)<EPS);
      assert(std::abs(b[i])<EPS);
    } else {
      assert(std::abs(a[i]-(2.1*i+0.1)<EPS));
      assert(std::abs(b[i]-(2*i)<EPS));
    }
    assert(std::abs(c[i]-(i+1.)/10)<EPS);
  }
  a = c * c;
  b = a / c;
  c = -b;
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-(i+1.)*(i+1.)/100)<EPS);
    assert(std::abs(b[i]-(i+1.)/10)<EPS);
    assert(std::abs(c[i]+(i+1.)/10)<EPS);
  }
  a = rsqrt(a);
  b = min(a,b);
  c = max(a,c);
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-10/(i+1.))<EPS);
    assert(std::abs(b[i]-std::min((i+1.)/10,double(a[i])))<EPS);
    assert(std::abs(c[i]-10/(i+1.))<EPS);
  }
  a = a[0]/5;
  b = sum(a);
  c = norm(a);
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-2)<EPS);
    assert(std::abs(b[i]-2*S)<EPS);
    assert(std::abs(c[i]-4*S)<EPS);
  }
  a = M_PI/2;
  b = 0;
  a = sin(a);
  b = cos(b);
  c = exp(a);
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-1)<EPS);
    assert(std::abs(b[i]-1)<EPS);
    assert(std::abs(c[i]-M_E)<EPS);
  }
  c = M_PI/2;
  sincos(a,b,c);
  for (int i=0; i<S; i++) {
    assert(std::abs(a[i]-1)<EPS);
    assert(std::abs(b[i])<EPS);
  }
  SIMD<S,T>::print();
}

int main(int , char ** ) {
#if defined __AVX512F__ || defined __MIC__
  test<16,float>();
  test<8,double>();
#endif
#ifdef __AVX__
  test<8,float>();
  test<4,double>();
#endif
#ifdef __SSE__
  test<4,float>();
  test<2,double>();
#endif
  return 0;
};
