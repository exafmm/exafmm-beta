#ifndef vec_h
#define vec_h
#include <ostream>
#include "../vectorclass/vectormath_trig.h"
#ifndef EXAFMM_VEC_NEWTON
#define EXAFMM_VEC_NEWTON 1
#endif

namespace exafmm {
#ifndef __CUDACC__
  template<int N, typename T>
  class vec {
  private:
    T data[N];
  public:
    vec(){}                                                     // Default constructor
    vec(const T & v) {                                          // Copy constructor (scalar)
      for (int i=0; i<N; i++) data[i] = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      for (int i=0; i<N; i++) data[i] = v[i];
    }
    ~vec(){}                                                    // Destructor
    const vec &operator=(const T v) {                           // Scalar assignment
      for (int i=0; i<N; i++) data[i] = v;
      return *this;
    }
    const vec &operator+=(const T v) {                          // Scalar compound assignment (add)
      for (int i=0; i<N; i++) data[i] += v;
      return *this;
    }
    const vec &operator-=(const T v) {                          // Scalar compound assignment (subtract)
      for (int i=0; i<N; i++) data[i] -= v;
      return *this;
    }
    const vec &operator*=(const T v) {                          // Scalar compound assignment (multiply)
      for (int i=0; i<N; i++) data[i] *= v;
      return *this;
    }
    const vec &operator/=(const T v) {                          // Scalar compound assignment (divide)
      for (int i=0; i<N; i++) data[i] /= v;
      return *this;
    }
    const vec &operator>=(const T v) {                          // Scalar compound assignment (greater than)
      for (int i=0; i<N; i++) data[i] >= v;
      return *this;
    }
    const vec &operator<=(const T v) {                          // Scalar compound assignment (less than)
      for (int i=0; i<N; i++) data[i] <= v;
      return *this;
    }
    const vec &operator&=(const T v) {                          // Scalar compound assignment (bitwise and)
      for (int i=0; i<N; i++) data[i] &= v;
      return *this;
    }
    const vec &operator|=(const T v) {                          // Scalar compound assignment (bitwise or)
      for (int i=0; i<N; i++) data[i] |= v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      for (int i=0; i<N; i++) data[i] = v[i];
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      for (int i=0; i<N; i++) data[i] += v[i];
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      for (int i=0; i<N; i++) data[i] -= v[i];
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      for (int i=0; i<N; i++) data[i] *= v[i];
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      for (int i=0; i<N; i++) data[i] /= v[i];
      return *this;
    }
    const vec &operator>=(const vec & v) {                      // Vector compound assignment (greater than)
      for (int i=0; i<N; i++) data[i] >= v[i];
      return *this;
    }
    const vec &operator<=(const vec & v) {                      // Vector compound assignment (less than)
      for (int i=0; i<N; i++) data[i] <= v[i];
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      for (int i=0; i<N; i++) data[i] &= v[i];
      return *this;
    }
    const vec &operator|=(const vec & v) {                      // Vector compound assignment (bitwise or)
      for (int i=0; i<N; i++) data[i] |= v[i];
      return *this;
    }
    vec operator+(const T v) const {                            // Scalar arithmetic (add)
      return vec(*this) += v;
    }
    vec operator-(const T v) const {                            // Scalar arithmetic (subtract)
      return vec(*this) -= v;
    }
    vec operator*(const T v) const {                            // Scalar arithmetic (multiply)
      return vec(*this) *= v;
    }
    vec operator/(const T v) const {                            // Scalar arithmetic (divide)
      return vec(*this) /= v;
    }
    vec operator>(const T v) const {                            // Scalar arithmetic (greater than)
      return vec(*this) >= v;
    }
    vec operator<(const T v) const {                            // Scalar arithmetic (less than)
      return vec(*this) <= v;
    }
    vec operator&(const T v) const {                            // Scalar arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    vec operator|(const T v) const {                            // Scalar arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return vec(*this) += v;
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return vec(*this) -= v;
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return vec(*this) *= v;
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return vec(*this) /= v;
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return vec(*this) >= v;
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return vec(*this) <= v;
    }
    vec operator&(const vec & v) const {                        // Vector arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    vec operator|(const vec & v) const {                        // Vector arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      vec temp;
      for (int i=0; i<N; i++) temp[i] = -data[i];
      return temp;
    }
    T &operator[](int i) {                                      // Indexing (lvalue)
      return data[i];
    }
    const T &operator[](int i) const {                          // Indexing (rvalue)
      return data[i];
    }
    operator       T* ()       {return data;}                   // Type-casting (lvalue)
    operator const T* () const {return data;}                   // Type-casting (rvalue)
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<N; i++) s << v[i] << ' ';
      return s;
    }
    friend T sum(const vec & v) {                               // Sum vector
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i];
      return temp;
    }
    friend T norm(const vec & v) {                              // L2 norm squared
      T temp = 0;
      for (int i=0; i<N; i++) temp += v[i] * v[i];
      return temp;
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
      return temp;
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
      return temp;
    }
    friend T min(const vec & v) {                               // Reduce minimum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp < v[i] ? temp : v[i];
      return temp;
    }
    friend T max(const vec & v) {                               // Reduce maximum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp > v[i] ? temp : v[i];
      return temp;
    }
    friend vec sin(const vec & v) {                             // Sine function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = sin(v[i]);
      return temp;
    }
    friend vec cos(const vec & v) {                             // Cosine function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = cos(v[i]);
      return temp;
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      for (int i=0; i<N; i++) {
	s[i] = sin(v[i]);
	c[i] = cos(v[i]);
      }
    }
    friend vec exp(const vec & v) {                             // Exponential function
      vec temp;
      for (int i=0; i<N; i++) temp[i] = exp(v[i]);
      return temp;
    }
    friend int wrap(vec & v, const vec & w) {                   // Wrap around periodic boundary
      int iw = 0;
      for (int i=0; i<N; i++) {
	if(v[i] < -w[i] / 2) {
	  v[i] += w[i];
	  iw |= 1 << i;
	}
	if(v[i] >  w[i] / 2) {
	  v[i] -= w[i];
	  iw |= 1 << i;
	}
      }
      return iw;
    }
    friend void unwrap(vec & v, const vec & w, const int & iw) {// Undo wrap around periodic boundary
      for (int i=0; i<N; i++) {
	if((iw >> i) & 1) v[i] += (v[i] > 0 ? -w[i] : w[i]);
      }
    }
  };
#else
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for CUDA")
#endif
#include "unroll.h"
  template<int N, typename T>
  class vec {
  private:
    T data[N];
  public:
    __host__ __device__ __forceinline__
    vec(){}                                                     // Default constructor
    __host__ __device__ __forceinline__
    vec(const T & v) {                                          // Copy constructor (scalar)
      Unroll<Ops::Assign<T>,T,N>::loop(data,v);
    }
    __host__ __device__ __forceinline__
    vec(const float4 & v) {                                     // Copy constructor (float4)
      data[0] = v.x;
      data[1] = v.y;
      data[2] = v.z;
      data[3] = v.w;
    }
    __host__ __device__ __forceinline__
    vec(const vec & v) {                                        // Copy constructor (vector)
      Unroll<Ops::Assign<T>,T,N>::loop(data,v);
    }
    __host__ __device__ __forceinline__
    vec(const float & x, const float & y, const float & z, const float & w) {// Copy constructor (4 floats)
      data[0] = x;
      data[1] = y;
      data[2] = z;
      data[3] = w;
    }
    __host__ __device__ __forceinline__
    vec(const float & x, const float & y, const float & z) {    // Copy constructor (3 floats)
      data[0] = x;
      data[1] = y;
      data[2] = z;
    }
    __host__ __device__ __forceinline__
    ~vec(){}                                                    // Destructor
    __host__ __device__ __forceinline__
    const vec &operator=(const T & v) {                         // Scalar assignment
      Unroll<Ops::Assign<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator+=(const T & v) {                        // Scalar compound assignment (add)
      Unroll<Ops::Add<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator-=(const T & v) {                        // Scalar compound assignment (subtract)
      Unroll<Ops::Sub<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator*=(const T & v) {                        // Scalar compound assignment (multiply)
      Unroll<Ops::Mul<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator/=(const T & v) {                        // Scalar compound assignment (divide)
      Unroll<Ops::Div<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator>=(const T & v) {                        // Scalar compound assignment (greater than)
      Unroll<Ops::Gt<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator<=(const T & v) {                        // Scalar compound assignment (less than)
      Unroll<Ops::Lt<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator&=(const T & v) {                        // Scalar compound assignment (bitwise and)
      Unroll<Ops::And<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator|=(const T & v) {                        // Scalar compound assignment (bitwise or)
      Unroll<Ops::Or<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator=(const vec & v) {                       // Vector assignment
      Unroll<Ops::Assign<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      Unroll<Ops::Add<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      Unroll<Ops::Sub<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      Unroll<Ops::Mul<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      Unroll<Ops::Div<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator>=(const vec & v) {                      // Vector compound assignment (greater than)
      Unroll<Ops::Gt<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator<=(const vec & v) {                      // Vector compound assignment (less than)
      Unroll<Ops::Lt<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      Unroll<Ops::And<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    const vec &operator|=(const vec & v) {                      // Vector compound assignment (bitwise or)
      Unroll<Ops::Or<T>,T,N>::loop(data,v);
      return *this;
    }
    __host__ __device__ __forceinline__
    vec operator+(const T & v) const {                          // Scalar arithmetic (add)
      return vec(*this) += v;
    }
    __host__ __device__ __forceinline__
    vec operator-(const T & v) const {                          // Scalar arithmetic (subtract)
      return vec(*this) -= v;
    }
    __host__ __device__ __forceinline__
    vec operator*(const T & v) const {                          // Scalar arithmetic (multiply)
      return vec(*this) *= v;
    }
    __host__ __device__ __forceinline__
    vec operator/(const T & v) const {                          // Scalar arithmetic (divide)
      return vec(*this) /= v;
    }
    __host__ __device__ __forceinline__
    vec operator>(const T & v) const {                          // Scalar arithmetic (greater than)
      return vec(*this) >= v;
    }
    __host__ __device__ __forceinline__
    vec operator<(const T & v) const {                          // Scalar arithmetic (less than)
      return vec(*this) <= v;
    }
    __host__ __device__ __forceinline__
    vec operator&(const T & v) const {                          // Scalar arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    __host__ __device__ __forceinline__
    vec operator|(const T & v) const {                          // Scalar arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    __host__ __device__ __forceinline__
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return vec(*this) += v;
    }
    __host__ __device__ __forceinline__
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return vec(*this) -= v;
    }
    __host__ __device__ __forceinline__
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return vec(*this) *= v;
    }
    __host__ __device__ __forceinline__
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return vec(*this) /= v;
    }
    __host__ __device__ __forceinline__
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return vec(*this) >= v;
    }
    __host__ __device__ __forceinline__
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return vec(*this) <= v;
    }
    __host__ __device__ __forceinline__
    vec operator&(const vec & v) const {                        // Vector arithmetic (bitwise and)
      return vec(*this) &= v;
    }
    __host__ __device__ __forceinline__
    vec operator|(const vec & v) const {                        // Vector arithmetic (bitwise or)
      return vec(*this) |= v;
    }
    __host__ __device__ __forceinline__
    vec operator-() const {                                     // Vector arithmetic (negation)
      vec temp;
      Unroll<Ops::Negate<T>,T,N>::loop(temp,data);
      return temp;
    }
    __host__ __device__ __forceinline__
    T &operator[](int i) {                                      // Indexing (lvalue)
      return data[i];
    }
    __host__ __device__ __forceinline__
    const T &operator[](int i) const {                          // Indexing (rvalue)
      return data[i];
    }
    __host__ __device__ __forceinline__
    operator       T* ()       {return data;}                   // Type-casting (lvalue)
    __host__ __device__ __forceinline__
    friend T min(const vec & v) {                               // Reduce minimum
      T temp;
      for (int i=0; i<N; i++) temp = temp < v[i] ? temp : v[i];
      return temp;
    }
    __host__ __device__ __forceinline__
    friend T max(const vec & v) {                               // Reduce maximum
      T temp;
      for (int i=0; i<N; i++) temp = temp > v[i] ? temp : v[i];
      return temp;
    }  __host__ __device__ __forceinline__
    operator const T* () const {return data;}                   // Type-casting (rvalue)
    __host__ __device__ __forceinline__
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<N; i++) s << v[i] << ' ';
      return s;
    }
    __host__ __device__ __forceinline__
    friend T sum(const vec & v) {                               // Sum vector
      return Unroll<Ops::Add<T>,T,N>::reduce(v);
    }
    __host__ __device__ __forceinline__
    friend T norm(const vec & v) {                              // L2 norm squared
      return sum(v * v);
    }
    __host__ __device__ __forceinline__
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
      return temp;
    }
    __host__ __device__ __forceinline__
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      vec temp;
      for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
      return temp;
    }
    __host__ __device__ __forceinline__
    friend T min(const vec & v) {                               // Reduce minimum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp < v[i] ? temp : v[i];
      return temp;
    }
    __host__ __device__ __forceinline__
    friend T max(const vec & v) {                               // Reduce maximum
      T temp = v[0];
      for (int i=1; i<N; i++) temp = temp > v[i] ? temp : v[i];
      return temp;
    }
    __device__ __forceinline__
    friend vec abs(const vec & v) {                             // Absolute value
      vec temp;
      Unroll<Ops::Abs<T>,T,N>::loop(temp,v);
      return temp;
    }
    __device__ __forceinline__
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
      vec temp;
      Unroll<Ops::Rsqrt<T>,T,N>::loop(temp,v);
      return temp;
    }
    __host__ __device__ __forceinline__
    friend vec sin(const vec & v) {                             // Sine function
      vec temp;
      Unroll<Ops::Sin<T>,T,N>::loop(temp,v);
      return temp;
    }
    __host__ __device__ __forceinline__
    friend vec cos(const vec & v) {                             // Cosine function
      vec temp;
      Unroll<Ops::Cos<T>,T,N>::loop(temp,v);
      return temp;
    }
    __host__ __device__ __forceinline__
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      Unroll<Ops::SinCos<T>,T,N>::loop(s,c,v);
    }
    __host__ __device__ __forceinline__
    friend vec exp(const vec & v) {                             // Exponential function
      vec temp;
      Unroll<Ops::Exp<T>,T,N>::loop(temp,v);
      return temp;
    }
    __host__ __device__ __forceinline__
    friend int wrap(vec & v, const T & w) {                     // Wrap around periodic boundary
      int iw = 0;
      for (int i=0; i<N; i++) {
	if(v[i] < -w / 2) {
	  v[i] += w;
	  iw |= 1 << i;
	}
	if(v[i] >  w / 2) {
	  v[i] -= w;
	  iw |= 1 << i;
	}
      }
      return iw;
    }
    __host__ __device__ __forceinline__
    friend void unwrap(vec & v, const T & w, const int & iw) {  // Undo wrap around periodic boundary
      for (int i=0; i<N; i++) {
	if((iw >> i) & 1) v[i] += (v[i] > 0 ? -w : w);
      }
    }
  };
#endif

#if __AVX512F__ | __MIC__
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for AVX512/MIC")
#endif
#include <immintrin.h>
  template<>
  class vec<16,float> {
  private:
    union {
      __m512 data;
      float array[16];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const float & v) {                                      // Copy constructor (scalar)
      data = _mm512_set1_ps(v);
    }
    vec(const __m512 & v) {                                     // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const float & a, const float & b, const float & c, const float & d,
	const float & e, const float & f, const float & g, const float & h,
	const float & i, const float & j, const float & k, const float & l,
	const float & m, const float & n, const float & o, const float & p) {// Copy constructor (component-wise)
      data = _mm512_setr_ps(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
    }
    ~vec(){}                                                    // Destructor
    operator __m512() const {                                   // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const float & v) {                     // Scalar assignment
      data = _mm512_set1_ps(v);
      return *this;
    }
    const vec &operator=(const __m512 & v) {                    // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm512_add_ps(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm512_sub_ps(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm512_mul_ps(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm512_div_ps(data,v);
      return *this;
    }
    const vec &operator&=(const __mmask16 & v) {                // Vector compound assignment (and)
      data = _mm512_mask_mov_ps(_mm512_setzero_ps(),v,data);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm512_add_ps(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm512_sub_ps(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm512_mul_ps(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm512_div_ps(data,v);
    }
    __mmask16 operator>(const vec & v) const {                  // Vector arithmetic (greater than)
      return _mm512_cmp_ps_mask(data,v,_MM_CMPINT_GT);
    }
    __mmask16 operator<(const vec & v) const {                  // Vector arithmetic (less than)
      return _mm512_cmp_ps_mask(data,v,_MM_CMPINT_LT);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm512_sub_ps(_mm512_setzero_ps(),data);
    }
    float &operator[](int i) {                                  // Indexing (lvalue)
      return array[i];
    }
    const float &operator[](int i) const {                      // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<16; i++) s << v[i] << ' ';
      return s;
    }
    friend float sum(const vec & v) {                           // Sum vector
      return _mm512_reduce_add_ps(v);
    }
    friend float norm(const vec & v) {                          // L2 norm squared
      __m512 temp = _mm512_mul_ps(v,v);
      return _mm512_reduce_add_ps(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm512_min_ps(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm512_max_ps(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
#ifdef __MIC__
      vec temp = _mm512_rsqrt23_ps(v);
#else
      vec temp = _mm512_rsqrt14_ps(v);
#endif
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm512_div_ps(one,_mm512_sqrt_ps(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm512_sin_ps(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm512_cos_ps(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm512_sincos_ps(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm512_exp_ps(v);
    }
  };

  template<>
  class vec<8,double> {
  private:
    union {
      __m512d data;
      double array[8];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const double & v) {                                     // Copy constructor (scalar)
      data = _mm512_set1_pd(v);
    }
    vec(const __m512d & v) {                                    // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const double & a, const double & b, const double & c, const double & d,
	const double & e, const double & f, const double & g, const double & h) {// Copy constructor (component-wise)
      data = _mm512_setr_pd(a,b,c,d,e,f,g,h);
    }
    ~vec(){}                                                    // Destructor
    operator __m512d() const {                                  // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const double & v) {                    // Scalar assignment
      data = _mm512_set1_pd(v);
      return *this;
    }
    const vec &operator=(const __m512d & v) {                   // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm512_add_pd(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm512_sub_pd(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm512_mul_pd(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm512_div_pd(data,v);
      return *this;
    }
    const vec &operator&=(const __mmask8 & v) {                 // Vector compound assignment (and)
      data = _mm512_mask_mov_pd(_mm512_setzero_pd(),v,data);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm512_add_pd(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm512_sub_pd(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm512_mul_pd(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm512_div_pd(data,v);
    }
    __mmask8 operator>(const vec & v) const {                   // Vector arithmetic (greater than)
      return _mm512_cmp_pd_mask(data,v,_MM_CMPINT_GT);
    }
    __mmask8 operator<(const vec & v) const {                   // Vector arithmetic (less than)
      return _mm512_cmp_pd_mask(data,v,_MM_CMPINT_LT);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm512_sub_pd(_mm512_setzero_pd(),data);
    }
    double &operator[](int i) {                                 // Indexing (lvalue)
      return array[i];
    }
    const double &operator[](int i) const {                     // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<8; i++) s << v[i] << ' ';
      return s;
    }
    friend double sum(const vec & v) {                          // Sum vector
      return _mm512_reduce_add_pd(v);
    }
    friend double norm(const vec & v) {                         // L2 norm squared
      __m512d temp = _mm512_mul_pd(v,v);;
      return _mm512_reduce_add_pd(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm512_min_pd(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm512_max_pd(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON
#ifdef __MIC__
      vec temp = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(v)));
#else
      vec temp = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(v)));
#endif
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm512_div_pd(one,_mm512_sqrt_pd(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm512_sin_pd(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm512_cos_pd(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm512_sincos_pd(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm512_exp_pd(v);
    }
  };
#endif

#ifdef __AVX__
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for AVX")
#endif
#include <immintrin.h>
  template<>
  class vec<8,float> {
  private:
    union {
      __m256 data;
      float array[8];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const float & v) {                                      // Copy constructor (scalar)
      data = _mm256_set1_ps(v);
    }
    vec(const __m256 & v) {                                     // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const float & a, const float & b, const float & c, const float & d,
        const float & e, const float & f, const float & g, const float & h) {// Copy constructor (component-wise)
      data = _mm256_setr_ps(a,b,c,d,e,f,g,h);
    }
    ~vec(){}                                                    // Destructor
    operator __m256() const {                                   // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const float & v) {                     // Scalar assignment
      data = _mm256_set1_ps(v);
      return *this;
    }
    const vec &operator=(const __m256 & v) {                    // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm256_add_ps(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm256_sub_ps(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm256_mul_ps(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm256_div_ps(data,v);
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = _mm256_and_ps(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm256_add_ps(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm256_sub_ps(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm256_mul_ps(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm256_div_ps(data,v);
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return _mm256_cmp_ps(data,v,_CMP_GT_OQ);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return _mm256_cmp_ps(data,v,_CMP_LT_OQ);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm256_sub_ps(_mm256_setzero_ps(),data);
    }
    float &operator[](int i) {                                  // Indexing (lvalue)
      return array[i];
    }
    const float &operator[](int i) const {                      // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<8; i++) s << v[i] << ' ';
      return s;
    }
    friend float sum(const vec & v) {                           // Sum vector
      return _mm256_reduce_add_ps(v);
    }
    friend float norm(const vec & v) {                          // L2 norm squared
      __m256 temp = _mm256_mul_ps(v,v);
      return _mm256_reduce_add_ps(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm256_min_ps(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm256_max_ps(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = _mm256_rsqrt_ps(v);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm256_div_ps(one,_mm256_sqrt_ps(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm256_sin_ps(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm256_cos_ps(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm256_sincos_ps(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm256_exp_ps(v);
    }
  };

  template<>
  class vec<4,double> {
  private:
    union {
      __m256d data;
      double array[4];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const double & v) {                                     // Copy constructor (scalar)
      data = _mm256_set1_pd(v);
    }
    vec(const __m256d & v) {                                    // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const double & a, const double & b, const double & c, const double & d) {// Copy constructor (component-wise)
      data = _mm256_setr_pd(a,b,c,d);
    }
    ~vec(){}                                                    // Destructor
    operator __m256d() const {                                  // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const double & v) {                    // Scalar assignment
      data = _mm256_set1_pd(v);
      return *this;
    }
    const vec &operator=(const __m256d & v) {                   // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm256_add_pd(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm256_sub_pd(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm256_mul_pd(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm256_div_pd(data,v);
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = _mm256_and_pd(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm256_add_pd(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm256_sub_pd(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm256_mul_pd(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm256_div_pd(data,v);
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return _mm256_cmp_pd(data,v,_CMP_GT_OQ);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return _mm256_cmp_pd(data,v,_CMP_LT_OQ);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm256_sub_pd(_mm256_setzero_pd(),data);
    }
    double &operator[](int i) {                                 // Indexing (lvalue)
      return array[i];
    }
    const double &operator[](int i) const {                     // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<4; i++) s << v[i] << ' ';
      return s;
    }
    friend double sum(const vec & v) {                          // Sum vector
      return _mm256_reduce_add_pd(v);
    }
    friend double norm(const vec & v) {                         // L2 norm squared
      __m256d temp = _mm256_mul_pd(v,v);
      return _mm256_reduce_add_pd(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm256_min_pd(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm256_max_pd(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(v)));
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm256_div_pd(one,_mm256_sqrt_pd(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm256_sin_pd(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm256_cos_pd(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm256_sincos_pd(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm256_exp_pd(v);
    }
  };
#endif

#ifdef __bgq__
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for BG/Q")
#endif
  template<>
  class vec<4,double> {
  private:
    vector4double data;
  public:
    vec(){}                                                     // Default constructor
    vec(const double & v) {                                     // Copy constructor (scalar)
      vector4double temp = {v};
      data = temp;
    }
    vec(const vector4double & v) {                              // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const double & a, const double & b, const double & c, const double & d) {// Copy constructor (component-wise)
      vector4double temp = {a,b,c,d};
      data = temp;
    }
    ~vec(){}                                                    // Destructor
    operator vector4double() const {                            // Type cast to SIMD register
      return data;
    }
    const vec &operator=(double & v) {                          // Scalar assignment
      vector4double temp = {v};
      data = temp;
      return *this;
    }
    const vec &operator=(const vector4double & v) {             // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = vec_add(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = vec_sub(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = vec_mul(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = vec_swdiv_nochk(data,v);
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = vec_and(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return vec_add(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return vec_sub(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return vec_mul(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return vec_swdiv_nochk(data,v);
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return vec_cmpgt(data,v);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return vec_cmplt(data,v);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return vec_sub((vector4double)(0),data);
    }
    double &operator[](int i) {                                 // Indexing (lvalue)
      return array[i];
    }
    const double &operator[](int i) const {                     // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<4; i++) s << vec_extract(v,i) << ' ';
      return s;
    }
    friend double sum(const vec & v) {                          // Sum vector
      double temp = 0;
      for (int i=0; i<4; i++) temp += v[i];
      return temp;
    }
    friend double norm(const vec & v) {                         // L2 norm squared
      double temp = 0;
      for (int i=0; i<4; i++) temp += v[i] * v[i];
      return temp;
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      vec temp;
      for (int i=0; i<4; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
      return temp;
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      vec temp;
      for (int i=0; i<4; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
      return temp;
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = vec_rsqrtes(v);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return one / vec_sqrt(v);
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return sind4(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return cosd4(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      sincosd4(v, s, c);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return expd4(v);
    }
  };
#endif

#ifdef __SSE__
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for SSE")
#endif
#include <pmmintrin.h>
  template<>
  class vec<4,float> {
  private:
    union {
      __m128 data;
      float array[4];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const float & v) {                                      // Copy constructor (scalar)
      data = _mm_set1_ps(v);
    }
    vec(const __m128 & v) {                                     // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const float & a, const float & b, const float & c, const float & d) {// Copy constructor (component-wise)
      data = _mm_setr_ps(a,b,c,d);
    }
    ~vec(){}                                                    // Destructor
    operator __m128() const {                                   // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const float & v) {                     // Scalar assignment
      data = _mm_set1_ps(v);
      return *this;
    }
    const vec &operator=(const __m128 & v) {                    // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm_add_ps(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm_sub_ps(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm_mul_ps(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm_div_ps(data,v);
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = _mm_and_ps(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm_add_ps(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm_sub_ps(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm_mul_ps(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm_div_ps(data,v);
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return _mm_cmpgt_ps(data,v);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return _mm_cmplt_ps(data,v);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm_sub_ps(_mm_setzero_ps(),data);
    }
    float &operator[](int i) {                                  // Indexing (lvalue)
      return array[i];
    }
    const float &operator[](int i) const {                      // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<4; i++) s << v[i] << ' ';
      return s;
    }
    friend float sum(const vec & v) {                           // Sum vector
      return _mm_reduce_add_ps(v);
    }
    friend float norm(const vec & v) {                          // L2 norm squared
      __m128 temp = _mm_mul_ps(v,v);
      return _mm_reduce_add_ps(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm_min_ps(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm_max_ps(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = _mm_rsqrt_ps(v);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm_div_ps(one,_mm_sqrt_ps(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm_sin_ps(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm_cos_ps(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm_sincos_ps(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm_exp_ps(v);
    }
  };

  template<>
  class vec<2,double> {
  private:
    union {
      __m128d data;
      double array[2];
    };
  public:
    vec(){}                                                     // Default constructor
    vec(const double & v) {                                     // Copy constructor (scalar)
      data = _mm_set1_pd(v);
    }
    vec(const __m128d & v) {                                    // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const double & a, const double & b) {                   // Copy constructor (component-wise)
      data = _mm_setr_pd(a,b);
    }
    ~vec(){}                                                    // Destructor
    operator __m128d() const {                                   // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const double & v) {                    // Scalar assignment
      data = _mm_set1_pd(v);
      return *this;
    }
    const vec &operator=(const __m128d & v) {                   // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm_add_pd(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm_sub_pd(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm_mul_pd(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm_div_pd(data,v);
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = _mm_and_pd(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm_add_pd(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm_sub_pd(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm_mul_pd(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm_div_pd(data,v);
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return _mm_cmpgt_pd(data,v);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return _mm_cmplt_pd(data,v);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm_sub_pd(_mm_setzero_pd(),data);
    }
    double &operator[](int i) {                                 // Indexing (lvalue)
      return array[i];
    }
    const double &operator[](int i) const {                     // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<2; i++) s << v[i] << ' ';
      return s;
    }
    friend double sum(const vec & v) {                          // Sum vector
      return _mm_reduce_add_pd(v);
    }
    friend double norm(const vec & v) {                         // L2 norm squared
      __m128d temp = _mm_mul_pd(v,v);
      return _mm_reduce_add_pd(temp);
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm_min_pd(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm_max_pd(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(v)));
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return _mm_div_pd(one,_mm_sqrt_pd(v));
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      return _mm_sin_pd(v);
    }
    friend vec cos(const vec & v) {                             // Cosine function
      return _mm_cos_pd(v);
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s = _mm_sincos_pd(&c.data, v);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      return _mm_exp_pd(v);
    }
  };
#endif

#if defined __sparc_v9__ & EXAFMM_USE_SIMD
#if EXAFMM_VEC_VERBOSE
#pragma message("Overloading vector operators for SPARC")
#endif
#include <emmintrin.h>

  template<>
  class vec<2,double> {
  private:
    __m128d data;
  public:
    vec(){}                                                     // Default constructor
    vec(const double & v) {                                     // Copy constructor (scalar)
      data = _mm_set_pd(v,v);
    }
    vec(const __m128d & v) {                                    // Copy constructor (SIMD register)
      data = v;
    }
    vec(const vec & v) {                                        // Copy constructor (vector)
      data = v;
    }
    vec(const double & a, const double & b) {                   // Copy constructor (component-wise)
      data = _mm_set_pd(b,a);
    }
    ~vec(){}                                                    // Destructor
    operator __m128d() const {                                  // Type cast to SIMD register
      return data;
    }
    const vec &operator=(const double & v) {                    // Scalar assignment
      data = _mm_set_pd(v,v);
      return *this;
    }
    const vec &operator=(const __m128d & v) {                   // SIMD register assignment
      data = v;
      return *this;
    }
    const vec &operator=(const vec & v) {                       // Vector assignment
      data = v;
      return *this;
    }
    const vec &operator+=(const vec & v) {                      // Vector compound assignment (add)
      data = _mm_add_pd(data,v);
      return *this;
    }
    const vec &operator-=(const vec & v) {                      // Vector compound assignment (subtract)
      data = _mm_sub_pd(data,v);
      return *this;
    }
    const vec &operator*=(const vec & v) {                      // Vector compound assignment (multiply)
      data = _mm_mul_pd(data,v);
      return *this;
    }
    const vec &operator/=(const vec & v) {                      // Vector compound assignment (divide)
      data = _mm_mul_pd(data,_fjsp_rcpa_v2r8(v));
      return *this;
    }
    const vec &operator&=(const vec & v) {                      // Vector compound assignment (bitwise and)
      data = _mm_and_pd(data,v);
      return *this;
    }
    vec operator+(const vec & v) const {                        // Vector arithmetic (add)
      return _mm_add_pd(data,v);
    }
    vec operator-(const vec & v) const {                        // Vector arithmetic (subtract)
      return _mm_sub_pd(data,v);
    }
    vec operator*(const vec & v) const {                        // Vector arithmetic (multiply)
      return _mm_mul_pd(data,v);
    }
    vec operator/(const vec & v) const {                        // Vector arithmetic (divide)
      return _mm_mul_pd(data,_fjsp_rcpa_v2r8(v));
    }
    vec operator>(const vec & v) const {                        // Vector arithmetic (greater than)
      return _mm_cmpgt_pd(data,v);
    }
    vec operator<(const vec & v) const {                        // Vector arithmetic (less than)
      return _mm_cmplt_pd(data,v);
    }
    vec operator-() const {                                     // Vector arithmetic (negation)
      return _mm_sub_pd(_mm_setzero_pd(),data);
    }
    double &operator[](int i) {                                 // Indexing (lvalue)
      return array[i];
    }
    const double &operator[](int i) const {                     // Indexing (rvalue)
      return array[i];
    }
    friend std::ostream &operator<<(std::ostream & s, const vec & v) {// Component-wise output stream
      for (int i=0; i<2; i++) s << v[i] << ' ';
      return s;
    }
    friend double sum(const vec & v) {                          // Sum vector
      return v[0] + v[1];
    }
    friend double norm(const vec & v) {                         // L2 norm squared
      return v[0] * v[0] + v[1] * v[1];
    }
    friend vec min(const vec & v, const vec & w) {              // Element-wise minimum
      return _mm_min_pd(v,w);
    }
    friend vec max(const vec & v, const vec & w) {              // Element-wise maximum
      return _mm_max_pd(v,w);
    }
    friend vec rsqrt(const vec & v) {                           // Reciprocal square root
#if EXAFMM_VEC_NEWTON                                           // Switch on Newton-Raphson correction
      vec temp = _fjsp_rsqrta_v2r8(v);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      temp *= (temp * temp * v - 3.0f) * (-0.5f);
      return temp;
#else
      vec one = 1;
      return one / _mm_sqrt_pd(v);
#endif
    }
    friend vec sin(const vec & v) {                             // Sine function
      vec temp;
      temp[0] = std::sin(v[0]);
      temp[1] = std::sin(v[1]);
      return temp;
    }
    friend vec cos(const vec & v) {                             // Cosine function
      vec temp;
      temp[0] = std::cos(v[0]);
      temp[1] = std::cos(v[1]);
      return temp;
    }
    friend void sincos(vec & s, vec & c, const vec & v) {       // Sine & cosine function
      s[0] = std::sin(v[0]);
      s[1] = std::sin(v[1]);
      c[0] = std::cos(v[0]);
      c[1] = std::cos(v[1]);
    }
    friend vec exp(const vec & v) {                             // Exponential function
      vec temp;
      temp[0] = std::exp(v[0]);
      temp[1] = std::exp(v[2]);
      return temp;
    }
  };
#endif

}
#endif
