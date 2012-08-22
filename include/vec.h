#ifndef vec_h
#define vec_h
#include <ostream>
//! Custom vector type for small vectors
template<int N, typename T>
class vec {
private:
  T a[N];
public:
  vec(){}                                                          // Default constructor
  vec(const T &b) {                                                // Copy constructor (scalar)
    for (int i=0; i<N; i++) a[i] = b;
  }
  vec(const vec &b) {                                              // Copy constructor (vector)
    for (int i=0; i<N; i++) a[i] = b[i];
  }
  ~vec(){}                                                         // Destructor
  const vec &operator=(const T b) {                                // Scalar assignment
    for (int i=0; i<N; i++) a[i] = b;
    return *this;
  }
  const vec &operator+=(const T b) {                               // Scalar compound assignment (add)
    for (int i=0; i<N; i++) a[i] += b;
    return *this;
  }
  const vec &operator-=(const T b) {                               // Scalar compound assignment (subtract)
    for (int i=0; i<N; i++) a[i] -= b;
    return *this;
  }
  const vec &operator*=(const T b) {                               // Scalar compound assignment (multiply)
    for (int i=0; i<N; i++) a[i] *= b;
    return *this;
  }
  const vec &operator/=(const T b) {                               // Scalar compound assignment (divide)
    for (int i=0; i<N; i++) a[i] /= b;
    return *this;
  }
  const vec &operator=(const vec &b) {                             // Vector assignment
    for (int i=0; i<N; i++) a[i] = b[i];
    return *this;
  }
  const vec &operator+=(const vec &b) {                            // Vector compound assignment (add)
    for (int i=0; i<N; i++) a[i] += b[i];
    return *this;
  }
  const vec &operator-=(const vec &b) {                            // Vector compound assignment (subtract)
    for (int i=0; i<N; i++) a[i] -= b[i];
    return *this;
  }
  const vec &operator*=(const vec &b) {                            // Vector compound assignment (multiply)
    for (int i=0; i<N; i++) a[i] *= b[i];
    return *this;
  }
  const vec &operator/=(const vec &b) {                            // Vector compound assignment (divide)
    for (int i=0; i<N; i++) a[i] /= b[i];
    return *this;
  }
  vec operator+(const T b) const {                                 // Scalar arithmetic (add)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] + b;
    return c;
  }
  vec operator-(const T b) const {                                 // Scalar arithmetic (subtract)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] - b;
    return c;
  }
  vec operator*(const T b) const {                                 // Scalar arithmetic (multiply)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] * b;
    return c;
  }
  vec operator/(const T b) const {                                 // Scalar arithmetic (divide)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] / b;
    return c;
  }
  vec operator+(const vec &b) const {                              // Vector arithmetic (add)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] + b[i];
    return c;
  }
  vec operator-(const vec &b) const {                              // Vector arithmetic (subtract)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] - b[i];
    return c;
  }
  vec operator*(const vec &b) const {                              // Vector arithmetic (multiply)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] * b[i];
    return c;
  }
  vec operator/(const vec &b) const {                              // Vector arithmetic (divide)
    vec c;
    for (int i=0; i<N; i++) c[i] = a[i] / b[i];
    return c;
  }
  T &operator[](int i) {                                           // Indexing (lvalue)
    return a[i];
  }
  const T &operator[](int i) const {                               // Indexing (rvalue)
    return a[i];
  }
  operator       T* ()       {return a;}                           // Type-casting (lvalue)
  operator const T* () const {return a;}                           // Type-casting (rvalue)
  friend std::ostream &operator<<(std::ostream &s, const vec &a) { // Component-wise output stream
    for (int i=0; i<N; i++) s << a[i] << ' ';
    return s;
  }
  friend T norm(const vec &b) {                                    // L2 norm squared
    T c=0;
    for (int i=0; i<N; i++) c += b[i] * b[i];
    return c;
  }
  friend vec min(const vec &b, const vec &c) {                     // Element-wise minimum
    vec d;
    for (int i=0; i<N; i++) d[i] = b[i] < c[i] ? b[i] : c[i];
    return d;
  }
  friend vec max(const vec &b, const vec &c) {                     // Element-wise maximum
    vec d;
    for (int i=0; i<N; i++) d[i] = b[i] > c[i] ? b[i] : c[i];
    return d;
  }
};

#endif
