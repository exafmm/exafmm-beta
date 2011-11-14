/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef vec_h
#define vec_h
#define for_i for( int i=0; i!=N; ++i )
//! Custom vector type for small vectors
template<int N, typename T>
class vec {
private:
  T a[N];
public:
  vec(){}                                                          // Default constructor
  vec(const T &b) {                                                // Copy constructor (scalar)
    for_i a[i] = b;
  }
  vec(const vec &b) {                                              // Copy constructor (vector)
    for_i a[i] = b[i];
  }
  ~vec(){}                                                         // Destructor
  const vec &operator=(const T b) {                                // Scalar assignment
    for_i a[i] = b;
    return *this;
  }
  const vec &operator+=(const T b) {                               // Scalar compound assignment (add)
    for_i a[i] += b;
    return *this;
  }
  const vec &operator-=(const T b) {                               // Scalar compound assignment (subtract)
    for_i a[i] -= b;
    return *this;
  }
  const vec &operator*=(const T b) {                               // Scalar compound assignment (multiply)
    for_i a[i] *= b;
    return *this;
  }
  const vec &operator/=(const T b) {                               // Scalar compound assignment (divide)
    for_i a[i] /= b;
    return *this;
  }
  const vec &operator=(const vec &b) {                             // Vector assignment
    for_i a[i] = b[i];
    return *this;
  }
  const vec &operator+=(const vec &b) {                            // Vector compound assignment (add)
    for_i a[i] += b[i];
    return *this;
  }
  const vec &operator-=(const vec &b) {                            // Vector compound assignment (subtract)
    for_i a[i] -= b[i];
    return *this;
  }
  const vec &operator*=(const vec &b) {                            // Vector compound assignment (multiply)
    for_i a[i] *= b[i];
    return *this;
  }
  const vec &operator/=(const vec &b) {                            // Vector compound assignment (divide)
    for_i a[i] /= b[i];
    return *this;
  }
  vec operator+(const vec &b) const {                              // Vector arithmetic (add)
    vec c;
    for_i c[i] = a[i] + b[i];
    return c;
  }
  vec operator-(const vec &b) const {                              // Vector arithmetic (subtract)
    vec c;
    for_i c[i] = a[i] - b[i];
    return c;
  }
  vec operator*(const vec &b) const {                              // Vector arithmetic (multiply)
    vec c;
    for_i c[i] = a[i] * b[i];
    return c;
  }
  vec operator/(const vec &b) const {                              // Vector arithmetic (divide)
    vec c;
    for_i c[i] = a[i] / b[i];
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
    for_i s<<a[i]<<' ';
    return s;
  }
  friend T norm(const vec &b) {                                    // L2 norm squared
    T c=0;
    for_i c+=b[i]*b[i];
    return c;
  }
};

#undef for_i
#endif
