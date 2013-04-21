#ifndef kahan_h
#define kahan_h
#include <iostream>
//! Operator overloading for Kahan summation
template<typename T>
struct kahan {
  T s;
  T c;
  kahan(){}                                                     // Default constructor
  kahan(const T &v) {                                           // Copy constructor (scalar)
    s = v;
    c = 0;
  }
  kahan(const kahan &v) {                                       // Copy constructor (structure)
    s = v.s;
    c = v.c;
  }
  ~kahan(){}                                                    // Destructor
  const kahan &operator=(const T v) {                           // Scalar assignment
    s = v;
    c = 0;
    return *this;
  }
  const kahan &operator+=(const T v) {                          // Scalar compound assignment (add)
    T y = v - c;
    T t = s + y;
    c = (t - s) - y;
    s = t;
    return *this;
  }
  const kahan &operator-=(const T v) {                          // Scalar compound assignment (subtract)
    T y = - v - c;
    T t = s + y;
    c = (t - s) - y;
    s = t; 
    return *this;
  }
  const kahan &operator*=(const T v) {                          // Scalar compound assignment (multiply)
    s *= v;
    c *= v;
    return *this;
  }
  const kahan &operator/=(const T v) {                          // Scalar compound assignment (divide)
    s /= v;
    c /= v;
    return *this;
  }
  const kahan &operator&=(const T v) {                          // Scalar compound assignment (bitwise and)
    s &= v;
    c &= v;
    return *this;
  }
  const kahan &operator|=(const T v) {                          // Scalar compound assignment (bitwise or)
    s |= v;
    c |= v;
    return *this;
  }
  const kahan &operator=(const kahan &v) {                      // Vector assignment
    s = v.s;
    c = v.c;
    return *this;
  }
  const kahan &operator+=(const kahan &v) {                     // Vector compound assignment (add)
    T y = v.s - c;
    T t = s + y;
    c = (t - s) - y;
    s = t;
    y = v.c - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
    return *this;
  }
  const kahan &operator-=(const kahan &v) {                     // Vector compound assignment (subtract)
    T y = - v.s - c;
    T t = s + y;
    c = (t - s) - y;
    s = t;
    y = - v.c - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
    return *this;
  }
  const kahan &operator*=(const kahan &v) {                     // Vector compound assignment (multiply)
    s *= (v.s + v.c);
    c *= (v.s + v.c);
    return *this;
  }
  const kahan &operator/=(const kahan &v) {                     // Vector compound assignment (divide)
    s /= (v.s + v.c);
    c /= (v.s + v.c);
    return *this;
  }
  const kahan &operator&=(const kahan &v) {                     // Vector compound assignment (bitwise and)
    s &= v.s;
    c &= v.c;
    return *this;
  }
  const kahan &operator|=(const kahan &v) {                     // Vector compound assignment (bitwise or)
    s |= v.s;
    c |= v.c;
    return *this;
  }
  kahan operator+(const T v) const {                            // Scalar arithmetic (add)
    kahan temp;
    T y = v - c;
    T t = s + y;
    temp.c = (t - s) - y;
    temp.s = t;
    return temp;
  }
  kahan operator-(const T v) const {                            // Scalar arithmetic (subtract)
    kahan temp;
    T y = - v - c;
    T t = s + y;
    temp.c = (t - s) - y;
    temp.s = t;
    return temp;
  }
  kahan operator*(const T v) const {                            // Scalar arithmetic (multiply)
    kahan temp; 
    temp.s = s * v;
    temp.c = c * v;
    return temp;
  }
  kahan operator/(const T v) const {                            // Scalar arithmetic (divide)
    kahan temp;
    temp.s = s / v;
    temp.c = c / v;
    return temp;
  }
  kahan operator>(const T v) const {                            // Scalar arithmetic (greater than)
    kahan temp;
    temp.s = s > v;
    temp.c = c > v;
    return temp;
  }
  kahan operator<(const T v) const {                            // Scalar arithmetic (less than)
    kahan temp;
    temp.s = s < v;
    temp.c = c < v;
    return temp;
  }
  kahan operator&(const T v) const {                            // Scalar arithmetic (bitwise and)
    kahan temp;
    temp.s = s & v;
    temp.c = c & v;
    return temp;
  }
  kahan operator|(const T v) const {                            // Scalar arithmetic (bitwise or)
    kahan temp;
    temp.s = s | v;
    temp.c = c | v;
    return temp;
  }
  kahan operator+(const kahan &v) const {                       // Vector arithmetic (add)
    kahan temp;
    temp.s = s + v.s;
    temp.c = c + v.c;
    return temp;
  }
  kahan operator-(const kahan &v) const {                       // Vector arithmetic (subtract)
    kahan temp;
    temp.s = s - v.s;
    temp.c = c - v.c;
    return temp;
  }
  kahan operator*(const kahan &v) const {                       // Vector arithmetic (multiply)
    kahan temp;
    temp.s = s * (v.s + v.c);
    temp.c = c * (v.s + v.c);
    return temp;
  }
  kahan operator/(const kahan &v) const {                       // Vector arithmetic (divide)
    kahan temp;
    temp.s = s / (v.s + v.c);
    temp.c = c / (v.s + v.c);
    return temp;
  }
  kahan operator>(const kahan &v) const {                       // Vector arithmetic (greater than)
    kahan temp;
    temp.s = s > v.s;
    temp.c = c > v.c;
    return temp;
  }
  kahan operator<(const kahan &v) const {                       // Vector arithmetic (less than)
    kahan temp;
    temp.s = s < v.s;
    temp.c = c < v.c;
    return temp;
  }
  kahan operator&(const kahan &v) const {                       // Vector arithmetic (bitwise and)
    kahan temp;
    temp.s = s & v.s;
    temp.c = c & v.c;
    return temp;
  }
  kahan operator|(const kahan &v) const {                       // Vector arithmetic (bitwise or)
    kahan temp;
    temp.s = s | v.s;
    temp.c = c | v.c;
    return temp;
  }
  kahan operator-() const {                                     // Vector arithmetic (negation)
    kahan temp;
    temp.s = -s;
    temp.c = -c;
    return temp;
  }
  operator       T ()       {return s+c;}                       // Type-casting (lvalue)
  operator const T () const {return s+c;}                       // Type-casting (rvalue)
  friend std::ostream &operator<<(std::ostream &s, const kahan &v) {// Output stream
    s << (v.s + v.c);
    return s;
  }
  friend std::istream &operator>>(std::istream &s, kahan &v) {  // Input stream
    s >> v.s;
    v.c = 0;
    return s;
  }
};
#endif
