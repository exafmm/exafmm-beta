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
  const kahan &operator>=(const T v) {                          // Scalar compound assignment (greater than)
    s >= v;
    c >= v;
    return *this;
  }
  const kahan &operator<=(const T v) {                          // Scalar compound assignment (less than)
    s <= v;
    c <= v;
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
  const kahan &operator>=(const kahan &v) {                     // Vector compound assignment (greater than)
    s >= v.s;
    c >= v.c;
    return *this;
  }
  const kahan &operator<=(const kahan &v) {                     // Vector compound assignment (less than)
    s <= (v.s + v.c);
    c <= (v.s + v.c);
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
    return kahan(*this) += v;
  }
  kahan operator-(const T v) const {                            // Scalar arithmetic (subtract)
    return kahan(*this) -= v;
  }
  kahan operator*(const T v) const {                            // Scalar arithmetic (multiply)
    return kahan(*this) *= v;
  }
  kahan operator/(const T v) const {                            // Scalar arithmetic (divide)
    return kahan(*this) /= v;
  }
  kahan operator>(const T v) const {                            // Scalar arithmetic (greater than)
    return kahan(*this) >= v;
  }
  kahan operator<(const T v) const {                            // Scalar arithmetic (less than)
    return kahan(*this) <= v;
  }
  kahan operator&(const T v) const {                            // Scalar arithmetic (bitwise and)
    return kahan(*this) &= v;
  }
  kahan operator|(const T v) const {                            // Scalar arithmetic (bitwise or)
    return kahan(*this) |= v;
  }
  kahan operator+(const kahan &v) const {                       // Vector arithmetic (add)
    return kahan(*this) += v;
  }
  kahan operator-(const kahan &v) const {                       // Vector arithmetic (subtract)
    return kahan(*this) -= v;
  }
  kahan operator*(const kahan &v) const {                       // Vector arithmetic (multiply)
    return kahan(*this) *= v;
  }
  kahan operator/(const kahan &v) const {                       // Vector arithmetic (divide)
    return kahan(*this) /= v;
  }
  kahan operator>(const kahan &v) const {                       // Vector arithmetic (greater than)
    return kahan(*this) >= v;
  }
  kahan operator<(const kahan &v) const {                       // Vector arithmetic (less than)
    return kahan(*this) <= v;
  }
  kahan operator&(const kahan &v) const {                       // Vector arithmetic (bitwise and)
    return kahan(*this) &= v;
  }
  kahan operator|(const kahan &v) const {                       // Vector arithmetic (bitwise or)
    return kahan(*this) |= v;
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
