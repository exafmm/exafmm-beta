#ifndef vec_h
#define vec_h
#define for_i for( int i=0; i!=N; ++i )
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
  friend T norm(const vec &b) {
    T c=0;
    for_i c+=b[i]*b[i];
    return c;
  }

};

#undef for_i
#endif
