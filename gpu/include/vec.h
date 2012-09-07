#ifndef vec_h
#define vec_h
#define for_i for( int i=0; i<N; ++i )
//! Custom vector type for small vectors
template<int N, typename T>
class vec {
private:
  T a[N];
public:
  __host__ __device__
  vec(){}                                                          // Default constructor
  __host__ __device__
  vec(const T &b) {                                                // Copy constructor (scalar)
    for_i a[i] = b;
  }
  __host__ __device__
  vec(const vec &b) {                                              // Copy constructor (vector)
    for_i a[i] = b[i];
  }
  __host__ __device__
  ~vec(){}                                                         // Destructor
  __host__ __device__
  const vec &operator=(const T b) {                                // Scalar assignment
    for_i a[i] = b;
    return *this;
  }
  __host__ __device__
  const vec &operator+=(const T b) {                               // Scalar compound assignment (add)
    for_i a[i] += b;
    return *this;
  }
  __host__ __device__
  const vec &operator-=(const T b) {                               // Scalar compound assignment (subtract)
    for_i a[i] -= b;
    return *this;
  }
  __host__ __device__
  const vec &operator*=(const T b) {                               // Scalar compound assignment (multiply)
    for_i a[i] *= b;
    return *this;
  }
  __host__ __device__
  const vec &operator/=(const T b) {                               // Scalar compound assignment (divide)
    for_i a[i] /= b;
    return *this;
  }
  __host__ __device__
  const vec &operator=(const vec &b) {                             // Vector assignment
    for_i a[i] = b[i];
    return *this;
  }
  __host__ __device__
  const vec &operator+=(const vec &b) {                            // Vector compound assignment (add)
    for_i a[i] += b[i];
    return *this;
  }
  __host__ __device__
  const vec &operator-=(const vec &b) {                            // Vector compound assignment (subtract)
    for_i a[i] -= b[i];
    return *this;
  }
  __host__ __device__
  const vec &operator*=(const vec &b) {                            // Vector compound assignment (multiply)
    for_i a[i] *= b[i];
    return *this;
  }
  __host__ __device__
  const vec &operator/=(const vec &b) {                            // Vector compound assignment (divide)
    for_i a[i] /= b[i];
    return *this;
  }
  __host__ __device__
  vec operator+(const T b) const {                                 // Scalar arithmetic (add)
    vec c;
    for_i c[i] = a[i] + b;
    return c;
  }
  __host__ __device__
  vec operator-(const T b) const {                                 // Scalar arithmetic (subtract)
    vec c;
    for_i c[i] = a[i] - b;
    return c;
  }
  __host__ __device__
  vec operator*(const T b) const {                                 // Scalar arithmetic (multiply)
    vec c;
    for_i c[i] = a[i] * b;
    return c;
  }
  __host__ __device__
  vec operator/(const T b) const {                                 // Scalar arithmetic (divide)
    vec c;
    for_i c[i] = a[i] / b;
    return c;
  }
  __host__ __device__
  vec operator+(const vec &b) const {                              // Vector arithmetic (add)
    vec c;
    for_i c[i] = a[i] + b[i];
    return c;
  }
  __host__ __device__
  vec operator-(const vec &b) const {                              // Vector arithmetic (subtract)
    vec c;
    for_i c[i] = a[i] - b[i];
    return c;
  }
  __host__ __device__
  vec operator*(const vec &b) const {                              // Vector arithmetic (multiply)
    vec c;
    for_i c[i] = a[i] * b[i];
    return c;
  }
  __host__ __device__
  vec operator/(const vec &b) const {                              // Vector arithmetic (divide)
    vec c;
    for_i c[i] = a[i] / b[i];
    return c;
  }
  __host__ __device__
  T &operator[](int i) {                                           // Indexing (lvalue)
    return a[i];
  }
  __host__ __device__
  const T &operator[](int i) const {                               // Indexing (rvalue)
    return a[i];
  }
  __host__ __device__
  operator       T* ()       {return a;}                           // Type-casting (lvalue)
  __host__ __device__
  operator const T* () const {return a;}                           // Type-casting (rvalue)
  friend std::ostream &operator<<(std::ostream &s, const vec &a) { // Component-wise output stream
    for_i s<<a[i]<<' ';
    return s;
  }
  __host__ __device__
  friend T norm(const vec &b) {                                    // L2 norm squared
    T c=0;
    for_i c+=b[i]*b[i];
    return c;
  }
};

#undef for_i
#endif
