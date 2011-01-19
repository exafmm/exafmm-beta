#ifndef vec_h
#define vec_h
#define for_i for( int i=0; i!=N; ++i )
template<int N, typename T>
class vec {
private:
  T a[N];
public:
  vec(){}                                                          // Default constructor
  vec(T const &b) {                                                // Copy constructor (scalar)
    for_i a[i] = b;
  }
  vec(vec const &b) {                                              // Copy constructor (vector)
    for_i a[i] = b[i];
  }
  ~vec(){}                                                         // Destructor
  vec const &operator=(T const b) {                                // Scalar assignment
    for_i a[i] = b;
    return *this;
  }
  vec const &operator+=(T const b) {                               // Scalar compound assignment (add)
    for_i a[i] += b;
    return *this;
  }
  vec const &operator-=(T const b) {                               // Scalar compound assignment (subtract)
    for_i a[i] -= b;
    return *this;
  }
  vec const &operator*=(T const b) {                               // Scalar compound assignment (multiply)
    for_i a[i] *= b;
    return *this;
  }
  vec const &operator/=(T const b) {                               // Scalar compound assignment (divide)
    for_i a[i] /= b;
    return *this;
  }
  vec const &operator=(vec const &b) {                             // Vector assignment
    for_i a[i] = b[i];
    return *this;
  }
  vec const &operator+=(vec const &b) {                            // Vector compound assignment (add)
    for_i a[i] += b[i];
    return *this;
  }
  vec const &operator-=(vec const &b) {                            // Vector compound assignment (subtract)
    for_i a[i] -= b[i];
    return *this;
  }
  vec const &operator*=(vec const &b) {                            // Vector compound assignment (multiply)
    for_i a[i] *= b[i];
    return *this;
  }
  vec const &operator/=(vec const &b) {                            // Vector compound assignment (divide)
    for_i a[i] /= b[i];
    return *this;
  }
  vec operator+(vec const &b) const {                              // Vector arithmetic (add)
    vec c;
    for_i c[i] = a[i] + b[i];
    return c;
  }
  vec operator-(vec const &b) const {                              // Vector arithmetic (subtract)
    vec c;
    for_i c[i] = a[i] - b[i];
    return c;
  }
  vec operator*(vec const &b) const {                              // Vector arithmetic (multiply)
    vec c;
    for_i c[i] = a[i] * b[i];
    return c;
  }
  vec operator/(vec const &b) const {                              // Vector arithmetic (divide)
    vec c;
    for_i c[i] = a[i] / b[i];
    return c;
  }
  T &operator[](int i) {                                           // Indexing (lvalue)
    return a[i];
  }
  T const &operator[](int i) const {                               // Indexing (rvalue)
    return a[i];
  }
  operator       T* ()       {return a;}                           // Type-casting (lvalue)
  operator const T* () const {return a;}                           // Type-casting (rvalue)
  friend std::ostream &operator<<(std::ostream &s, vec const &a) { // Component-wise output stream
    for_i s<<a[i]<<' ';
    return s;
  }
  friend T norm(vec const &b) {
    T c=0;
    for_i c+=b[i]*b[i];
    return c;
  }
};
#undef for_i
#endif
