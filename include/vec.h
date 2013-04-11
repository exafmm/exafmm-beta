#ifndef vec_h
#define vec_h
#include <ostream>
//! Custom vector type for small vectors with template specialization for MIC, AVX, SSE intrinsics
template<int N, typename T>
class vec {
private:
  T data[N];
public:
  vec(){}                                                        // Default constructor
  vec(const T &v) {                                              // Copy constructor (scalar)
    for (int i=0; i<N; i++) data[i] = v;
  }
  vec(const vec &v) {                                            // Copy constructor (vector)
    for (int i=0; i<N; i++) data[i] = v[i];
  }
  ~vec(){}                                                       // Destructor
  const vec &operator=(const T v) {                              // Scalar assignment
    for (int i=0; i<N; i++) data[i] = v;
    return *this;
  }
  const vec &operator+=(const T v) {                             // Scalar compound assignment (add)
    for (int i=0; i<N; i++) data[i] += v;
    return *this;
  }
  const vec &operator-=(const T v) {                             // Scalar compound assignment (subtract)
    for (int i=0; i<N; i++) data[i] -= v;
    return *this;
  }
  const vec &operator*=(const T v) {                             // Scalar compound assignment (multiply)
    for (int i=0; i<N; i++) data[i] *= v;
    return *this;
  }
  const vec &operator/=(const T v) {                             // Scalar compound assignment (divide)
    for (int i=0; i<N; i++) data[i] /= v;
    return *this;
  }
  const vec &operator&=(const T v) {                             // Scalar compound assignment (bitwise and)
    for (int i=0; i<N; i++) data[i] &= v;
    return *this;
  }
  const vec &operator|=(const T v) {                             // Scalar compound assignment (bitwise or)
    for (int i=0; i<N; i++) data[i] |= v;
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    for (int i=0; i<N; i++) data[i] = v[i];
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    for (int i=0; i<N; i++) data[i] += v[i];
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    for (int i=0; i<N; i++) data[i] -= v[i];
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    for (int i=0; i<N; i++) data[i] *= v[i];
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    for (int i=0; i<N; i++) data[i] /= v[i];
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    for (int i=0; i<N; i++) data[i] &= v[i];
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    for (int i=0; i<N; i++) data[i] |= v[i];
    return *this;
  }
  vec operator+(const T v) const {                               // Scalar arithmetic (add)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] + v;
    return temp;
  }
  vec operator-(const T v) const {                               // Scalar arithmetic (subtract)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] - v;
    return temp;
  }
  vec operator*(const T v) const {                               // Scalar arithmetic (multiply)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] * v;
    return temp;
  }
  vec operator/(const T v) const {                               // Scalar arithmetic (divide)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] / v;
    return temp;
  }
  vec operator>(const T v) const {                               // Scalar arithmetic (greater than)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] > v;
    return temp;
  }
  vec operator<(const T v) const {                               // Scalar arithmetic (less than)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] < v;
    return temp;
  }
  vec operator&(const T v) const {                               // Scalar arithmetic (bitwise and)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] & v;
    return temp;
  }
  vec operator|(const T v) const {                               // Scalar arithmetic (bitwise or)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] | v;
    return temp;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] + v[i];
    return temp;
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] - v[i];
    return temp;
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] * v[i];
    return temp;
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] / v[i];
    return temp;
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] > v[i];
    return temp;
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] < v[i];
    return temp;
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] & v[i];
    return temp;
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = data[i] | v[i];
    return temp;
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    vec temp;
    for (int i=0; i<N; i++) temp[i] = -data[i];
    return temp;
  }
  T &operator[](int i) {                                         // Indexing (lvalue)
    return data[i];
  }
  const T &operator[](int i) const {                             // Indexing (rvalue)
    return data[i];
  }
  operator       T* ()       {return data;}                      // Type-casting (lvalue)
  operator const T* () const {return data;}                      // Type-casting (rvalue)
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<N; i++) s << v[i] << ' ';
    return s;
  }
  friend T sum(const vec &v) {                                   // Sum vector
    T temp = 0;
    for (int i=0; i<N; i++) temp += v[i];
    return temp;
  }
  friend T norm(const vec &v) {                                  // L2 norm squared
    T temp = 0;
    for (int i=0; i<N; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    vec temp;
    for (int i=0; i<N; i++) temp[i] = v[i] < w[i] ? v[i] : w[i];
    return temp;
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    vec temp;
    for (int i=0; i<N; i++) temp[i] = v[i] > w[i] ? v[i] : w[i];
    return temp;
  }
};

#if __MIC__
#include <immintrin.h>
template<>
class vec<16,float> {
private:
  __m512 data;
public:
  vec(){}                                                        // Default constructor
  vec(const float v) {                                           // Copy constructor scalar
    data = _mm512_set1_ps(v);
  }
  vec(const __m512 v) {                                          // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  ~vec(){}                                                       // Destructor
  vec(const float a, const float b, const float c, const float d,
      const float e, const float f, const float g, const float h,
      const float i, const float j, const float k, const float l,
      const float m, const float n, const float o, const float p) {// Copy constructor (component-wise)
    data = _mm512_setr_ps(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
  }
  ~vec(){}                                                       // Destructor
  const vec &operator=(const float v) {                          // Scalar assignment
    data = _mm512_set1_ps(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm512_add_ps(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm512_sub_ps(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm512_mul_ps(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm512_div_ps(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm512_and_ps(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm512_or_ps(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm512_add_ps(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm512_sub_ps(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm512_mul_ps(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm512_div_ps(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm512_cmp_ps(data,v.data,_CMP_GT_OQ));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm512_cmp_ps(data,v.data,_CMP_LT_OQ));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm512_and_ps(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm512_or_ps(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm512_xor_ps(data,_mm512_castsi512_ps(_mm512_set1_epi32(0x80000000)));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<16; i++) s << v[i] << ' ';
    return s;
  }
  friend float sum(const vec &v) {                               // Sum vector
    float temp = 0;
    for (int i=0; i<16; i++) temp += v[i];
    return temp;
  }
  friend float norm(const vec &v) {                              // L2 norm squared
    float temp = 0;
    for (int i=0; i<16; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm512_min_ps(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm512_max_ps(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    return vec(_mm512_rsqrt23_ps(v.data));
  }
};

template<>
class vec<8,double> {
private:
  __m512d data;
public:
  vec(){}                                                        // Default constructor
  vec(const double v) {                                          // Copy constructor scalar
    data = _mm512_set1_pd(v);
  }
  vec(const __m512d v) {                                         // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  ~vec(){}                                                       // Destructor
  vec(const double a, const double b, const double c, const double d,
      const double e, const double f, const double g, const double h) {// Copy constructor (component-wise)
    data = _mm512_setr_pd(a,b,c,d,e,f,g,h);
  }
  ~vec(){}                                                       // Destructor
  const vec &operator=(const double v) {                         // Scalar assignment
    data = _mm512_set1_pd(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm512_add_pd(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm512_sub_pd(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm512_mul_pd(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm512_div_pd(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm512_and_pd(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm512_or_pd(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm512_add_pd(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm512_sub_pd(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm512_mul_pd(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm512_div_pd(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm512_cmp_pd(data,v.data,_CMP_GT_OQ));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm512_cmp_pd(data,v.data,_CMP_LT_OQ));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm512_and_pd(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm512_or_pd(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm512_xor_pd(data,_mm512_castsi512_pd(_mm512_set1_epi64(0x8000000000000000)));
  }
  double &operator[](int i) {                                    // Indexing (lvalue)
    return ((double*)&data)[i];
  }
  const double &operator[](int i) const {                        // Indexing (rvalue)
    return ((double*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<8; i++) s << v[i] << ' ';
    return s;
  }
  friend double sum(const vec &v) {                              // Sum vector
    double temp = 0;
    for (int i=0; i<8; i++) temp += v[i];
    return temp;
  }
  friend double norm(const vec &v) {                             // L2 norm squared
    double temp = 0;
    for (int i=0; i<8; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm512_min_pd(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm512_max_pd(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    vec<16,float> in(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],0,0,0,0,0,0,0,0);
    vec<16,float> temp = rsqrt(in);
    temp *= (temp * temp * in - 3.0f) * (-0.5f);
    vec<8,double> out(temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],temp[7]);
    return out;
  }
};
#endif

#if __AVX__
#include <immintrin.h>
template<>
class vec<8,float> {
private:
  __m256 data;
public:
  vec(){}                                                        // Default constructor
  vec(const float v) {                                           // Copy constructor scalar
    data = _mm256_set1_ps(v);
  }
  vec(const __m256 v) {                                          // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  ~vec(){}                                                       // Destructor
  vec(const float a, const float b, const float c, const float d,
      const float e, const float f, const float g, const float h) {// Copy constructor (component-wise)
    data = _mm256_setr_ps(a,b,c,d,e,f,g,h);
  }
  const vec &operator=(const float v) {                          // Scalar assignment
    data = _mm256_set1_ps(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm256_add_ps(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm256_sub_ps(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm256_mul_ps(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm256_div_ps(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm256_and_ps(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm256_or_ps(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm256_add_ps(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm256_sub_ps(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm256_mul_ps(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm256_div_ps(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm256_cmp_ps(data,v.data,_CMP_GT_OQ));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm256_cmp_ps(data,v.data,_CMP_LT_OQ));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm256_and_ps(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm256_or_ps(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm256_xor_ps(data,_mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<8; i++) s << v[i] << ' ';
    return s;
  }
  friend float sum(const vec &v) {                               // Sum vector
    float temp = 0;
    for (int i=0; i<8; i++) temp += v[i];
    return temp;
  }
  friend float norm(const vec &v) {                              // L2 norm squared
    float temp = 0;
    for (int i=0; i<8; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm256_min_ps(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm256_max_ps(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    return vec(_mm256_rsqrt_ps(v.data));
  }
};

template<>
class vec<4,double> {
private:
  __m256d data;
public:
  vec(){}                                                        // Default constructor
  vec(const double v) {                                          // Copy constructor scalar
    data = _mm256_set1_pd(v);
  }
  vec(const __m256d v) {                                         // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  ~vec(){}                                                       // Destructor
  vec(const double a, const double b, const double c, const double d) {// Copy constructor (component-wise)
    data = _mm256_setr_pd(a,b,c,d);
  }
  const vec &operator=(const double v) {                         // Scalar assignment
    data = _mm256_set1_pd(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm256_add_pd(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm256_sub_pd(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm256_mul_pd(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm256_div_pd(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm256_and_pd(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm256_or_pd(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm256_add_pd(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm256_sub_pd(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm256_mul_pd(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm256_div_pd(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm256_cmp_pd(data,v.data,_CMP_GT_OQ));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm256_cmp_pd(data,v.data,_CMP_LT_OQ));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm256_and_pd(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm256_or_pd(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm256_xor_pd(data,_mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000)));
  }
  double &operator[](int i) {                                    // Indexing (lvalue)
    return ((double*)&data)[i];
  }
  const double &operator[](int i) const {                        // Indexing (rvalue)
    return ((double*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<4; i++) s << v[i] << ' ';
    return s;
  }
  friend double sum(const vec &v) {                              // Sum vector
    double temp = 0;
    for (int i=0; i<4; i++) temp += v[i];
    return temp;
  }
  friend double norm(const vec &v) {                             // L2 norm squared
    double temp = 0;
    for (int i=0; i<4; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm256_min_pd(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm256_max_pd(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    //vec<8,float> in(v[0],v[1],v[2],v[3],0,0,0,0);
    //vec<8,float> temp = rsqrt(in);
    //temp *= (temp * temp * in - 3.0f) * (-0.5f);
    //vec<4,double> out(temp[0],temp[1],temp[2],temp[3]);
    //return out;
    vec one = 1;
    return vec(_mm256_div_pd(one.data,_mm256_sqrt_pd(v.data)));
  }
};
#endif

#if __SSE__
#include <xmmintrin.h>
template<>
class vec<4,float> {
private:
  __m128 data;
public:
  vec(){}                                                        // Default constructor
  vec(const float v) {                                           // Copy constructor scalar
    data = _mm_set1_ps(v);
  }
  vec(const __m128 v) {                                          // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  vec(const float a, const float b, const float c, const float d) {// Copy constructor (component-wise)
    data = _mm_setr_ps(a,b,c,d);
  }
  ~vec(){}                                                       // Destructor
  const vec &operator=(const float v) {                          // Scalar assignment
    data = _mm_set1_ps(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm_add_ps(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm_sub_ps(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm_mul_ps(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm_div_ps(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm_and_ps(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm_or_ps(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm_add_ps(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm_sub_ps(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm_mul_ps(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm_div_ps(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm_cmpgt_ps(data,v.data));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm_cmplt_ps(data,v.data));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm_and_ps(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm_or_ps(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm_xor_ps(data,_mm_castsi128_ps(_mm_set1_epi32(0x80000000)));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<4; i++) s << v[i] << ' ';
    return s;
  }
  friend float sum(const vec &v) {                               // Sum vector
    float temp = 0;
    for (int i=0; i<4; i++) temp += v[i];
    return temp;
  }
  friend float norm(const vec &v) {                              // L2 norm squared
    float temp = 0;
    for (int i=0; i<4; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm_min_ps(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm_max_ps(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    return vec(_mm_rsqrt_ps(v.data));
  }
};

template<>
class vec<2,double> {
private:
  __m128d data;
public:
  vec(){}                                                        // Default constructor
  vec(const double v) {                                          // Copy constructor scalar
    data = _mm_set1_pd(v);
  }
  vec(const __m128d v) {                                         // Copy constructor SIMD register
    data = v;
  }
  vec(const vec &v) {                                            // Copy constructor vector
    data = v.data;
  }
  vec(const double a, const double b) {                          // Copy constructor (component-wise)
    data = _mm_setr_pd(a,b);
  }
  ~vec(){}                                                       // Destructor
  const vec &operator=(const double v) {                         // Scalar assignment
    data = _mm_set1_pd(v);
    return *this;
  }
  const vec &operator=(const vec &v) {                           // Vector assignment
    data = v.data;
    return *this;
  }
  const vec &operator+=(const vec &v) {                          // Vector compound assignment (add)
    data = _mm_add_pd(data,v.data);
    return *this;
  }
  const vec &operator-=(const vec &v) {                          // Vector compound assignment (subtract)
    data = _mm_sub_pd(data,v.data);
    return *this;
  }
  const vec &operator*=(const vec &v) {                          // Vector compound assignment (multiply)
    data = _mm_mul_pd(data,v.data);
    return *this;
  }
  const vec &operator/=(const vec &v) {                          // Vector compound assignment (divide)
    data = _mm_div_pd(data,v.data);
    return *this;
  }
  const vec &operator&=(const vec &v) {                          // Vector compound assignment (bitwise and)
    data = _mm_and_pd(data,v.data);
    return *this;
  }
  const vec &operator|=(const vec &v) {                          // Vector compound assignment (bitwise or)
    data = _mm_or_pd(data,v.data);
    return *this;
  }
  vec operator+(const vec &v) const {                            // Vector arithmetic (add)
    return vec(_mm_add_pd(data,v.data));
  }
  vec operator-(const vec &v) const {                            // Vector arithmetic (subtract)
    return vec(_mm_sub_pd(data,v.data));
  }
  vec operator*(const vec &v) const {                            // Vector arithmetic (multiply)
    return vec(_mm_mul_pd(data,v.data));
  }
  vec operator/(const vec &v) const {                            // Vector arithmetic (divide)
    return vec(_mm_div_pd(data,v.data));
  }
  vec operator>(const vec &v) const {                            // Vector arithmetic (greater than)
    return vec(_mm_cmpgt_pd(data,v.data));
  }
  vec operator<(const vec &v) const {                            // Vector arithmetic (less than)
    return vec(_mm_cmplt_pd(data,v.data));
  }
  vec operator&(const vec &v) const {                            // Vector arithmetic (bitwise and)
    return vec(_mm_and_pd(data,v.data));
  }
  vec operator|(const vec &v) const {                            // Vector arithmetic (bitwise or)
    return vec(_mm_or_pd(data,v.data));
  }
  vec operator-() const {                                        // Vector arithmetic (negation)
    return _mm_xor_pd(data,_mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000)));
  }
  double &operator[](int i) {                                    // Indexing (lvalue)
    return ((double*)&data)[i];
  }
  const double &operator[](int i) const {                        // Indexing (rvalue)
    return ((double*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const vec &v) {// Component-wise output stream
    for (int i=0; i<2; i++) s << v[i] << ' ';
    return s;
  }
  friend double sum(const vec &v) {                              // Sum vector
    double temp = 0;
    for (int i=0; i<2; i++) temp += v[i];
    return temp;
  }
  friend double norm(const vec &v) {                             // L2 norm squared
    double temp = 0;
    for (int i=0; i<2; i++) temp += v[i] * v[i];
    return temp;
  }
  friend vec min(const vec &v, const vec &w) {                   // Element-wise minimum
    return vec(_mm_min_pd(v.data,w.data));
  }
  friend vec max(const vec &v, const vec &w) {                   // Element-wise maximum
    return vec(_mm_max_pd(v.data,w.data));
  }
  friend vec rsqrt(const vec &v) {                               // Reciprocal square root
    vec<4,float> in(v[0],v[1],0,0);
    vec<4,float> temp = rsqrt(in);
    temp *= (temp * temp * in - 3.0f) * (-0.5f);
    vec<2,double> out(temp[0],temp[1]);
    return out;
  }
};
#endif

#endif
