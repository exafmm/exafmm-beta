#ifndef simd_h
#define simd_h
#include <iostream>

#if __MIC__
#include <immintrin.h>
class fvec16 {
private:
  __m512 data;
public:
  fvec16() {}                                                    // Default constructor
  fvec16(const float v) {                                        // Copy constructor scalar
    data = _mm512_set1_ps(v);
  }
  fvec16(const __m512 v) {                                       // Copy constructor SIMD register
    data = v;
  }
  fvec16(const fvec16 &v) {                                      // Copy constructor vector
    data = v.data;
  }
  fvec16(const float a, const float b, const float c, const float d,
         const float e, const float f, const float g, const float h,
         const float i, const float j, const float k, const float l,
         const float m, const float n, const float o, const float p) {// Copy constructor (component-wise)
    data = _mm512_setr_ps(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p);
  }
  ~fvec16(){}                                                    // Destructor
  const fvec16 &operator=(const float v) {                       // Scalar assignment
    data = _mm512_set1_ps(v);
    return *this;
  }
  const fvec16 &operator=(const fvec16 &v) {                     // Vector assignment
    data = v.data;
    return *this;
  }
  const fvec16 &operator+=(const fvec16 &v) {                    // Vector compound assignment (add)
    data = _mm512_add_ps(data,v.data);
    return *this;
  }
  const fvec16 &operator-=(const fvec16 &v) {                    // Vector compound assignment (subtract)
    data = _mm512_sub_ps(data,v.data);
    return *this;
  }
  const fvec16 &operator*=(const fvec16 &v) {                    // Vector compound assignment (multiply)
    data = _mm512_mul_ps(data,v.data);
    return *this;
  }
  const fvec16 &operator/=(const fvec16 &v) {                    // Vector compound assignment (divide)
    data = _mm512_div_ps(data,v.data);
    return *this;
  }
  const fvec16 &operator&=(const fvec16 &v) {                    // Vector compound assignment (bitwise and)
    data = _mm512_and_ps(data,v.data);
    return *this;
  }
  const fvec16 &operator|=(const fvec16 &v) {                    // Vector compound assignment (bitwise or)
    data = _mm512_or_ps(data,v.data);
    return *this;
  }
  fvec16 operator+(const fvec16 &v) const {                      // Vector arithmetic (add)
    return fvec16(_mm512_add_ps(data,v.data));
  }
  fvec16 operator-(const fvec16 &v) const {                      // Vector arithmetic (subtract)
    return fvec16(_mm512_sub_ps(data,v.data));
  }
  fvec16 operator*(const fvec16 &v) const {                      // Vector arithmetic (multiply)
    return fvec16(_mm512_mul_ps(data,v.data));
  }
  fvec16 operator/(const fvec16 &v) const {                      // Vector arithmetic (divide)
    return fvec16(_mm512_div_ps(data,v.data));
  }
  fvec16 operator>(const fvec16 &v) const {                      // Vector arithmetic (greater than)
    return fvec16(_mm512_cmp_ps(data,v.data,_CMP_GT_OQ));
  }
  fvec16 operator<(const fvec16 &v) const {                      // Vector arithmetic (less than)
    return fvec16(_mm512_cmp_ps(data,v.data,_CMP_LT_OQ));
  }
  fvec16 operator&(const fvec16 &v) const {                      // Vector arithmetic (bitwise and)
    return fvec16(_mm512_and_ps(data,v.data));
  }
  fvec16 operator|(const fvec16 &v) const {                      // Vector arithmetic (bitwise or)
    return fvec16(_mm512_or_ps(data,v.data));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const fvec16 &v) {// Component-wise output stream
    for (int i=0; i<8; i++) s << v[i] << ' ';
    return s;
  }
  friend fvec16 rsqrt(const fvec16 &v) {                         // reciprocal square root
    return fvec16(_mm512_rsqrt23_ps(v.data));
  }
};
#endif

#if __AVX__
#include <immintrin.h>
class fvec8 {
private:
  __m256 data;
public:
  fvec8() {}                                                     // Default constructor
  fvec8(const float v) {                                         // Copy constructor scalar
    data = _mm256_set1_ps(v);
  }
  fvec8(const __m256 v) {                                        // Copy constructor SIMD register
    data = v;
  }
  fvec8(const fvec8 &v) {                                        // Copy constructor vector
    data = v.data;
  }
  fvec8(const float a, const float b, const float c, const float d,
        const float e, const float f, const float g, const float h) {// Copy constructor (component-wise)
    data = _mm256_setr_ps(a,b,c,d,e,f,g,h);
  }
  ~fvec8(){}                                                     // Destructor
  const fvec8 &operator=(const float v) {                        // Scalar assignment
    data = _mm256_set1_ps(v);
    return *this;
  }
  const fvec8 &operator=(const fvec8 &v) {                       // Vector assignment
    data = v.data;
    return *this;
  }
  const fvec8 &operator+=(const fvec8 &v) {                      // Vector compound assignment (add)
    data = _mm256_add_ps(data,v.data);
    return *this;
  }
  const fvec8 &operator-=(const fvec8 &v) {                      // Vector compound assignment (subtract)
    data = _mm256_sub_ps(data,v.data);
    return *this;
  }
  const fvec8 &operator*=(const fvec8 &v) {                      // Vector compound assignment (multiply)
    data = _mm256_mul_ps(data,v.data);
    return *this;
  }
  const fvec8 &operator/=(const fvec8 &v) {                      // Vector compound assignment (divide)
    data = _mm256_div_ps(data,v.data);
    return *this;
  }
  const fvec8 &operator&=(const fvec8 &v) {                      // Vector compound assignment (bitwise and)
    data = _mm256_and_ps(data,v.data);
    return *this;
  }
  const fvec8 &operator|=(const fvec8 &v) {                      // Vector compound assignment (bitwise or)
    data = _mm256_or_ps(data,v.data);
    return *this;
  }
  fvec8 operator+(const fvec8 &v) const {                        // Vector arithmetic (add)
    return fvec8(_mm256_add_ps(data,v.data));
  }
  fvec8 operator-(const fvec8 &v) const {                        // Vector arithmetic (subtract)
    return fvec8(_mm256_sub_ps(data,v.data));
  }
  fvec8 operator*(const fvec8 &v) const {                        // Vector arithmetic (multiply)
    return fvec8(_mm256_mul_ps(data,v.data));
  }
  fvec8 operator/(const fvec8 &v) const {                        // Vector arithmetic (divide)
    return fvec8(_mm256_div_ps(data,v.data));
  }
  fvec8 operator>(const fvec8 &v) const {                        // Vector arithmetic (greater than)
    return fvec8(_mm256_cmp_ps(data,v.data,_CMP_GT_OQ));
  }
  fvec8 operator<(const fvec8 &v) const {                        // Vector arithmetic (less than)
    return fvec8(_mm256_cmp_ps(data,v.data,_CMP_LT_OQ));
  }
  fvec8 operator&(const fvec8 &v) const {                        // Vector arithmetic (bitwise and)
    return fvec8(_mm256_and_ps(data,v.data));
  }
  fvec8 operator|(const fvec8 &v) const {                        // Vector arithmetic (bitwise or)
    return fvec8(_mm256_or_ps(data,v.data));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const fvec8 &v) {// Component-wise output stream
    for (int i=0; i<8; i++) s << v[i] << ' ';
    return s;
  }
  friend fvec8 rsqrt(const fvec8 &v) {                           // reciprocal square root
    return fvec8(_mm256_rsqrt_ps(v.data));
  }
};
#endif

#if __SSE__
#include <xmmintrin.h>
class fvec4 {
private:
  __m128 data;
public:
  fvec4() {}                                                     // Default constructor
  fvec4(const float v) {                                         // Copy constructor scalar
    data = _mm_set1_ps(v);
  }
  fvec4(const __m128 v) {                                        // Copy constructor SIMD register
    data = v;
  }
  fvec4(const fvec4 &v) {                                        // Copy constructor vector
    data = v.data;
  }
  fvec4(const float a, const float b, const float c, const float d) {// Copy constructor (component-wise)
    data = _mm_setr_ps(a,b,c,d);
  }
  ~fvec4(){}                                                     // Destructor
  const fvec4 &operator=(const float v) {                        // Scalar assignment
    data = _mm_set1_ps(v);
    return *this;
  }
  const fvec4 &operator=(const fvec4 &v) {                       // Vector assignment
    data = v.data;
    return *this;
  }
  const fvec4 &operator+=(const fvec4 &v) {                      // Vector compound assignment (add)
    data = _mm_add_ps(data,v.data);
    return *this;
  }
  const fvec4 &operator-=(const fvec4 &v) {                      // Vector compound assignment (subtract)
    data = _mm_sub_ps(data,v.data);
    return *this;
  }
  const fvec4 &operator*=(const fvec4 &v) {                      // Vector compound assignment (multiply)
    data = _mm_mul_ps(data,v.data);
    return *this;
  }
  const fvec4 &operator/=(const fvec4 &v) {                      // Vector compound assignment (divide)
    data = _mm_div_ps(data,v.data);
    return *this;
  }
  const fvec4 &operator&=(const fvec4 &v) {                      // Vector compound assignment (bitwise and)
    data = _mm_add_ps(data,v.data);
    return *this;
  }
  const fvec4 &operator|=(const fvec4 &v) {                      // Vector compound assignment (bitwise or)
    data = _mm_or_ps(data,v.data);
    return *this;
  }
  fvec4 operator+(const fvec4 &v) const {                        // Vector arithmetic (add)
    return fvec4(_mm_add_ps(data,v.data));
  }
  fvec4 operator-(const fvec4 &v) const {                        // Vector arithmetic (subtract)
    return fvec4(_mm_sub_ps(data,v.data));
  }
  fvec4 operator*(const fvec4 &v) const {                        // Vector arithmetic (multiply)
    return fvec4(_mm_mul_ps(data,v.data));
  }
  fvec4 operator/(const fvec4 &v) const {                        // Vector arithmetic (divide)
    return fvec4(_mm_div_ps(data,v.data));
  }
  fvec4 operator>(const fvec4 &v) const {                        // Vector arithmetic (greater than)
    return fvec4(_mm_cmpgt_ps(data,v.data));
  }
  fvec4 operator<(const fvec4 &v) const {                        // Vector arithmetic (less than)
    return fvec4(_mm_cmplt_ps(data,v.data));
  }
  fvec4 operator&(const fvec4 &v) const {                        // Vector arithmetic (bitwise and)
    return fvec4(_mm_and_ps(data,v.data));
  }
  fvec4 operator|(const fvec4 &v) const {                        // Vector arithmetic (bitwise or)
    return fvec4(_mm_or_ps(data,v.data));
  }
  float &operator[](int i) {                                     // Indexing (lvalue)
    return ((float*)&data)[i];
  }
  const float &operator[](int i) const {                         // Indexing (rvalue)
    return ((float*)&data)[i];
  }
  friend std::ostream &operator<<(std::ostream &s, const fvec4 &v) {// Component-wise output stream
    for (int i=0; i<4; i++) s << v[i] << ' ';
    return s;
  }
  friend fvec4 rsqrt(const fvec4 &v) {                           // reciprocal square root
    return fvec4(_mm_rsqrt_ps(v.data));
  }
};
#endif

#endif
