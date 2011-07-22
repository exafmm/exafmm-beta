#ifndef vec_h
#define vec_h
#define for_i for( int i=0; i<N; i++ )
template<int N, typename T>
class vec {
private:
  T a[N];
public:
  vec(){}
  vec(T const &b) {
    for_i a[i] = b;
  }
  vec(vec const &b) {
    for_i a[i] = b[i];
  }
  ~vec(){}
  vec const &operator=(T b) {
    for_i a[i] = b;
    return *this;
  }
  vec const &operator*=(T b) {
    for_i a[i] *= b;
    return *this;
  }
  vec const &operator/=(T b) {
    for_i a[i] /= b;
    return *this;
  }
  vec const &operator=(vec const &b) {
    for_i a[i] = b[i];
    return *this;
  }
  vec const &operator+=(vec const &b) {
    for_i a[i] += b[i];
    return *this;
  }
  vec const &operator-=(vec const &b) {
    for_i a[i] -= b[i];
    return *this;
  }
  T &operator[](int i) {
    return a[i];
  }
  T const &operator[](int i) const {
    return a[i];
  }
  bool operator==(T const &b) const {
    bool c = true;
    for_i c = c && (a[i] == b);
    return c;
  }
  bool operator!=(T const &b) const {
    bool c = true;
    for_i c = c && (a[i] != b);
    return c;
  }
  vec operator+(vec const &b) const {
    vec c;
    for_i c[i] = a[i] + b[i];
    return c;
  }
  vec operator-(vec const &b) const {
    vec c;
    for_i c[i] = a[i] - b[i];
    return c;
  }
  vec operator*(vec const &b) const {
    vec c;
    for_i c[i] = a[i] * b[i];
    return c;
  }
  vec operator/(vec const &b) const {
    vec c;
    for_i c[i] = a[i] / b[i];
    return c;
  }
  operator       T* ()       {return a;}
  operator const T* () const {return a;}
  friend std::ostream &operator<<(std::ostream &s, vec const &a) {
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
