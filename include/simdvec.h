#ifndef simdvec_h
#define simdvec_h

namespace exafmm {
  template<typename T, int D, int N>
  struct SIMD {
    static inline T setBody(B_iter, int) {
      T v;
      return v;
    }
    static inline T setIndex(int) {
      T v;
      return v;
    }
  };
  template<typename T, int D>
  struct SIMD<T,D,16> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i   ].X[D],B[i+1 ].X[D],B[i+2 ].X[D],B[i+3 ].X[D],
	  B[i+4 ].X[D],B[i+5 ].X[D],B[i+6 ].X[D],B[i+7 ].X[D],
	  B[i+8 ].X[D],B[i+9 ].X[D],B[i+10].X[D],B[i+11].X[D],
	  B[i+12].X[D],B[i+13].X[D],B[i+14].X[D],B[i+15].X[D]);
      return v;
    }
    static inline T setIndex(int i) {
      T v(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7,i+8,i+9,i+10,i+11,i+12,i+13,i+14,i+15);
      return v;
    }
  };
  template<typename T, int D>
  struct SIMD<T,D,8> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i  ].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D],
	  B[i+4].X[D],B[i+5].X[D],B[i+6].X[D],B[i+7].X[D]);
      return v;
    }
    static inline T setIndex(int i) {
      T v(i,i+1,i+2,i+3,i+4,i+5,i+6,i+7);
      return v;
    }
  };
  template<typename T, int D>
  struct SIMD<T,D,4> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i].X[D],B[i+1].X[D],B[i+2].X[D],B[i+3].X[D]);
      return v;
    }
    static inline T setIndex(int i) {
      T v(i,i+1,i+2,i+3);
      return v;
    }
  };
  template<typename T, int D>
  struct SIMD<T,D,2> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i].X[D],B[i+1].X[D]);
      return v;
    }
    static inline T setIndex(int i) {
      T v(i,i+1);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,3,16> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i   ].SRC,B[i+1 ].SRC,B[i+2 ].SRC,B[i+3 ].SRC,
	  B[i+4 ].SRC,B[i+5 ].SRC,B[i+6 ].SRC,B[i+7 ].SRC,
	  B[i+8 ].SRC,B[i+9 ].SRC,B[i+10].SRC,B[i+11].SRC,
	  B[i+12].SRC,B[i+13].SRC,B[i+14].SRC,B[i+15].SRC);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,3,8> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i  ].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC,
	  B[i+4].SRC,B[i+5].SRC,B[i+6].SRC,B[i+7].SRC);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,3,4> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,3,2> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i].SRC,B[i+1].SRC);
      return v;
    }
  };

  kreal_t transpose(ksimdvec v, int i) {
#if EXAFMM_USE_KAHAN
    kreal_t temp;
    temp.s = v.s[i];
    temp.c = v.c[i];
    return temp;
#else
    return v[i];
#endif
  }
}
#endif
