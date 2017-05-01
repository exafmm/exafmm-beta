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
#if EXAFMM_HELMHOLTZ
  template<typename T>
  struct SIMD<T,4,16> {
    static inline T setBody(B_iter B, int i) {      
#if EXAFMM_ACOUSTICS      
    T v(std::real(B[i   ].SRC * B[i].QWEIGHT),std::real(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
	  std::real(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::real(B[i+3 ].SRC * B[i+3 ].QWEIGHT),
	  std::real(B[i+4 ].SRC * B[i+4 ].QWEIGHT),std::real(B[i+5 ].SRC * B[i+5 ].QWEIGHT),
	  std::real(B[i+6 ].SRC * B[i+6 ].QWEIGHT),std::real(B[i+7 ].SRC * B[i+7 ].QWEIGHT),
	  std::real(B[i+8 ].SRC * B[i+8 ].QWEIGHT),std::real(B[i+9 ].SRC * B[i+9 ].QWEIGHT),
	  std::real(B[i+10].SRC * B[i+10 ].QWEIGHT),std::real(B[i+11].SRC * B[i+11 ].QWEIGHT),
	  std::real(B[i+12].SRC * B[i+12 ].QWEIGHT),std::real(B[i+13].SRC * B[i+13 ].QWEIGHT),
	  std::real(B[i+14].SRC * B[i+14 ].QWEIGHT),std::real(B[i+15].SRC * B[i+15 ].QWEIGHT));
#else 
      T v(std::real(B[i   ].SRC),std::real(B[i+1 ].SRC),
    std::real(B[i+2 ].SRC),std::real(B[i+3 ].SRC),
    std::real(B[i+4 ].SRC),std::real(B[i+5 ].SRC),
    std::real(B[i+6 ].SRC),std::real(B[i+7 ].SRC),
    std::real(B[i+8 ].SRC),std::real(B[i+9 ].SRC),
    std::real(B[i+10].SRC),std::real(B[i+11].SRC),
    std::real(B[i+12].SRC),std::real(B[i+13].SRC),
    std::real(B[i+14].SRC),std::real(B[i+15].SRC));
#endif      
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,4,8> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS
        T v(std::real(B[i   ].SRC * B[i].QWEIGHT),std::real(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
    std::real(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::real(B[i+3 ].SRC * B[i+3 ].QWEIGHT),
    std::real(B[i+4 ].SRC * B[i+4 ].QWEIGHT),std::real(B[i+5 ].SRC * B[i+5 ].QWEIGHT),
    std::real(B[i+6 ].SRC * B[i+6 ].QWEIGHT),std::real(B[i+7 ].SRC * B[i+7 ].QWEIGHT));
#else      
      T v(std::real(B[i  ].SRC),std::real(B[i+1].SRC),
	  std::real(B[i+2].SRC),std::real(B[i+3].SRC),
	  std::real(B[i+4].SRC),std::real(B[i+5].SRC),
	  std::real(B[i+6].SRC),std::real(B[i+7].SRC));
#endif
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,4,4> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS
         T v(std::real(B[i   ].SRC * B[i].QWEIGHT),std::real(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
    std::real(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::real(B[i+3 ].SRC * B[i+3 ].QWEIGHT));     
#else      
      T v(std::real(B[i  ].SRC),std::real(B[i+1].SRC),
	  std::real(B[i+2].SRC),std::real(B[i+3].SRC));
#endif      
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,4,2> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS      
    T v(std::real(B[i   ].SRC * B[i].QWEIGHT),std::real(B[i+1 ].SRC * B[i+1 ].QWEIGHT));
#else
    T v(std::real(B[i].SRC),std::real(B[i+1].SRC));      
#endif      
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,5,16> {
    static inline T setBody(B_iter B, int i) {      
#if EXAFMM_ACOUSTICS      
    T v(std::imag(B[i   ].SRC * B[i].QWEIGHT),std::imag(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
    std::imag(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::imag(B[i+3 ].SRC * B[i+3 ].QWEIGHT),
    std::imag(B[i+4 ].SRC * B[i+4 ].QWEIGHT),std::imag(B[i+5 ].SRC * B[i+5 ].QWEIGHT),
    std::imag(B[i+6 ].SRC * B[i+6 ].QWEIGHT),std::imag(B[i+7 ].SRC * B[i+7 ].QWEIGHT),
    std::imag(B[i+8 ].SRC * B[i+8 ].QWEIGHT),std::imag(B[i+9 ].SRC * B[i+9 ].QWEIGHT),
    std::imag(B[i+10].SRC * B[i+10 ].QWEIGHT),std::imag(B[i+11].SRC * B[i+11 ].QWEIGHT),
    std::imag(B[i+12].SRC * B[i+12 ].QWEIGHT),std::imag(B[i+13].SRC * B[i+13 ].QWEIGHT),
    std::imag(B[i+14].SRC * B[i+14 ].QWEIGHT),std::imag(B[i+15].SRC * B[i+15 ].QWEIGHT));
#else 
      T v(std::imag(B[i   ].SRC),std::imag(B[i+1 ].SRC),
    std::imag(B[i+2 ].SRC),std::imag(B[i+3 ].SRC),
    std::imag(B[i+4 ].SRC),std::imag(B[i+5 ].SRC),
    std::imag(B[i+6 ].SRC),std::imag(B[i+7 ].SRC),
    std::imag(B[i+8 ].SRC),std::imag(B[i+9 ].SRC),
    std::imag(B[i+10].SRC),std::imag(B[i+11].SRC),
    std::imag(B[i+12].SRC),std::imag(B[i+13].SRC),
    std::imag(B[i+14].SRC),std::imag(B[i+15].SRC));
#endif      
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,5,8> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS
        T v(std::imag(B[i   ].SRC * B[i].QWEIGHT),std::imag(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
    std::imag(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::imag(B[i+3 ].SRC * B[i+3 ].QWEIGHT),
    std::imag(B[i+4 ].SRC * B[i+4 ].QWEIGHT),std::imag(B[i+5 ].SRC * B[i+5 ].QWEIGHT),
    std::imag(B[i+6 ].SRC * B[i+6 ].QWEIGHT),std::imag(B[i+7 ].SRC * B[i+7 ].QWEIGHT));
#else      
      T v(std::imag(B[i  ].SRC),std::imag(B[i+1].SRC),
    std::imag(B[i+2].SRC),std::imag(B[i+3].SRC),
    std::imag(B[i+4].SRC),std::imag(B[i+5].SRC),
    std::imag(B[i+6].SRC),std::imag(B[i+7].SRC));
#endif
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,5,4> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS
         T v(std::imag(B[i   ].SRC * B[i].QWEIGHT),std::imag(B[i+1 ].SRC * B[i+1 ].QWEIGHT),
    std::imag(B[i+2 ].SRC * B[i+2 ].QWEIGHT),std::imag(B[i+3 ].SRC * B[i+3 ].QWEIGHT));     
#else      
      T v(std::imag(B[i  ].SRC),std::imag(B[i+1].SRC),
    std::imag(B[i+2].SRC),std::imag(B[i+3].SRC));
#endif      
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,5,2> {
    static inline T setBody(B_iter B, int i) {
#if EXAFMM_ACOUSTICS      
    T v(std::imag(B[i   ].SRC * B[i].QWEIGHT),std::imag(B[i+1 ].SRC * B[i+1 ].QWEIGHT));
#else
    T v(std::imag(B[i].SRC),std::imag(B[i+1].SRC));      
#endif      
      return v;
    }
  };

#if EXAFMM_ACOUSTICS
  template<typename T>
  struct SIMD<T,6,16> {
    static inline T setBody(B_iter B, int i) {      
      T v(B[i   ].PATCH,B[i+1 ].PATCH,
          B[i+2 ].PATCH,B[i+3 ].PATCH,
          B[i+4 ].PATCH,B[i+5 ].PATCH,
          B[i+6 ].PATCH,B[i+7 ].PATCH,
          B[i+8 ].PATCH,B[i+9 ].PATCH,
          B[i+10].PATCH,B[i+11].PATCH,
          B[i+12].PATCH,B[i+13].PATCH,
          B[i+14].PATCH,B[i+15].PATCH);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,6,8> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i   ].PATCH,B[i+1 ].PATCH,
        B[i+2 ].PATCH,B[i+3 ].PATCH,
        B[i+4 ].PATCH,B[i+5 ].PATCH,
        B[i+6 ].PATCH,B[i+7 ].PATCH);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,6,4> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i   ].PATCH,B[i+1 ].PATCH,
          B[i+2 ].PATCH,B[i+3 ].PATCH);
      return v;
    }
  };
  template<typename T>
  struct SIMD<T,6,2> {
    static inline T setBody(B_iter B, int i) {
      T v(B[i   ].PATCH,B[i+1 ].PATCH);      
      return v;
    }
  };
#endif

#endif

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

  kcomplex_t transpose(ksimdvec v_r, ksimdvec v_i, int i) {
#if EXAFMM_USE_KAHAN
    kcomplex_t temp;
    temp.s = complex_t(v_r.s[i], v_i.s[i]);
    temp.c = complex_t(v_r.c[i], v_i.c[i]);
    return temp;
#else
    return kcomplex_t(v_r[i], v_i[i]);
#endif
  }
}
#endif
