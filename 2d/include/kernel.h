#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "log2f4.h"
#include "types.h"

class Kernel {
 private:
  real_t EPS2;                                                  //!< Softening parameter (squared)

 protected:
  vec2 Xperiodic;                                               //!< Coordinate offset for periodic B.C.

 public:
//!< Constructor
  Kernel(real_t _EPS2) : EPS2(_EPS2), Xperiodic(0) {}

//!< P2P kernel between cells Ci and Cj 
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const {
    B_iter Bi = Ci->BODY;
    B_iter Bj = Cj->BODY;
    int ni = Ci->NDBODY;
    int nj = Cj->NDBODY;
    int i = 0;
#if USE_SIMD
    __m128 zero = _mm_set1_ps(0.0);
    __m128 one = _mm_set1_ps(1.0);
    __m128 base = _mm_set1_ps(M_E);
    base = log2f4(base);
    base = _mm_div_ps(one, base);
    for ( ; i<=ni-4; i+=4) {
      __m128 pot = zero;
      __m128 xi = _mm_setr_ps(Bi[i].X[0], Bi[i+1].X[0], Bi[i+2].X[0], Bi[i+3].X[0]);
      __m128 yi = _mm_setr_ps(Bi[i].X[1], Bi[i+1].X[1], Bi[i+2].X[1], Bi[i+3].X[1]);
      __m128 mi = _mm_setr_ps(Bi[i].SRC,  Bi[i+1].SRC,  Bi[i+2].SRC,  Bi[i+3].SRC);
      __m128 R2 = _mm_set1_ps(EPS2);

      __m128 xj = _mm_set1_ps(Xperiodic[0]);
      xi = _mm_sub_ps(xi, xj);
      __m128 yj = _mm_set1_ps(Xperiodic[1]);
      yi = _mm_sub_ps(yi, yj);

      __m128 x2 = _mm_set1_ps(Bj[0].X[0]);
      x2 = _mm_sub_ps(x2, xi);
      __m128 y2 = _mm_set1_ps(Bj[0].X[1]);
      y2 = _mm_sub_ps(y2, yi);
      __m128 mj = _mm_set1_ps(Bj[0].SRC);

      xj = x2;
      x2 = _mm_mul_ps(x2, x2);
      R2 = _mm_add_ps(R2, x2);
      yj = y2;
      y2 = _mm_mul_ps(y2, y2);
      R2 = _mm_add_ps(R2, y2);
      __m128 logR;

      x2 = _mm_set1_ps(Bj[1].X[0]);
      y2 = _mm_set1_ps(Bj[1].X[1]);
      for (int j=0; j<nj-2; j++) {
	logR = _mm_rsqrt_ps(R2);
        logR = log2f4(logR);
        R2 = _mm_cmpgt_ps(R2, zero);
        logR = _mm_and_ps(logR, R2);
	R2 = _mm_set1_ps(EPS2);
	x2 = _mm_sub_ps(x2, xi);
	y2 = _mm_sub_ps(y2, yi);

        logR = _mm_mul_ps(logR, base);
        logR = _mm_mul_ps(logR, mi);
        mj = _mm_mul_ps(mj, logR);
        pot = _mm_add_ps(pot, mj);
        mj = _mm_hadd_ps(mj, mj);
        mj = _mm_hadd_ps(mj, mj);
	if (mutual) Bj[j].TRG += ((float*)&mj)[0];
	mj = _mm_set1_ps(Bj[j+1].SRC);

	xj = x2;
	x2 = _mm_mul_ps(x2, x2);
	R2 = _mm_add_ps(R2, x2);
	x2 = _mm_set1_ps(Bj[j+2].X[0]);
	yj = y2;
	y2 = _mm_mul_ps(y2, y2);
	R2 = _mm_add_ps(R2, y2);
	y2 = _mm_set1_ps(Bj[j+2].X[1]);
      }
      if ( nj > 1 ) {
	logR = _mm_rsqrt_ps(R2);
        logR = log2f4(logR);
        R2 = _mm_cmpgt_ps(R2, zero);
        logR = _mm_and_ps(logR, R2);
        R2 = _mm_set1_ps(EPS2);
        x2 = _mm_sub_ps(x2, xi);
        y2 = _mm_sub_ps(y2, yi);

	logR = _mm_mul_ps(logR, base);
        logR = _mm_mul_ps(logR, mi);
	mj = _mm_mul_ps(mj, logR);
        pot = _mm_add_ps(pot, mj);
        mj = _mm_hadd_ps(mj, mj);
        mj = _mm_hadd_ps(mj, mj);
        if (mutual) Bj[nj-2].TRG += ((float*)&mj)[0];
        mj = _mm_set1_ps(Bj[nj-1].SRC);

	xj = x2;
        x2 = _mm_mul_ps(x2, x2);
        R2 = _mm_add_ps(R2, x2);
        yj = y2;
        y2 = _mm_mul_ps(y2, y2);
	R2 = _mm_add_ps(R2, y2);
      }
      logR = _mm_rsqrt_ps(R2);
      logR = log2f4(logR);
      R2 = _mm_cmpgt_ps(R2, zero);
      logR = _mm_and_ps(logR, R2);
      R2 = _mm_set1_ps(EPS2);
      x2 = _mm_sub_ps(x2, xi);
      y2 = _mm_sub_ps(y2, yi);

      logR = _mm_mul_ps(logR, base);
      logR = _mm_mul_ps(logR, mi);
      mj = _mm_mul_ps(mj, logR);
      pot = _mm_add_ps(pot, mj);

      mj = _mm_hadd_ps(mj, mj);
      mj = _mm_hadd_ps(mj, mj);
      if (mutual) Bj[nj-1].TRG += ((float*)&mj)[0];

      for (int k=0; k<4; k++) {
	Bi[i+k].TRG += ((float*)&pot)[k];
      }
    }
#endif
    for (; i<ni; i++) {
      real_t pot = 0;
      for (int j=0; j<nj; j++) {
	vec2 dX = Bi[i].X - Bj[j].X - Xperiodic;
	real_t R2 = norm(dX) + EPS2;
	if (R2 != 0) {
	  real_t invR = 1 / sqrt(R2);
	  real_t logR = Bi[i].SRC * Bj[j].SRC * log(invR);
	  pot += logR;
	  if (mutual) {
	    Bj[j].TRG += logR;
	  }
	}
      }
      Bi[i].TRG += pot;
    }
  }

//!< P2P kernel for cell C
  void P2P(C_iter C) const {
    B_iter B = C->BODY;
    int n = C->NDBODY;
    int i = 0;
#if USE_SIMD
    __m128 zero = _mm_set1_ps(0.0);
    __m128 one = _mm_set1_ps(1.0);
    __m128 base = _mm_set1_ps(M_E);
    base = log2f4(base);
    base = _mm_div_ps(one, base);
    for ( ; i<=n-4; i+=4) {
      __m128 pot = zero;
      __m128 index = _mm_setr_ps(i, i+1, i+2, i+3);
      __m128 xi = _mm_setr_ps(B[i].X[0], B[i+1].X[0], B[i+2].X[0], B[i+3].X[0]);
      __m128 yi = _mm_setr_ps(B[i].X[1], B[i+1].X[1], B[i+2].X[1], B[i+3].X[1]);
      __m128 mi = _mm_setr_ps(B[i].SRC,  B[i+1].SRC,  B[i+2].SRC,  B[i+3].SRC);
      __m128 R2 = _mm_set1_ps(EPS2);

      __m128 xj = _mm_set1_ps(Xperiodic[0]);
      xi = _mm_sub_ps(xi, xj);
      __m128 yj = _mm_set1_ps(Xperiodic[1]);
      yi = _mm_sub_ps(yi, yj);

      __m128 x2 = _mm_set1_ps(B[i+1].X[0]);
      x2 = _mm_sub_ps(x2, xi);
      __m128 y2 = _mm_set1_ps(B[i+1].X[1]);
      y2 = _mm_sub_ps(y2, yi);
      __m128 mj = _mm_set1_ps(B[i+1].SRC);

      xj = x2;
      x2 = _mm_mul_ps(x2, x2);
      R2 = _mm_add_ps(R2, x2);
      yj = y2;
      y2 = _mm_mul_ps(y2, y2);
      R2 = _mm_add_ps(R2, y2);
      __m128 logR, flag;

      x2 = _mm_set1_ps(B[i+2].X[0]);
      y2 = _mm_set1_ps(B[i+2].X[1]);
      for (int j=i+1; j<n-2; j++) {
	logR = _mm_rsqrt_ps(R2);
        logR = log2f4(logR);
        flag = _mm_cmplt_ps(index, _mm_set1_ps(j));
        logR = _mm_and_ps(logR, flag);
        flag = _mm_cmpgt_ps(R2, zero);
        logR = _mm_and_ps(logR, flag);
	R2 = _mm_set1_ps(EPS2);
	x2 = _mm_sub_ps(x2, xi);
	y2 = _mm_sub_ps(y2, yi);

        logR = _mm_mul_ps(logR, base);
        logR = _mm_mul_ps(logR, mi);
        mj = _mm_mul_ps(mj, logR);
        pot = _mm_add_ps(pot, mj);
        mj = _mm_hadd_ps(mj, mj);
        mj = _mm_hadd_ps(mj, mj);
	B[j].TRG += ((float*)&mj)[0];
	mj = _mm_set1_ps(B[j+1].SRC);

	xj = x2;
	x2 = _mm_mul_ps(x2, x2);
	R2 = _mm_add_ps(R2, x2);
	x2 = _mm_set1_ps(B[j+2].X[0]);
	yj = y2;
	y2 = _mm_mul_ps(y2, y2);
	R2 = _mm_add_ps(R2, y2);
	y2 = _mm_set1_ps(B[j+2].X[1]);
      }
      if ( n-2 > i ) {
	logR = _mm_rsqrt_ps(R2);
        logR = log2f4(logR);
	flag = _mm_cmplt_ps(index, _mm_set1_ps(n-2));
        logR = _mm_and_ps(logR, flag);
        flag = _mm_cmpgt_ps(R2, zero);
        logR = _mm_and_ps(logR, flag);
        R2 = _mm_set1_ps(EPS2);
        x2 = _mm_sub_ps(x2, xi);
        y2 = _mm_sub_ps(y2, yi);

	logR = _mm_mul_ps(logR, base);
        logR = _mm_mul_ps(logR, mi);
	mj = _mm_mul_ps(mj, logR);
        pot = _mm_add_ps(pot, mj);
        mj = _mm_hadd_ps(mj, mj);
        mj = _mm_hadd_ps(mj, mj);
        B[n-2].TRG += ((float*)&mj)[0];
        mj = _mm_set1_ps(B[n-1].SRC);

	xj = x2;
        x2 = _mm_mul_ps(x2, x2);
        R2 = _mm_add_ps(R2, x2);
        yj = y2;
        y2 = _mm_mul_ps(y2, y2);
	R2 = _mm_add_ps(R2, y2);
      }
      logR = _mm_rsqrt_ps(R2);
      logR = log2f4(logR);
      flag = _mm_cmplt_ps(index, _mm_set1_ps(n-1));
      logR = _mm_and_ps(logR, flag);
      flag = _mm_cmpgt_ps(R2, zero);
      logR = _mm_and_ps(logR, flag);
      R2 = _mm_set1_ps(EPS2);
      x2 = _mm_sub_ps(x2, xi);
      y2 = _mm_sub_ps(y2, yi);

      logR = _mm_mul_ps(logR, base);
      logR = _mm_mul_ps(logR, mi);
      mj = _mm_mul_ps(mj, logR);
      pot = _mm_add_ps(pot, mj);

      mj = _mm_hadd_ps(mj, mj);
      mj = _mm_hadd_ps(mj, mj);
      B[n-1].TRG += ((float*)&mj)[0];

      for (int k=0; k<4; k++) {
	B[i+k].TRG += ((float*)&pot)[k];
      }
    }
#endif
    for (; i<n; i++) {
      real_t pot = 0;
      for (int j=i+1; j<n; j++) {
	vec2 dX = B[i].X - B[j].X;
	real_t R2 = norm(dX) + EPS2;
	if (R2 != 0) {
	  real_t invR = 1 / sqrt(R2);
	  real_t logR = B[i].SRC * B[j].SRC * log(invR);
	  pot += logR;
	  B[j].TRG += logR;
	}
      }
      B[i].TRG += pot;
    }
  }

//!< P2M kernel for cell C
  void P2M(C_iter C) const {
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {
      vec2 dX = B->X - C->X;
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);
      vecP M;
      C->M[0] += B->SRC;
      for (int n=1; n<P; n++) {
        powZ *= Z;
        C->M[n] += -powZ * real_t(B->SRC / n);
      }
    }
  }

//!< M2M kernel for one parent cell Ci
  void M2M(C_iter Ci, C_iter C0) const {
    for (C_iter Cj=C0+Ci->CHILD; Cj!=C0+Ci->CHILD+Ci->NCHILD; Cj++) {
      vec2 dX = Cj->X - Ci->X;
      complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0), powZnk(1.0, 0.0), invZ(powZn/Z);
      Ci->M[0] += Cj->M[0];
      for (int n=1; n<P; n++) {
        powZn *= Z;
	Ci->M[n] -= Cj->M[0] * powZn / real_t(n);
        powZnk = powZn;
        real_t Cnk = 1.0;
	for (int k=1; k<=n; k++) {
          powZnk *= invZ;
	  Ci->M[n] += Cj->M[k] * powZnk * Cnk;
          Cnk *= real_t(n - k) / k;
	}
      }
    }
  }

//!< M2L kernel between cells Ci and Cj
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const {
    vec2 dX = Ci->X - Cj->X - Xperiodic;
    complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0), powZnk(1.0, 0.0), invZ(powZn/Z);
    Ci->L[0] += Cj->M[0] * log(Z);
    for (int k=1; k<P; k++) {
      powZn *= invZ;
      Ci->L[0] += Cj->M[k] * powZn;
    }
    powZn = complex_t(1.0, 0.0);
    for (int n=1; n<P; n++) {
      powZn *= -invZ;
      Ci->L[n] -= Cj->M[0] * powZn / real_t(n);
      powZnk = powZn;
      real_t Cnk = 1.0;
      for (int k=1; k<P; k++) {
        powZnk *= invZ;
	Ci->L[n] += Cj->M[k] * powZnk * Cnk;
        Cnk *= real_t(n + 1) / k;
      }
    }
    if (mutual) {
      Cj->L[0] += Ci->M[0] * log(-Z);
      powZn = complex_t(1.0, 0.0);
      for (int k=1; k<P; k++) {
        powZn *= -invZ;
	Cj->L[0] += Ci->M[k] * powZn;
      }
      powZn = complex_t(1.0, 0.0);
      for (int n=1; n<P; n++) {
        powZn *= invZ;
	Cj->L[n] -= Ci->M[0] * powZn / real_t(n);
        powZnk = powZn;
        real_t Cnk = 1.0;
	for (int k=1; k<P; k++) {
          powZnk *= -invZ;
	  Cj->L[n] += Ci->M[k] * powZnk * Cnk;
          Cnk *= real_t(n + 1) / k;
	}
      }
    }
  }

//!< L2L kernel for one child cell Ci
  void L2L(C_iter Ci, C_iter C0) const {
    C_iter Cj = C0 + Ci->PARENT;
    vec2 dX = Ci->X - Cj->X;
    complex_t Z(dX[0],dX[1]);
    for (int n=0; n<P; n++) {
      complex_t powZ(1.0, 0.0);
      real_t Cnk = 1.0;
      for (int k=n; k<P; k++) {
	Ci->L[n] += Cj->L[k] * powZ * Cnk;
        powZ *= Z;
        Cnk *= real_t(k + 1) / (k - n + 1);
      }
    }
  }

//!< L2P kernel for cell Ci
  void L2P(C_iter Ci) const {
    for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {
      vec2 dX = B->X - Ci->X;
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);
      B->TRG /= B->SRC;
      for (int n=0; n<P; n++) {
        B->TRG -= std::real(Ci->L[n] * powZ);
        powZ *= Z;
      }
    }
  }
};
#endif
