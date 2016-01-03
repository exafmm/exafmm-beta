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
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {         // Loop over bodies
      vec2 dX = B->X - C->X;                                    //  Get distance vector
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);                 //  Convert to complex plane
      C->M[0] += B->SRC;                                        //  Add constant term
      for (int n=1; n<P; n++) {                                 //  Loop over coefficients
        powZ *= Z / real_t(n);                                  //   Store z^n / n!
        C->M[n] += powZ * B->SRC;                               //   Add to coefficient
      }                                                         //  End loop
    }                                                           // End loop
  }

//!< M2M kernel for one parent cell Ci
  void M2M(C_iter Ci, C_iter C0) const {
    for (C_iter Cj=C0+Ci->CHILD;
	 Cj!=C0+Ci->CHILD+Ci->NCHILD; Cj++) {                   // Loop over child cells
      vec2 dX = Cj->X - Ci->X;                                  //  Get distance vector
      complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0),
	powZnk(1.0, 0.0), invZ(powZn/Z);                        //  Convert to complex plane
      for (int k=0; k<P; k++) {                                 //  Loop over coefficients
	complex_t powZ(1.0, 0.0);                               //   z^0 = 1
	Ci->M[k] += Cj->M[k];                                   //   Add constant term
	for (int kml=1; kml<=k; kml++) {                        //   Loop over k-l
	  powZ *= Z / real_t(kml);                              //    Store z^(k-l) / (k-l)!
	  Ci->M[k] += Cj->M[k-kml] * powZ;                      //    Add to coefficient
	}                                                       //   End loop
      }                                                         //  End loop
    }                                                           // End loop
  }

//!< M2L kernel between cells Ci and Cj
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const {
    vec2 dX = Ci->X - Cj->X - Xperiodic;                        // Get distance vector
    complex_t Z(dX[0],dX[1]), powZn(1.0, 0.0),
      powZnk(1.0, 0.0), invZ(powZn/Z);                          // Convert to complex plane
    Ci->L[0] += -Cj->M[0] * log(Z);                             // Log term (for 0th order)
    Ci->L[0] += Cj->M[1] * invZ;                                // Constant term
    powZn = invZ;                                               // 1/z term
    for (int k=2; k<P; k++) {                                   // Loop over coefficients
      powZn *= real_t(k-1)*invZ;                                //  Store (k-1)! / z^k
      Ci->L[0] += Cj->M[k] * powZn;                             //  Add to coefficient
    }                                                           // End loop
    Ci->L[1] += -Cj->M[0] * invZ;                               // Constant term (for 1st order)
    powZn = invZ;                                               // 1/z term
    for (int k=1; k<P; k++) {                                   // Loop over coefficients
      powZn *= real_t(1+k-1)*invZ;                              //  Store (k)! / z^k
      Ci->L[1] += -Cj->M[k] * powZn;                            //  Add to coefficient
    }                                                           // End loop
    real_t Cnk = -1;                                            // Fix sign term
    for (int n=2; n<P; n++) {                                   // Loop over 
      Cnk *= -1;                                                //  Flip sign
      powZnk *= invZ;                                           //  Store 1 / z^n
      powZn = Cnk*powZnk;                                       //  Combine terms
      for (int k=0; k<P; k++) {                                 //  Loop over
	powZn *= real_t(n+k-1)*invZ;                            //   (n+k-1)! / z^k
	Ci->L[n] += Cj->M[k] * powZn;                           //   Add to coefficient
      }                                                         //  End loop
      powZnk *= real_t(n-1);                                    //  Store (n-1)! / z^n
    }                                                           // End loop
    if (mutual) {                                               // If mutual interaction
      M2L(Cj,Ci,false);                                         //  Lazy recursion
    }                                                           // End if
  }

//!< L2L kernel for one child cell Ci
  void L2L(C_iter Ci, C_iter C0) const {
    C_iter Cj = C0 + Ci->PARENT;                                // Get parent cell
    vec2 dX = Ci->X - Cj->X;                                    // Get distance vector
    complex_t Z(dX[0],dX[1]);                                   // Convert to complex plane
    for (int l=0; l<P; l++) {                                   // Loop over coefficients
      complex_t powZ(1.0, 0.0);                                 //  z^0 = 1
      Ci->L[l] += Cj->L[l];                                     //  Add constant term
      for (int k=1; k<P-l; k++) {                               //  Loop over coefficients
	powZ *= Z / real_t(k);                                  //   Store z^k / k!
	Ci->L[l] += Cj->L[l+k] * powZ;                          //   Add to coefficient
      }                                                         //  End loop
    }                                                           // End loop
  }

//!< L2P kernel for cell Ci
  void L2P(C_iter Ci) const {
    for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {      // Loop over bodies
      vec2 dX = B->X - Ci->X;                                   //  Get distance vector
      complex_t Z(dX[0],dX[1]), powZ(1.0, 0.0);                 //  Convert to complex plane
      B->TRG /= B->SRC;                                         //  Normalize result
      B->TRG += std::real(Ci->L[0]);                            //  Add constant term
      for (int n=1; n<P; n++) {                                 //  Loop over coefficients
        powZ *= Z / real_t(n);                                  //   Store z^n / n!
        B->TRG += std::real(Ci->L[n] * powZ);                   //   Add real part to solution
      }                                                         //  End loop
    }                                                           // End loop
  }
};
#endif
