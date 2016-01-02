#include "kernel.h"
#if EXAFMM_USE_SIMD
#include "simdvec.h"
#endif

namespace exafmm {
  namespace kernel {
    real_t eps2;
    vec3 Xperiodic;

    void P2P(C_iter Ci, C_iter Cj, bool mutual) {
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      for ( ; i<=ni-NSIMD; i+=NSIMD) {
	simdvec zero = 0.0;
	ksimdvec pot = zero;
	ksimdvec ax = zero;
	ksimdvec ay = zero;
	ksimdvec az = zero;

	simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(Bi,i);
	simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(Bi,i);
	simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(Bi,i);
	simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(Bi,i);

	simdvec xj = Xperiodic[0];
	xi -= xj;
	simdvec yj = Xperiodic[1];
	yi -= yj;
	simdvec zj = Xperiodic[2];
	zi -= zj;

	for (int j=0; j<nj; j++) {
	  simdvec dx = Bj[j].X[0];
	  dx -= xi;
	  simdvec dy = Bj[j].X[1];
	  dy -= yi;
	  simdvec dz = Bj[j].X[2];
	  dz -= zi;
	  simdvec mj = Bj[j].SRC;

	  simdvec R2 = eps2;
	  xj = dx;
	  R2 += dx * dx;
	  yj = dy;
	  R2 += dy * dy;
	  zj = dz;
	  R2 += dz * dz;
	  simdvec invR = rsqrt(R2); 
	  invR &= R2 > zero;

	  mj *= invR * mi;
	  pot += mj;
	  if (mutual) Bj[j].TRG[0] += sum(mj);
	  invR = invR * invR * mj; 

	  xj *= invR;
	  ax += xj;
	  if (mutual) Bj[j].TRG[1] -= sum(xj);

	  yj *= invR;
	  ay += yj;
	  if (mutual) Bj[j].TRG[2] -= sum(yj);

	  zj *= invR;
	  az += zj;
	  if (mutual) Bj[j].TRG[3] -= sum(zj);
	}
	for (int k=0; k<NSIMD; k++) {
	  Bi[i+k].TRG[0] += transpose(pot, k);
	  Bi[i+k].TRG[1] += transpose(ax, k);
	  Bi[i+k].TRG[2] += transpose(ay, k);
	  Bi[i+k].TRG[3] += transpose(az, k);
	}
      }
#endif
      for ( ; i<ni; i++) {
	kreal_t pot = 0; 
	kreal_t ax = 0;
	kreal_t ay = 0;
	kreal_t az = 0;
	for (int j=0; j<nj; j++) {
	  vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t invR2 = 1.0 / R2;
	    real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
	    dX *= invR2 * invR;
	    pot += invR;
	    ax += dX[0];
	    ay += dX[1];
	    az += dX[2];
	    if (mutual) {
	      Bj[j].TRG[0] += invR;
	      Bj[j].TRG[1] += dX[0];
	      Bj[j].TRG[2] += dX[1];
	      Bj[j].TRG[3] += dX[2];
	    }
	  }
	}
	Bi[i].TRG[0] += pot;
	Bi[i].TRG[1] -= ax;
	Bi[i].TRG[2] -= ay;
	Bi[i].TRG[3] -= az;
      }
    }

    void P2P(C_iter C) {
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      for ( ; i<=n-NSIMD; i+=NSIMD) {
	simdvec zero = 0;
	ksimdvec pot = zero;
	ksimdvec ax = zero;
	ksimdvec ay = zero;
	ksimdvec az = zero;

	simdvec index = SIMD<simdvec,0,NSIMD>::setIndex(i);
	simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(B,i);
	simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(B,i);
	simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(B,i);
	simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(B,i);
	for (int j=i+1; j<n; j++) {
	  simdvec dx = B[j].X[0];
	  dx -= xi;
	  simdvec dy = B[j].X[1];
	  dy -= yi;
	  simdvec dz = B[j].X[2];
	  dz -= zi;
	  simdvec mj = B[j].SRC;

	  simdvec R2 = eps2;
	  simdvec xj = dx;
	  R2 += dx * dx;
	  simdvec yj = dy;
	  R2 += dy * dy;
	  simdvec zj = dz;
	  R2 += dz * dz;
	  simdvec invR = rsqrt(R2);
	  invR &= index < j;
	  invR &= R2 > zero;

	  mj *= invR * mi;
	  pot += mj;
	  B[j].TRG[0] += sum(mj);
	  invR = invR * invR * mj;

	  xj *= invR;
	  ax += xj;
	  B[j].TRG[1] -= sum(xj);

	  yj *= invR;
	  ay += yj;
	  B[j].TRG[2] -= sum(yj);

	  zj *= invR;
	  az += zj;
	  B[j].TRG[3] -= sum(zj);
	}
	for (int k=0; k<NSIMD; k++) {
	  B[i+k].TRG[0] += transpose(pot, k);
	  B[i+k].TRG[1] += transpose(ax, k);
	  B[i+k].TRG[2] += transpose(ay, k);
	  B[i+k].TRG[3] += transpose(az, k);
	}
      }
#endif
      for ( ; i<n; i++) {
	kreal_t pot = 0;
	kreal_t ax = 0;
	kreal_t ay = 0;
	kreal_t az = 0;
	for (int j=i+1; j<n; j++) {
	  vec3 dX = B[j].X - B[i].X;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t invR2 = 1.0 / R2;
	    real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
	    dX *= invR2 * invR;
	    pot += invR;
	    ax += dX[0];
	    ay += dX[1];
	    az += dX[2];
	    B[j].TRG[0] += invR;
	    B[j].TRG[1] -= dX[0];
	    B[j].TRG[2] -= dX[1];
	    B[j].TRG[3] -= dX[2];
	  }
	}
	B[i].TRG[0] += pot;
	B[i].TRG[1] += ax;
	B[i].TRG[2] += ay;
	B[i].TRG[3] += az;
      }
    }
  }
}
