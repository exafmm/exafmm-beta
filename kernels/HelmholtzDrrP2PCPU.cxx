#include "kernel.h"
#include "simdvec.h"

namespace exafmm {
  namespace kernel {
    real_t eps2;
    complex_t wavek;
    vec3 Xperiodic;

    const complex_t I(0.,1.);

    void P2P(C_iter Ci, C_iter Cj, bool mutual) {
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      simdvec wave_rvec = wave_r;
      simdvec wave_ivec = wave_i;
      for ( ; i<=ni-NSIMD; i+=NSIMD) {
	simdvec zero = 0.0;
	simdvec one = 1.0;
	ksimdvec pot_r = zero;
	ksimdvec pot_i = zero;
	ksimdvec ax_r = zero;
	ksimdvec ax_i = zero;
	ksimdvec ay_r = zero;
	ksimdvec ay_i = zero;
	ksimdvec az_r = zero;
	ksimdvec az_i = zero;

	simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(Bi,i);
	simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(Bi,i);
	simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(Bi,i);
	simdvec mi_r = SIMD<simdvec,4,NSIMD>::setBody(Bi,i);
	simdvec mi_i = SIMD<simdvec,5,NSIMD>::setBody(Bi,i);

	simdvec dx = Xperiodic[0];
	xi -= dx;
	simdvec dy = Xperiodic[1];
	yi -= dy;
	simdvec dz = Xperiodic[2];
	zi -= dz;

	for (int j=0; j<nj; j++) {
	  dx = Bj[j].X[0];
	  dx -= xi;
	  dy = Bj[j].X[1];
	  dy -= yi;
	  dz = Bj[j].X[2];
	  dz -= zi;

	  simdvec R2 = eps2;
	  R2 += dx * dx;
	  simdvec mj_r = std::real(Bj[j].SRC);
	  R2 += dy * dy;
	  simdvec mj_i = std::imag(Bj[j].SRC);
	  R2 += dz * dz;
	  simdvec invR = rsqrt(R2);
	  simdvec R = one / invR;
	  invR &= R2 > zero;
	  R &= R2 > zero;

	  simdvec tmp = mi_r * mj_r - mi_i * mj_i;
	  mj_i = mi_r * mj_i + mi_i * mj_r;
	  mj_r = tmp;
	  tmp = invR / exp(wave_ivec * R);
	  simdvec coef_r = cos(wave_rvec * R) * tmp;
	  simdvec coef_i = sin(wave_rvec * R) * tmp;
	  tmp = mj_r * coef_r - mj_i * coef_i;
	  coef_i = mj_r * coef_i + mj_i * coef_r;
	  coef_r = tmp;
	  mj_r = (one + wave_ivec * R) * invR * invR;
	  mj_i = - wave_rvec * invR;
	  pot_r += coef_r;
	  pot_i += coef_i;
	  if (mutual) Bj[j].TRG[0] += kcomplex_t(sum(coef_r), sum(coef_i));
	  tmp = mj_r * coef_r - mj_i * coef_i;
	  coef_i = mj_r * coef_i + mj_i * coef_r;
	  coef_r = tmp;
	  ax_r += coef_r * dx;
	  ax_i += coef_i * dx;
	  if (mutual) Bj[j].TRG[1] += kcomplex_t(sum(coef_r * dx), sum(coef_i * dx));
	  ay_r += coef_r * dy;
	  ay_i += coef_i * dy;
	  if (mutual) Bj[j].TRG[2] += kcomplex_t(sum(coef_r * dy), sum(coef_i * dy));
	  az_r += coef_r * dz;
	  az_i += coef_i * dz;
	  if (mutual) Bj[j].TRG[3] += kcomplex_t(sum(coef_r * dz), sum(coef_i * dz));
	}
	for (int k=0; k<NSIMD; k++) {
	  Bi[i+k].TRG[0] += transpose(pot_r, pot_i, k);
	  Bi[i+k].TRG[1] -= transpose(ax_r, ax_i, k);
	  Bi[i+k].TRG[2] -= transpose(ay_r, ay_i, k);
	  Bi[i+k].TRG[3] -= transpose(az_r, az_i, k);
	}
      }
#endif
      for ( ; i<ni; i++) {
	real_t pot_r = 0.0;
	real_t pot_i = 0.0;
	real_t ax_r = 0.0;
	real_t ax_i = 0.0;
	real_t ay_r = 0.0;
	real_t ay_i = 0.0;
	real_t az_r = 0.0;
	real_t az_i = 0.0;
	real_t mi_r = std::real(Bi[i].SRC);
	real_t mi_i = std::imag(Bi[i].SRC);
	for (int j=0; j<nj; j++) {
	  real_t mj_r = std::real(Bj[j].SRC);
	  real_t mj_i = std::imag(Bj[j].SRC);
	  vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t R = sqrt(R2);
	    real_t src2_r = mi_r * mj_r - mi_i * mj_i;
	    real_t src2_i = mi_r * mj_i + mi_i * mj_r;
	    real_t expikr = std::exp(wave_i * R) * R;
	    real_t expikr_r = std::cos(wave_r * R) / expikr;
	    real_t expikr_i = std::sin(wave_r * R) / expikr;
	    real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
	    real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
	    real_t kr_r = (1 + wave_i * R) / R2;
	    real_t kr_i = - wave_r / R;
	    real_t coef2_r = kr_r * coef1_r - kr_i * coef1_i;
	    real_t coef2_i = kr_r * coef1_i + kr_i * coef1_r;
	    pot_r += coef1_r;
	    pot_i += coef1_i;
	    ax_r += coef2_r * dX[0];
	    ax_i += coef2_i * dX[0];
	    ay_r += coef2_r * dX[1];
	    ay_i += coef2_i * dX[1];
	    az_r += coef2_r * dX[2];
	    az_i += coef2_i * dX[2];
	    if (mutual) {
	      Bj[j].TRG[0] += complex_t(coef1_r, coef1_i);
	      Bj[j].TRG[1] += complex_t(coef2_r, coef2_i) * dX[0];
	      Bj[j].TRG[2] += complex_t(coef2_r, coef2_i) * dX[1];
	      Bj[j].TRG[3] += complex_t(coef2_r, coef2_i) * dX[2];
	    }
	  }
	}
	Bi[i].TRG[0] += complex_t(pot_r, pot_i);
	Bi[i].TRG[1] += complex_t(ax_r, ax_i);
	Bi[i].TRG[2] += complex_t(ay_r, ay_i);
	Bi[i].TRG[3] += complex_t(az_r, az_i);
      }
    }

    void P2P(C_iter C) {
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      simdvec wave_rvec = wave_r;
      simdvec wave_ivec = wave_i;
      for ( ; i<=n-NSIMD; i+=NSIMD) {
	simdvec zero = 0.0;
	simdvec one = 1.0;
	ksimdvec pot_r = zero;
	ksimdvec pot_i = zero;
	ksimdvec ax_r = zero;
	ksimdvec ax_i = zero;
	ksimdvec ay_r = zero;
	ksimdvec ay_i = zero;
	ksimdvec az_r = zero;
	ksimdvec az_i = zero;

	simdvec index = SIMD<simdvec,0,NSIMD>::setIndex(i);
	simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(B,i);
	simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(B,i);
	simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(B,i);
	simdvec mi_r = SIMD<simdvec,4,NSIMD>::setBody(B,i);
	simdvec mi_i = SIMD<simdvec,5,NSIMD>::setBody(B,i);
	for (int j=i+1; j<n; j++) {
	  simdvec dx = B[j].X[0];
	  dx -= xi;
	  simdvec dy = B[j].X[1];
	  dy -= yi;
	  simdvec dz = B[j].X[2];
	  dz -= zi;

	  simdvec R2 = eps2;
	  R2 += dx * dx;
	  simdvec mj_r = std::real(B[j].SRC);
	  R2 += dy * dy;
	  simdvec mj_i = std::imag(B[j].SRC);
	  R2 += dz * dz;
	  simdvec invR = rsqrt(R2);
	  simdvec R = one / invR;
	  invR &= index < j;
	  invR &= R2 > zero;
	  R &= index < j;
	  R &= R2 > zero;

	  simdvec tmp = mi_r * mj_r - mi_i * mj_i;
	  mj_i = mi_r * mj_i + mi_i * mj_r;
	  mj_r = tmp;
	  tmp = invR / exp(wave_ivec * R);
	  simdvec coef_r = cos(wave_rvec * R) * tmp;
	  simdvec coef_i = sin(wave_rvec * R) * tmp;
	  tmp = mj_r * coef_r - mj_i * coef_i;
	  coef_i = mj_r * coef_i + mj_i * coef_r;
	  coef_r = tmp;
	  mj_r = (one + wave_ivec * R) * invR * invR;
	  mj_i = - wave_rvec * invR;
	  pot_r += coef_r;
	  pot_i += coef_i;
	  B[j].TRG[0] += kcomplex_t(sum(coef_r), sum(coef_i));
	  tmp = mj_r * coef_r - mj_i * coef_i;
	  coef_i = mj_r * coef_i + mj_i * coef_r;
	  coef_r = tmp;
	  ax_r += coef_r * dx;
	  ax_i += coef_i * dx;
	  B[j].TRG[1] += kcomplex_t(sum(coef_r * dx), sum(coef_i * dx));
	  ay_r += coef_r * dy;
	  ay_i += coef_i * dy;
	  B[j].TRG[2] += kcomplex_t(sum(coef_r * dy), sum(coef_i * dy));
	  az_r += coef_r * dz;
	  az_i += coef_i * dz;
	  B[j].TRG[3] += kcomplex_t(sum(coef_r * dz), sum(coef_i * dz));
	}
	for (int k=0; k<NSIMD; k++) {
	  B[i+k].TRG[0] += transpose(pot_r, pot_i, k);
	  B[i+k].TRG[1] -= transpose(ax_r, ax_i, k);
	  B[i+k].TRG[2] -= transpose(ay_r, ay_i, k);
	  B[i+k].TRG[3] -= transpose(az_r, az_i, k);
	}
      }
#endif
      for ( ; i<n; i++) {
	kreal_t pot_r = 0;
	kreal_t pot_i = 0;
	kreal_t ax_r = 0;
	kreal_t ax_i = 0;
	kreal_t ay_r = 0;
	kreal_t ay_i = 0;
	kreal_t az_r = 0;
	kreal_t az_i = 0;
	real_t mi_r = std::real(B[i].SRC);
	real_t mi_i = std::imag(B[i].SRC);
	for (int j=i+1; j<n; j++) {
	  real_t mj_r = std::real(B[j].SRC);
	  real_t mj_i = std::imag(B[j].SRC);
	  vec3 dX = B[j].X - B[i].X;
	  real_t R2 = norm(dX) + eps2;
	  real_t R = sqrt(R2);
	  real_t src2_r = mi_r * mj_r - mi_i * mj_i;
	  real_t src2_i = mi_r * mj_i + mi_i * mj_r;
	  real_t expikr = std::exp(wave_i * R) * R;
	  real_t expikr_r = std::cos(wave_r * R) / expikr;
	  real_t expikr_i = std::sin(wave_r * R) / expikr;
	  real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
	  real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
	  real_t kr_r = (1 + wave_i * R) / R2;
	  real_t kr_i = - wave_r / R;
	  real_t coef2_r = kr_r * coef1_r - kr_i * coef1_i;
	  real_t coef2_i = kr_r * coef1_i + kr_i * coef1_r;
	  pot_r += coef1_r;
	  pot_i += coef1_i;
	  ax_r += coef2_r * dX[0];
	  ax_i += coef2_i * dX[0];
	  ay_r += coef2_r * dX[1];
	  ay_i += coef2_i * dX[1];
	  az_r += coef2_r * dX[2];
	  az_i += coef2_i * dX[2];
	  B[j].TRG[0] += complex_t(coef1_r, coef1_i);
	  B[j].TRG[1] += complex_t(coef2_r, coef2_i) * dX[0];
	  B[j].TRG[2] += complex_t(coef2_r, coef2_i) * dX[1];
	  B[j].TRG[3] += complex_t(coef2_r, coef2_i) * dX[2];
	}
	B[i].TRG[0] += complex_t(pot_r, pot_i);
	B[i].TRG[1] -= complex_t(ax_r, ax_i);
	B[i].TRG[2] -= complex_t(ay_r, ay_i);
	B[i].TRG[3] -= complex_t(az_r, az_i);
      }
    }
  }
}
