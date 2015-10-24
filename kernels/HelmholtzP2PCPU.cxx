#include "kernel.h"
#include "simdvec.h"
using namespace exafmm;

real_t kernel::eps2;
complex_t kernel::wavek;
vec3 kernel::Xperiodic;

const complex_t I(0.,1.);

void kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) {
  assert(mutual == 0);
  real_t wave_r = std::real(wavek);
  real_t wave_i = std::imag(wavek);
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NBODY;
  int nj = Cj->NBODY;
  int i = 0;
  simdvec wave_rvec = wave_r;
  simdvec wave_ivec = wave_i;
#if EXAFMM_USE_SIMD
#if 0
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

      simdvec R2 = eps2;
      R2 += dx * dx;
      simdvec mj_r = std::real(Bj[j].SRC);
      R2 += dy * dy;
      simdvec mj_i = std::imag(Bj[j].SRC);
      R2 += dz * dz;
      simdvec invR = rsqrt(R2);
      invR &= R2 > zero;
      simdvec R = one / invR;

      simdvec tmp = mi_r * mj_r - mi_i * mj_i;
      mj_i = mi_r * mj_i + mi_i * mj_r;
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
      tmp = mj_r * coef_r - mj_i * coef_i;
      coef_i = mj_r * coef_i + mj_i * coef_r;
      coef_r = tmp;
      ax_r += coef_r * dx;
      ax_i += coef_i * dx;
      ay_r += coef_r * dy;
      ay_i += coef_i * dy;
      az_r += coef_r * dz;
      az_i += coef_i * dz;
    }
    for (int k=0; k<NSIMD; k++) {
      Bi[i+k].TRG[0] += transpose(pot_r, pot_i, k);
      Bi[i+k].TRG[1] -= transpose(ax_r, ax_i, k);
      Bi[i+k].TRG[2] -= transpose(ay_r, ay_i, k);
      Bi[i+k].TRG[3] -= transpose(az_r, az_i, k);
    }
  }
#else
  for ( ; i<=ni-NSIMD; i++) {
    real_t pot_r = 0.0;
    real_t pot_i = 0.0;
    real_t ax_r = 0.0;
    real_t ax_i = 0.0;
    real_t ay_r = 0.0;
    real_t ay_i = 0.0;
    real_t az_r = 0.0;
    real_t az_i = 0.0;
    real_t xi = Bi[i].X[0] - Xperiodic[0];
    real_t yi = Bi[i].X[1] - Xperiodic[1];
    real_t zi = Bi[i].X[2] - Xperiodic[2];
    real_t mi_r = std::real(Bi[i].SRC);
    real_t mi_i = std::imag(Bi[i].SRC);
    for (int j=0; j<nj; j++) {
      real_t R2 = eps2;
      real_t dx = Bj[j].X[0];
      dx -= xi;
      real_t dy = Bj[j].X[1];
      dy -= yi;
      real_t dz = Bj[j].X[2];
      dz -= zi;
      real_t mj_r = std::real(Bj[j].SRC);
      R2 += dx * dx;
      real_t mj_i = std::imag(Bj[j].SRC);
      R2 += dy * dy;
      R2 += dz * dz;
      if (R2 != 0) {
	real_t invR = 1 / sqrt(R2);
	real_t R = 1 / invR;
	real_t tmp = mi_r * mj_r - mi_i * mj_i;
	mj_i = mi_r * mj_i + mi_i * mj_r;
	mj_r = tmp;
	tmp = invR / exp(wave_i * R);
	real_t coef_r = cos(wave_r * R) * tmp;
	real_t coef_i = sin(wave_r * R) * tmp;
	tmp = mj_r * coef_r - mj_i * coef_i;
	coef_i = mj_r * coef_i + mj_i * coef_r;
	coef_r = tmp;
	mj_r = (1 + wave_i * R) * invR * invR;
	mj_i = - wave_r * invR;
	pot_r += coef_r;
	pot_i += coef_i;
	tmp = mj_r * coef_r - mj_i * coef_i;
	coef_i = mj_r * coef_i + mj_i * coef_r;
	coef_r = tmp;
	ax_r += coef_r * dx;
	ax_i += coef_i * dx;
	ay_r += coef_r * dy;
	ay_i += coef_i * dy;
	az_r += coef_r * dz;
	az_i += coef_i * dz;
      }
    }
    Bi[i].TRG[0] += complex_t(pot_r, pot_i);
    Bi[i].TRG[1] -= complex_t(ax_r, ax_i);
    Bi[i].TRG[2] -= complex_t(ay_r, ay_i);
    Bi[i].TRG[3] -= complex_t(az_r, az_i);
#endif    
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
	real_t expikr = exp(wave_i * R) * R;
	real_t expikr_r = cos(wave_r * R) / expikr;
	real_t expikr_i = sin(wave_r * R) / expikr;
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
      }
    }
    Bi[i].TRG[0] += complex_t(pot_r, pot_i);
    Bi[i].TRG[1] += complex_t(ax_r, ax_i);
    Bi[i].TRG[2] += complex_t(ay_r, ay_i);
    Bi[i].TRG[3] += complex_t(az_r, az_i);
  }
}

void kernel::P2P(C_iter C) {
  for (B_iter Bi=C->BODY; Bi!=C->BODY+C->NBODY; Bi++) {
    complex_t pot = 0.0;
    complex_t ax = 0.0;
    complex_t ay = 0.0;
    complex_t az = 0.0;
    for (B_iter Bj=C->BODY; Bj!=C->BODY+C->NBODY; Bj++) {
      vec3 dX = Bi->X - Bj->X - Xperiodic;
      real_t R2 = norm(dX) + eps2;
      if (R2 != 0) {
	real_t R = sqrt(R2);
	complex_t coef1 = Bi->SRC * Bj->SRC * exp(I * wavek * R) / R;
	complex_t coef2 = (real_t(1.0) - I * wavek * R) * coef1 / R2;
	pot += coef1;
	ax += coef2 * dX[0];
	ay += coef2 * dX[1];
	az += coef2 * dX[2];
      }
    }
    Bi->TRG[0] += pot;
    Bi->TRG[1] += ax;
    Bi->TRG[2] += ay;
    Bi->TRG[3] += az;
  }
}
