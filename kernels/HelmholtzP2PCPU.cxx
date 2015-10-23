#include "kernel.h"
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
#if EXAFMM_USE_SIMD
  for ( ; i<=ni-NSIMD; i++) {
#if 0
    simdvec zero = 0.0;
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
    simdvec R2 = eps2;

    simdvec xj = Xperiodic[0];
    xi -= xj;
    simdvec yj = Xperiodic[1];
    yi -= yj;
    simdvec zj = Xperiodic[2];
    zi -= zj;

    simdvec x2 = Bj[0].X[0];
    x2 -= xi;
    simdvec y2 = Bj[0].X[1];
    y2 -= yi;
    simdvec z2 = Bj[0].X[2];
    z2 -= zi;
    simdvec mj_r = std::real(Bj[0].SRC);
    simdvec mj_i = std::imag(Bj[0].SRC);

    xj = x2;
    R2 += x2 * x2;
    yj = y2;
    R2 += y2 * y2;
    zj = z2;
    R2 += z2 * z2;
    simdvec invR;

    x2 = Bj[1].X[0];
    y2 = Bj[1].X[1];
    z2 = Bj[1].X[2];

    for (int j=0; j<nj-2; j++) {
      invR = rsqrt(R2);
      invR &= R2 > zero;
      R2 = eps2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;
    }
#endif
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
      real_t dx = Bj[j].X[0];
      real_t dy = Bj[j].X[1];
      real_t dz = Bj[j].X[2];
      real_t mj_r = std::real(Bj[j].SRC);
      real_t mj_i = std::imag(Bj[j].SRC);
      real_t R2 = eps2;
      dx -= xi;
      dy -= yi;
      dz -= zi;
      real_t xj = dx;
      R2 += dx * dx;
      real_t yj = dy;
      R2 += dy * dy;
      real_t zj = dz;
      R2 += dz * dz;
      if (R2 != 0) {
	real_t invR = 1/sqrt(R2);
	real_t tmp = mi_r * mj_r - mi_i * mj_i;
	mj_i = mi_r * mj_i + mi_i * mj_r;
	mj_r = tmp;
	tmp = invR / exp(wave_i / invR);
	real_t expikr_r = cos(wave_r / invR) * tmp;
	real_t expikr_i = sin(wave_r / invR) * tmp;
	real_t coef1_r = mj_r * expikr_r - mj_i * expikr_i;
	real_t coef1_i = mj_r * expikr_i + mj_i * expikr_r;
	real_t kr_r = (1 + wave_i / invR) * invR * invR;
	real_t kr_i = - wave_r * invR;
	real_t coef2_r = kr_r * coef1_r - kr_i * coef1_i;
	real_t coef2_i = kr_r * coef1_i + kr_i * coef1_r;
	pot_r += coef1_r;
	pot_i += coef1_i;
	ax_r -= coef2_r * dx;
	ax_i -= coef2_i * dx;
	ay_r -= coef2_r * dy;
	ay_i -= coef2_i * dy;
	az_r -= coef2_r * dz;
	az_i -= coef2_i * dz;
      }
    }    
    Bi[i].TRG[0] += complex_t(pot_r, pot_i);
    Bi[i].TRG[1] += complex_t(ax_r, ax_i);
    Bi[i].TRG[2] += complex_t(ay_r, ay_i);
    Bi[i].TRG[3] += complex_t(az_r, az_i);
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
