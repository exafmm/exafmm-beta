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
  for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {
    real_t pot_r = 0.0;
    real_t pot_i = 0.0;
    real_t ax_r = 0.0;
    real_t ax_i = 0.0;
    real_t ay_r = 0.0;
    real_t ay_i = 0.0;
    real_t az_r = 0.0;
    real_t az_i = 0.0;
    real_t BiSRC_r = std::real(Bi->SRC);
    real_t BiSRC_i = std::imag(Bi->SRC);
    for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {
      real_t BjSRC_r = std::real(Bj->SRC);
      real_t BjSRC_i = std::imag(Bj->SRC);
      vec3 dX = Bi->X - Bj->X - Xperiodic;
      real_t R2 = norm(dX) + eps2;
      if (R2 != 0) {
	real_t R = sqrt(R2);
	real_t src2_r = BiSRC_r * BjSRC_r - BiSRC_i * BjSRC_i;
	real_t src2_i = BiSRC_r * BjSRC_i + BiSRC_i * BjSRC_r;
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
    Bi->TRG[0] += complex_t(pot_r, pot_i);
    Bi->TRG[1] += complex_t(ax_r, ax_i);
    Bi->TRG[2] += complex_t(ay_r, ay_i);
    Bi->TRG[3] += complex_t(az_r, az_i);
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
