#include "kernel.h"

real_t kernel::eps2;
complex_t kernel::wavek;
vec3 kernel::Xperiodic;

const complex_t I(0.,1.);

void kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) {
  assert(mutual == 0);
  for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {
    complex_t pot = 0.0;
    complex_t ax = 0.0;
    complex_t ay = 0.0;
    complex_t az = 0.0;
    for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {
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
