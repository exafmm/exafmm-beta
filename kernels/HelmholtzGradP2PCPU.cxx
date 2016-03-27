#include "kernel.h"
#include "simdvec.h"

namespace exafmm {
  namespace kernel {
    real_t eps2;
    complex_t wavek;
    vec3 Xperiodic;

    const complex_t I(0.,1.);

    void P2P(C_iter Ci, C_iter Cj, bool mutual) {
      complex_t one(1.0,0.0);
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;
      for ( ; i<ni; i++) {
	complex_t pot(0.0, 0.0);
	complex_t ax(0.0, 0.0);
	complex_t ay(0.0, 0.0);
	complex_t az(0.0, 0.0);
	complex_t mi = Bi[i].SRC;
	for (int j=0; j<nj; j++) {
	  complex_t mj = Bj[j].SRC;
	  vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t R = sqrt(R2);
	    complex_t src2 = mi * mj;
	    complex_t expikr = std::exp(I * wavek * R) / R;
	    complex_t coef1 = src2 * expikr;
	    complex_t coef2 = (one - I * wavek * R) / R2 * coef1;
	    pot += coef1;
	    ax += coef2 * dX[0];
	    ay += coef2 * dX[1];
	    az += coef2 * dX[2];
	  }
	}
	Bi[i].TRG[0] += pot;
	Bi[i].TRG[1] += ax;
	Bi[i].TRG[2] += ay;
	Bi[i].TRG[3] += az;
      }
    }

    void P2P(C_iter C) {
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
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
