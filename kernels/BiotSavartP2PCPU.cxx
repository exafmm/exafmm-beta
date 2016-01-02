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
      for ( ; i<ni; i++) {
	kreal_t ax = 0;
	kreal_t ay = 0;
	kreal_t az = 0;
	for (int j=0; j<nj; j++) {
	  vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t invR2 = 1.0 / R2;
	    real_t S2 = 2 * Bj[j].SRC[3] * Bj[j].SRC[3];
	    real_t RS = R2 / S2;
	    real_t cutoff = invR2 * std::sqrt(invR2) * (erf( std::sqrt(RS) ) - std::sqrt(4 / M_PI * RS) * std::exp(-RS));
	    ax += (dX[1] * Bj[j].SRC[2] - dX[2] * Bj[j].SRC[1]) * cutoff;
	    ay += (dX[2] * Bj[j].SRC[0] - dX[0] * Bj[j].SRC[2]) * cutoff;
	    az += (dX[0] * Bj[j].SRC[1] - dX[1] * Bj[j].SRC[0]) * cutoff;
	  }
	}
	Bi[i].TRG[0] = 1;
	Bi[i].TRG[1] += ax;
	Bi[i].TRG[2] += ay;
	Bi[i].TRG[3] += az; 
      }
    }

    void P2P(C_iter C) {
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
      for ( ; i<n; i++) {
	kreal_t ax = 0;
	kreal_t ay = 0;
	kreal_t az = 0;
	for (int j=i+1; j<n; j++) {
	  vec3 dX = B[j].X - B[i].X;
	  real_t R2 = norm(dX) + eps2;
	  if (R2 != 0) {
	    real_t invR2 = 1.0 / R2;
	    real_t S2 = 2 * B[j].SRC[3] * B[j].SRC[3];
	    real_t RS = R2 / S2;
	    real_t cutoff = invR2 * std::sqrt(invR2) * (erf( std::sqrt(RS) ) - std::sqrt(4 / M_PI * RS) * std::exp(-RS));
	    ax += (dX[1] * B[j].SRC[2] - dX[2] * B[j].SRC[1]) * cutoff;
	    ay += (dX[2] * B[j].SRC[0] - dX[0] * B[j].SRC[2]) * cutoff;
	    az += (dX[0] * B[j].SRC[1] - dX[1] * B[j].SRC[0]) * cutoff;
	  }
	}
	B[i].TRG[0] = 1;
	B[i].TRG[1] += ax;
	B[i].TRG[2] += ay;
	B[i].TRG[3] += az;
      }
    }
  }
}
