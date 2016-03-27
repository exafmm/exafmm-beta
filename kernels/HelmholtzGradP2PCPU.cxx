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
      complex_t one(1.0,0.0);
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
      for ( ; i<n; i++) {
        complex_t pot(0.0, 0.0);
        complex_t ax(0.0, 0.0);
        complex_t ay(0.0, 0.0);
        complex_t az(0.0, 0.0);
        complex_t mi = B[i].SRC;
        for (int j=i+1; j<n; j++) {
          complex_t mj = B[j].SRC;
          vec3 dX = B[j].X - B[i].X;
          real_t R2 = norm(dX) + eps2;
          real_t R = sqrt(R2);
          complex_t src2 = mi * mj;
          complex_t expikr = std::exp(I * wavek * R) / R;
          complex_t coef1 = src2 * expikr;
          complex_t coef2 = (one - I * wavek * R) / R2 * coef1;
          pot += coef1;
          ax += coef2 * dX[0];
          ay += coef2 * dX[1];
          az += coef2 * dX[2];
          B[j].TRG[0] += coef1;
          B[j].TRG[1] += coef2 * dX[0];
          B[j].TRG[2] += coef2 * dX[1];
          B[j].TRG[3] += coef2 * dX[2];
        }
        B[i].TRG[0] += pot;
        B[i].TRG[1] -= ax;
        B[i].TRG[2] -= ay;
        B[i].TRG[3] -= az;
      }
    }
  }
}
