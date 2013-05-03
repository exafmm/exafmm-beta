#include "kernel.h"
#include "matern.h"

const real_t DENOM = 1 / (std::pow(2,NU-1) * tgamma(NU));     //!< Denominator in Matern kernel

void Kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) const {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  for (int i=0; i<ni; i++) {
    kreal_t pot = 0; 
    for (int j=0; j<nj; j++) {
      vec3 dX = (Bi[i].X - Bj[j].X - Xperiodic) / SIGMA;
      real_t R = std::sqrt(2 * NU * norm(dX));
      if (R != 0) {
        real_t phi = Bi[i].SRC * Bj[j].SRC * std::pow(R,NU) * cyl_bessel_k(NU,R) * DENOM;
        pot += phi;
        if (mutual) {
          Bj[j].TRG[0] += phi;
        }
      } else {
        pot += Bi[i].SRC * Bj[j].SRC;
      }
    }
    Bi[i].TRG[0] += pot;
  }
}

void Kernel::P2P(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  for (int i=0; i<n; i++) {
    kreal_t pot = 0;
    for (int j=i; j<n; j++) {
      vec3 dX = (B[i].X - B[j].X) / SIGMA;
      real_t R = std::sqrt(2 * NU * norm(dX));
      if (R != 0) {
        real_t phi = B[i].SRC * B[j].SRC * std::pow(R,NU) * cyl_bessel_k(NU,R) * DENOM;
        pot += phi;
        B[j].TRG[0] += phi;
      } else {
        pot += B[i].SRC * B[j].SRC;
      }
    }
    B[i].TRG[0] += pot;
  }
}
