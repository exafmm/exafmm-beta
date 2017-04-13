#include "kernel.h"
#include "simdvec.h"

namespace exafmm {
  namespace kernel {
    real_t eps2;
    complex_t wavek;
    vec3 Xperiodic;
    int nhdgqp;                                          
    int nipp;                                            
    std::vector<std::vector<real_t> > ipolator_near;     
    std::vector<real_t> ws;       
    real_t nearpd;                       

    const complex_t I(0.,1.);
    void P2P(C_iter Ci, C_iter Cj, bool mutual) {
      real_t wave_r = std::real(wavek);
      real_t wave_i = std::imag(wavek);
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;
      for ( ; i<ni; i++) {
        real_t pot_r = 0.0;
        real_t pot_i = 0.0;
        real_t ax_r = 0.0;
        real_t ax_i = 0.0;
        real_t ay_r = 0.0;
        real_t ay_i = 0.0;
        real_t az_r = 0.0;
        real_t az_i = 0.0;
        real_t mi_r = std::real(Bi[i].SRC * Bi[i].QWEIGHT);
        real_t mi_i = std::imag(Bi[i].SRC * Bi[i].QWEIGHT);
        for (int j=0; j<nj; j++) {  
          vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
          real_t R2 = norm(dX) + eps2;
          if(Bi[i].PATCH == Bj[j].PATCH) continue;
          real_t mj_r = std::real(Bj[j].SRC * Bj[j].QWEIGHT);
          real_t mj_i = std::imag(Bj[j].SRC * Bj[j].QWEIGHT);
          if (R2 != 0) {
            real_t R = sqrt(R2);
            if(R < nearpd) {
              real_t coef1_r = 0;
              real_t coef1_i = 0;
              for (int ll = 0; ll < nhdgqp; ++ll) {
                vec3 dX_near = Bi[i].X - Bj[j].GAUSS_NEAR[ll] - Xperiodic;
                complex_t mj_near = (real_t)0.5 * Bj[j].SRC * ws[ll] * ipolator_near[Bj[j].POINT_LOC][ll]/(real_t)4.0/(real_t)M_PI;
                mj_r = std::real(mj_near);
                mj_i = std::imag(mj_near);
                real_t RR = sqrt(norm(dX_near));
                real_t src2_r = mi_r * mj_r - mi_i * mj_i;
                real_t src2_i = mi_r * mj_i + mi_i * mj_r;
                real_t expikr_r_= std::exp(wave_r * RR) / RR;     
                real_t expikr_r = std::cos(wave_i * RR) * expikr_r_;
                real_t expikr_i = std::sin(wave_i * RR) * expikr_r_;   
                coef1_r = src2_r * expikr_r - src2_i * expikr_i;
                coef1_i =coef1_i = src2_r * expikr_i + src2_i * expikr_r;
                pot_r += coef1_r;
                pot_i += coef1_i;
              }
              continue;
            }
            real_t src2_r = mi_r * mj_r - mi_i * mj_i;
            real_t src2_i = mi_r * mj_i + mi_i * mj_r;      
            real_t expikr_r_= std::exp(wave_r * R) / R;     
            real_t expikr_r = std::cos(wave_i * R) * expikr_r_;
            real_t expikr_i = std::sin(wave_i * R) * expikr_r_;   
            real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
            real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
            pot_r += coef1_r;
            pot_i += coef1_i;
            if (mutual) {
              Bj[j].TRG[0] += complex_t(coef1_r, coef1_i);
            }
          }
        }
        Bi[i].TRG[0] += complex_t(pot_r, pot_i);
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
        real_t mi_r = std::real(B[i].SRC * B[i].QWEIGHT);
        real_t mi_i = std::imag(B[i].SRC * B[i].QWEIGHT);
        for (int j=i+1; j<n; j++) {
          if(B[i].PATCH == B[j].PATCH) continue;  
          real_t mj_r = std::real(B[j].SRC * B[j].QWEIGHT);
          real_t mj_i = std::imag(B[j].SRC * B[j].QWEIGHT);
          vec3 dX = B[j].X - B[i].X;
          real_t R2 = norm(dX) + eps2;
          real_t R = sqrt(R2);
          if(R < nearpd) {
            real_t coef1_r = 0;
            real_t coef1_i = 0;
            for (int ll = 0; ll < nhdgqp; ++ll) {
              vec3 dX_near = B[i].X - B[j].GAUSS_NEAR[ll] - Xperiodic;
              complex_t mj_near = (real_t)0.5 * B[j].SRC * ws[ll] * ipolator_near[B[j].POINT_LOC][ll]/(real_t)4.0/(real_t)M_PI;
              mj_r = std::real(mj_near);
              mj_i = std::imag(mj_near);
              real_t RR = sqrt(norm(dX_near));
              real_t src2_r = mi_r * mj_r - mi_i * mj_i;
              real_t src2_i = mi_r * mj_i + mi_i * mj_r;
              real_t expikr_r_= std::exp(wave_r * RR) / RR;     
              real_t expikr_r = std::cos(wave_i * RR) * expikr_r_;
              real_t expikr_i = std::sin(wave_i * RR) * expikr_r_;   
              coef1_r = src2_r * expikr_r - src2_i * expikr_i;
              coef1_i = src2_r * expikr_i + src2_i * expikr_r;
              pot_r += coef1_r;
              pot_i += coef1_i;
            }
            continue;
          }
          real_t src2_r = mi_r * mj_r - mi_i * mj_i;
          real_t src2_i = mi_r * mj_i + mi_i * mj_r;
          real_t expikr_r_= std::exp(wave_r * R) / R;     
          real_t expikr_r = std::cos(wave_i * R) * expikr_r_;
          real_t expikr_i = std::sin(wave_i * R) * expikr_r_;   
          real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
          real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
          pot_r += coef1_r;
          pot_i += coef1_i;
          B[j].TRG[0] += complex_t(coef1_r, coef1_i);
        }
        B[i].TRG[0] += complex_t(pot_r, pot_i);
      }
    }
  }
}
