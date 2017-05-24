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
#if EXAFMM_USE_SIMD
      simdvec wave_rvec = wave_r;
      simdvec wave_ivec = wave_i;
      for ( ; i<=ni-NSIMD; i+=NSIMD) {
        simdvec zero = 0.0;
        simdvec one = 1.0;
        simdvec tnear = nearpd;
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
        simdvec patchi = SIMD<simdvec,6,NSIMD>::setBody(Bi,i);


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
          simdvec patchj = Bj[j].PATCH;
          simdvec R2 = eps2;
          R2 += dx * dx;
          simdvec mj_r = std::real(Bj[j].SRC * Bj[j].QWEIGHT);
          R2 += dy * dy;
          simdvec mj_i = std::imag(Bj[j].SRC * Bj[j].QWEIGHT);
          R2 += dz * dz;
          simdvec invR = rsqrt(R2);
          simdvec R = one / invR;
          invR &= R2 > zero;
          R &= R2 > zero;
          invR &= R > tnear;
          R &= R > tnear;
          simdvec tmp = invR;
          invR &= (patchi > patchj);
          tmp  &= (patchi < patchj);
          invR += tmp;
          tmp = R;
          R &= (patchi > patchj);
          tmp &= (patchi < patchj);
          R += tmp;
          tmp = mi_r * mj_r - mi_i * mj_i;
          mj_i = mi_r * mj_i + mi_i * mj_r;
          mj_r = tmp;
          tmp = invR / exp(wave_ivec * R);
          simdvec coef_r = cos(wave_rvec * R) * tmp;
          simdvec coef_i = sin(wave_rvec * R) * tmp;
          tmp = mj_r * coef_r - mj_i * coef_i;
          coef_i = mj_r * coef_i + mj_i * coef_r;
          coef_r = tmp;
          pot_r += coef_r;
          pot_i += coef_i;
          if (mutual) Bj[j].TRG[0] += kcomplex_t(sum(coef_r), sum(coef_i));
        }
        for (int k=0; k<NSIMD; k++) {
          Bi[i+k].TRG[0] += transpose(pot_r, pot_i, k);
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
        real_t mi_r = std::real(Bi[i].SRC * Bi[i].QWEIGHT);
        real_t mi_i = std::imag(Bi[i].SRC * Bi[i].QWEIGHT);
        for (int j=0; j<nj; j++) {  
          vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
          real_t R2 = norm(dX) + eps2;
          if(Bi[i].PATCH != Bj[j].PATCH) {
            real_t mj_r = std::real(Bj[j].SRC * Bj[j].QWEIGHT);
            real_t mj_i = std::imag(Bj[j].SRC * Bj[j].QWEIGHT);
            if (R2 != 0) {
              real_t R = sqrt(R2);
              if(R <= nearpd) {
                real_t coef1_r = 0;
                real_t coef1_i = 0;
                for (int ll = 0; ll < nhdgqp; ++ll) {
                  vec3 dX_near = Bi[i].X - Bj[j].GAUSS_NEAR[ll] - Xperiodic;
                  complex_t mj_near = (real_t)0.5 * ws[ll] * ipolator_near[Bj[j].POINT_LOC][ll]/((real_t)4.0*(real_t)M_PI);
                  mj_r = std::real(mj_near);
                  mj_i = std::imag(mj_near);
                  real_t RR = sqrt(norm(dX_near));
                  real_t src2_r = mi_r * mj_r - mi_i * mj_i;
                  real_t src2_i = mi_r * mj_i + mi_i * mj_r;
                  real_t expikr_r_= std::exp(wave_i * RR) * RR;     
                  real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
                  real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
                  coef1_r += src2_r * expikr_r - src2_i * expikr_i;
                  coef1_i += src2_r * expikr_i + src2_i * expikr_r;                
                }              
                mj_r = std::real(Bj[j].SRC);
                mj_i = std::imag(Bj[j].SRC);
                pot_r += mj_r * coef1_r - mj_i * coef1_i;
                pot_i += mj_r * coef1_i + mj_i * coef1_r;
              } else {
                real_t src2_r = mi_r * mj_r - mi_i * mj_i;
                real_t src2_i = mi_r * mj_i + mi_i * mj_r;      
                real_t expikr_r_= std::exp(wave_i * R) * R;     
                real_t expikr_r = std::cos(wave_r * R) / expikr_r_;
                real_t expikr_i = std::sin(wave_r * R) / expikr_r_;   
                real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
                real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
                pot_r += coef1_r;
                pot_i += coef1_i;
                if (mutual) {
                  Bj[j].TRG[0] += complex_t(coef1_r, coef1_i);
                }
              }
            }
          }
        }
        Bi[i].TRG[0] += complex_t(pot_r, pot_i);
      }
#if EXAFMM_USE_SIMD 
      for (i=0 ; i<=ni-NSIMD; i++) { 
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
          real_t R = sqrt(R2);
          if(R2 != 0 && R < nearpd && Bi[i].PATCH != Bj[j].PATCH) {
            real_t coef1_r = 0;
            real_t coef1_i = 0;
            for (int ll = 0; ll < nhdgqp; ++ll) {
              vec3 dX_near = Bi[i].X - Bj[j].GAUSS_NEAR[ll] - Xperiodic;
              complex_t mj_near = (real_t)0.5 * Bj[j].SRC * ws[ll] * ipolator_near[Bj[j].POINT_LOC][ll]/((real_t)4.0*(real_t)M_PI);
              real_t mj_r = std::real(mj_near);
              real_t mj_i = std::imag(mj_near);
              real_t RR = sqrt(norm(dX_near));
              real_t src2_r = mi_r * mj_r - mi_i * mj_i;
              real_t src2_i = mi_r * mj_i + mi_i * mj_r;
              real_t expikr_r_= std::exp(wave_i * RR) * RR;     
              real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
              real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
              coef1_r = src2_r * expikr_r - src2_i * expikr_i;
              coef1_i = src2_r * expikr_i + src2_i * expikr_r;
              pot_r += coef1_r;
              pot_i += coef1_i;
            }
            if (mutual) {
              Bj[j].TRG[0] += complex_t(coef1_r, coef1_i);
            }
        }        
      }
      Bi[i].TRG[0] += complex_t(pot_r, pot_i);
    }
#endif
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
        simdvec tnear = nearpd;
        ksimdvec pot_r = zero;
        ksimdvec pot_i = zero;
        ksimdvec ax_r = zero;
        ksimdvec ax_i = zero;
        ksimdvec ay_r = zero;
        ksimdvec ay_i = zero;
        ksimdvec az_r = zero;
        ksimdvec az_i = zero;

        simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(B,i);
        simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(B,i);
        simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(B,i);
        simdvec mi_r = SIMD<simdvec,4,NSIMD>::setBody(B,i);
        simdvec mi_i = SIMD<simdvec,5,NSIMD>::setBody(B,i);
        simdvec patchi = SIMD<simdvec,6,NSIMD>::setBody(B,i);


        simdvec dx = Xperiodic[0];
        xi -= dx;
        simdvec dy = Xperiodic[1];
        yi -= dy;
        simdvec dz = Xperiodic[2];
        zi -= dz;

        for (int j=0; j<n; j++) {
          dx = B[j].X[0];
          dx -= xi;
          dy = B[j].X[1];
          dy -= yi;
          dz = B[j].X[2];
          dz -= zi;
          simdvec patchj = B[j].PATCH;
          simdvec R2 = eps2;
          R2 += dx * dx;
          simdvec mj_r = std::real(B[j].SRC * B[j].QWEIGHT);
          R2 += dy * dy;
          simdvec mj_i = std::imag(B[j].SRC * B[j].QWEIGHT);
          R2 += dz * dz;
          simdvec invR = rsqrt(R2);
          simdvec R = one / invR;
          invR &= R2 > zero;
          R &= R2 > zero;
          invR &= R > tnear;
          R &= R > tnear;
          simdvec tmp = invR;
          invR &= (patchi > patchj);
          tmp  &= (patchi < patchj);
          invR += tmp;
          tmp = R;
          R &= (patchi > patchj);
          tmp &= (patchi < patchj);
          R += tmp;
          tmp = mi_r * mj_r - mi_i * mj_i;
          mj_i = mi_r * mj_i + mi_i * mj_r;
          mj_r = tmp;
          tmp = invR / exp(wave_ivec * R);
          simdvec coef_r = cos(wave_rvec * R) * tmp;
          simdvec coef_i = sin(wave_rvec * R) * tmp;
          tmp = mj_r * coef_r - mj_i * coef_i;
          coef_i = mj_r * coef_i + mj_i * coef_r;
          coef_r = tmp;
          pot_r += coef_r;
          pot_i += coef_i;
          B[j].TRG[0] += kcomplex_t(sum(coef_r), sum(coef_i));
        }
        for (int k=0; k<NSIMD; k++) {
          B[i+k].TRG[0] += transpose(pot_r, pot_i, k);
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
        real_t mi_r = std::real(B[i].SRC * B[i].QWEIGHT);
        real_t mi_i = std::imag(B[i].SRC * B[i].QWEIGHT);
        for (int j=i+1; j<n; j++) {
          if(B[i].PATCH != B[j].PATCH) { 
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
                complex_t mj_near = (real_t)0.5 * ws[ll] * ipolator_near[B[j].POINT_LOC][ll]/((real_t)4.0*(real_t)M_PI);
                mj_r = std::real(mj_near);
                mj_i = std::imag(mj_near);
                real_t RR = sqrt(norm(dX_near));
                real_t src2_r = mi_r * mj_r - mi_i * mj_i;
                real_t src2_i = mi_r * mj_i + mi_i * mj_r;
                real_t expikr_r_= std::exp(wave_i * RR) * RR;     
                real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
                real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
                coef1_r += src2_r * expikr_r - src2_i * expikr_i;
                coef1_i += src2_r * expikr_i + src2_i * expikr_r;
                // pot_r += coef1_r;
                // pot_i += coef1_i;
              }
                mj_r = std::real(B[j].SRC);
                mj_i = std::imag(B[j].SRC);
                pot_r += mj_r * coef1_r - mj_i * coef1_i;
                pot_i += mj_r * coef1_i + mj_i * coef1_r;
            } else {
              real_t src2_r = mi_r * mj_r - mi_i * mj_i;
              real_t src2_i = mi_r * mj_i + mi_i * mj_r;
              real_t expikr_r_= std::exp(wave_i * R) * R;     
              real_t expikr_r = std::cos(wave_r * R) / expikr_r_;
              real_t expikr_i = std::sin(wave_r * R) / expikr_r_;   
              real_t coef1_r = src2_r * expikr_r - src2_i * expikr_i;
              real_t coef1_i = src2_r * expikr_i + src2_i * expikr_r;
              pot_r += coef1_r;
              pot_i += coef1_i;
              B[j].TRG[0] += complex_t(coef1_r, coef1_i);
            }
          }
        }
        B[i].TRG[0] += complex_t(pot_r, pot_i);
      }
#if EXAFMM_USE_SIMD 
      for (i=0 ; i<=n-NSIMD; i++) { 
        real_t pot_r = 0.0;
        real_t pot_i = 0.0;
        real_t ax_r = 0.0;
        real_t ax_i = 0.0;
        real_t ay_r = 0.0;
        real_t ay_i = 0.0;
        real_t az_r = 0.0;
        real_t az_i = 0.0;
        real_t mi_r = std::real(B[i].SRC * B[i].QWEIGHT);
        real_t mi_i = std::imag(B[i].SRC * B[i].QWEIGHT);
        for (int j=0; j<n; j++) {  
          vec3 dX = B[i].X - B[j].X - Xperiodic;
          real_t R2 = norm(dX) + eps2;
          real_t R = sqrt(R2);
          if(R2 != 0 && R < nearpd && B[i].PATCH != B[j].PATCH) {
            real_t coef1_r = 0;
            real_t coef1_i = 0;
            for (int ll = 0; ll < nhdgqp; ++ll) {
              vec3 dX_near = B[i].X - B[j].GAUSS_NEAR[ll] - Xperiodic;
              complex_t mj_near = (real_t)0.5 * B[j].SRC * ws[ll] * ipolator_near[B[j].POINT_LOC][ll]/((real_t)4.0*(real_t)M_PI);
              real_t mj_r = std::real(mj_near);
              real_t mj_i = std::imag(mj_near);
              real_t RR = sqrt(norm(dX_near));
              real_t src2_r = mi_r * mj_r - mi_i * mj_i;
              real_t src2_i = mi_r * mj_i + mi_i * mj_r;
              real_t expikr_r_= std::exp(wave_i * RR) * RR;     
              real_t expikr_r = std::cos(wave_r * RR) / expikr_r_;
              real_t expikr_i = std::sin(wave_r * RR) / expikr_r_;   
              coef1_r = src2_r * expikr_r - src2_i * expikr_i;
              coef1_i = src2_r * expikr_i + src2_i * expikr_r;
              pot_r += coef1_r;
              pot_i += coef1_i;
            }
            B[j].TRG[0] += complex_t(coef1_r, coef1_i);
        }        
      }
      B[i].TRG[0] += complex_t(pot_r, pot_i);
    }
#endif
    }
  }
}
