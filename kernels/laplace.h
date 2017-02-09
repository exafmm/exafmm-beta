#ifndef laplace_h
#define laplace_h
#include "spherical.h"
#if EXAFMM_USE_SIMD
#include "simdvec.h"
#endif

namespace exafmm {
  class Kernel {
  public:
    const int P;
    const int NTERM;
    real_t    eps2;
    complex_t wavek;
    vec3      Xperiodic;

    Kernel(real_t _eps2, complex_t _wavek) : P(Pmax), NTERM(P*(P+1)/2), eps2(_eps2), wavek(_wavek) {
      Xperiodic = 0;
    }

    void init() {}
    void finalize() {}

    void normalize(Bodies & bodies) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
        B->TRG /= B->SRC;
      }
    }

    void P2P(C_iter Ci, C_iter Cj) {
      B_iter Bi = Ci->BODY;
      B_iter Bj = Cj->BODY;
      int ni = Ci->NBODY;
      int nj = Cj->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      for ( ; i<=ni-NSIMD; i+=NSIMD) {
        simdvec zero = 0.0;
        ksimdvec pot = zero;
        ksimdvec ax = zero;
        ksimdvec ay = zero;
        ksimdvec az = zero;

        simdvec xi = SIMD<simdvec,B_iter,0,NSIMD>::setBody(Bi,i);
        simdvec yi = SIMD<simdvec,B_iter,1,NSIMD>::setBody(Bi,i);
        simdvec zi = SIMD<simdvec,B_iter,2,NSIMD>::setBody(Bi,i);
        simdvec mi = SIMD<simdvec,B_iter,3,NSIMD>::setBody(Bi,i);

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
          simdvec mj = Bj[j].SRC;

          simdvec R2 = eps2;
          xj = dx;
          R2 += dx * dx;
          yj = dy;
          R2 += dy * dy;
          zj = dz;
          R2 += dz * dz;
          simdvec invR = rsqrt(R2);
          invR &= R2 > zero;

          mj *= invR * mi;
          pot += mj;
          invR = invR * invR * mj;

          xj *= invR;
          ax += xj;

          yj *= invR;
          ay += yj;

          zj *= invR;
          az += zj;
        }
        for (int k=0; k<NSIMD; k++) {
          Bi[i+k].TRG[0] += transpose(pot, k);
          Bi[i+k].TRG[1] += transpose(ax, k);
          Bi[i+k].TRG[2] += transpose(ay, k);
          Bi[i+k].TRG[3] += transpose(az, k);
        }
      }
#endif
      for ( ; i<ni; i++) {
        kreal_t pot = 0;
        kreal_t ax = 0;
        kreal_t ay = 0;
        kreal_t az = 0;
        for (int j=0; j<nj; j++) {
          vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
          real_t R2 = norm(dX) + eps2;
          if (R2 != 0) {
            real_t invR2 = 1.0 / R2;
            real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
            dX *= invR2 * invR;
            pot += invR;
            ax += dX[0];
            ay += dX[1];
            az += dX[2];
          }
        }
        Bi[i].TRG[0] += pot;
        Bi[i].TRG[1] -= ax;
        Bi[i].TRG[2] -= ay;
        Bi[i].TRG[3] -= az;
      }
    }

    void P2M(C_iter C) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
        vec3 dX = B->X - C->X;
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(P, rho, alpha, beta, Ynm, YnmTheta);
        for (int n=0; n<P; n++) {
          for (int m=0; m<=n; m++) {
            int nm  = n * n + n - m;
            int nms = n * (n + 1) / 2 + m;
            C->M[nms] += B->SRC * Ynm[nm];
          }
        }
      }
    }

    void M2M(C_iter Ci, C_iter C0) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (C_iter Cj=C0+Ci->ICHILD; Cj!=C0+Ci->ICHILD+Ci->NCHILD; Cj++) {
        vec3 dX = Ci->X - Cj->X;
        real_t rho, alpha, beta;
        cart2sph(dX, rho, alpha, beta);
        evalMultipole(P, rho, alpha, beta, Ynm, YnmTheta);
        for (int j=0; j<P; j++) {
          for (int k=0; k<=j; k++) {
            int jks = j * (j + 1) / 2 + k;
            complex_t M = 0;
            for (int n=0; n<=j; n++) {
              for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
                int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                int nm    = n * n + n - m;
                M += Cj->M[jnkms] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
              }
              for (int m=k; m<=std::min(n,j+k-n); m++) {
                int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                int nm    = n * n + n - m;
                M += std::conj(Cj->M[jnkms]) * Ynm[nm] * real_t(oddOrEven(k+n+m));
              }
            }
            Ci->M[jks] += M;
          }
        }
      }
    }

    void M2L(C_iter Ci, C_iter Cj) {
      complex_t Ynmi[P*P], Ynmj[P*P];
      vec3 dX = Ci->X - Cj->X - Xperiodic;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalLocal(P, rho, alpha, beta, Ynmi);
      for (int j=0; j<P; j++) {
        real_t Cnm = oddOrEven(j);
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          complex_t Li = 0, Lj = 0;
          for (int n=0; n<P-j; n++) {
            for (int m=-n; m<0; m++) {
              int nms  = n * (n + 1) / 2 - m;
              int jnkm = (j + n) * (j + n) + j + n + m - k;
              Li += std::conj(Cj->M[nms]) * Cnm * Ynmi[jnkm];
            }
            for (int m=0; m<=n; m++) {
              int nms  = n * (n + 1) / 2 + m;
              int jnkm = (j + n) * (j + n) + j + n + m - k;
              real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
              Li += Cj->M[nms] * Cnm2 * Ynmi[jnkm];
            }
          }
          Ci->L[jks] += Li;
        }
      }
    }

    void L2L(C_iter Ci, C_iter C0) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      C_iter Cj = C0 + Ci->IPARENT;
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(P, rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          complex_t L = 0;
          for (int n=j; n<P; n++) {
            for (int m=j+k-n; m<0; m++) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nms  = n * (n + 1) / 2 - m;
              L += std::conj(Cj->L[nms]) * Ynm[jnkm] * real_t(oddOrEven(k));
            }
            for (int m=0; m<=n; m++) {
              if( n-j >= abs(m-k) ) {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nms  = n * (n + 1) / 2 + m;
                L += Cj->L[nms] * Ynm[jnkm] * real_t(oddOrEven((m-k)*(m<k)));
              }
            }
          }
          Ci->L[jks] += L;
        }
      }
    }

    void L2P(C_iter Ci) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NBODY; B++) {
        vec3 dX = B->X - Ci->X + EPS;
        vec3 spherical = 0;
        vec3 cartesian = 0;
        real_t r, theta, phi;
        cart2sph(dX, r, theta, phi);
        evalMultipole(P, r, theta, phi, Ynm, YnmTheta);
        B->TRG /= B->SRC;
        for (int n=0; n<P; n++) {
          int nm  = n * n + n;
          int nms = n * (n + 1) / 2;
          B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
          spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
          spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
          for( int m=1; m<=n; m++) {
            nm  = n * n + n + m;
            nms = n * (n + 1) / 2 + m;
            B->TRG[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]);
            spherical[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) / r * n;
            spherical[1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
            spherical[2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;
          }
        }
        sph2cart(r, theta, phi, spherical, cartesian);
        B->TRG[1] += cartesian[0];
        B->TRG[2] += cartesian[1];
        B->TRG[3] += cartesian[2];
      }
    }
  };
}
#endif
