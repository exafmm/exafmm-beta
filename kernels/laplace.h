#ifndef laplace_h
#define laplace_h
#include "spherical.h"
#if EXAFMM_USE_SIMD
#include "simdvec.h"
#endif

namespace exafmm {
  template<int _P>
  class LaplaceKernel : public KernelBase {
  public:
    static const Equation equation = Laplace;                   //!< Set equation to Laplace
    static const int P = _P;                                    //!< Set order of expansion
    static const int NTERM = P*(P+1)/2;                         //!< # of terms in Laplace expansion
    typedef vec<NTERM,complex_t> vecP;                          //!< Vector type for expansion terms
    typedef std::vector<Body<equation> > Bodies;                //!< Vector of bodies for Laplace
    typedef typename Bodies::iterator B_iter;                   //!< Iterator of body vector
    typedef std::vector<Cell<B_iter,vecP> > Cells;              //!< Vector of cells for Laplace
    typedef typename Cells::iterator C_iter;                    //!< Iterator of cell vector
    using KernelBase::Xperiodic;

    static void init() {}
    static void finalize() {}

    static void normalize(Bodies & bodies) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
        B->TRG /= B->SRC;
      }
    }

    static void P2P(C_iter Ci, C_iter Cj, bool mutual) {
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
          if (mutual) Bj[j].TRG[0] += sum(mj);
          invR = invR * invR * mj;

          xj *= invR;
          ax += xj;
          if (mutual) Bj[j].TRG[1] -= sum(xj);

          yj *= invR;
          ay += yj;
          if (mutual) Bj[j].TRG[2] -= sum(yj);

          zj *= invR;
          az += zj;
          if (mutual) Bj[j].TRG[3] -= sum(zj);
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
            if (mutual) {
              Bj[j].TRG[0] += invR;
              Bj[j].TRG[1] += dX[0];
              Bj[j].TRG[2] += dX[1];
              Bj[j].TRG[3] += dX[2];
            }
          }
        }
        Bi[i].TRG[0] += pot;
        Bi[i].TRG[1] -= ax;
        Bi[i].TRG[2] -= ay;
        Bi[i].TRG[3] -= az;
      }
    }

    static void P2P(C_iter C) {
      B_iter B = C->BODY;
      int n = C->NBODY;
      int i = 0;
#if EXAFMM_USE_SIMD
      for ( ; i<=n-NSIMD; i+=NSIMD) {
        simdvec zero = 0;
        simdvec one = 1;
        ksimdvec pot = zero;
        ksimdvec ax = zero;
        ksimdvec ay = zero;
        ksimdvec az = zero;

        simdvec index = SIMD<simdvec,B_iter,0,NSIMD>::setIndex(i);
        simdvec xi = SIMD<simdvec,B_iter,0,NSIMD>::setBody(B,i);
        simdvec yi = SIMD<simdvec,B_iter,1,NSIMD>::setBody(B,i);
        simdvec zi = SIMD<simdvec,B_iter,2,NSIMD>::setBody(B,i);
        simdvec mi = SIMD<simdvec,B_iter,3,NSIMD>::setBody(B,i);
        for (int j=i+1; j<n; j++) {
          simdvec dx = B[j].X[0];
          dx -= xi;
          simdvec dy = B[j].X[1];
          dy -= yi;
          simdvec dz = B[j].X[2];
          dz -= zi;
          simdvec mj = B[j].SRC;

          simdvec R2 = eps2;
          simdvec xj = dx;
          R2 += dx * dx;
          simdvec yj = dy;
          R2 += dy * dy;
          simdvec zj = dz;
          R2 += dz * dz;
          simdvec invR = one;
          invR &= R2 > zero;
          R2 += one - invR;
          invR = rsqrt(R2);
          invR &= index < j;
          invR &= R2 > zero;

          mj *= invR * mi;
          pot += mj;
          B[j].TRG[0] += sum(mj);
          invR = invR * invR * mj;

          xj *= invR;
          ax += xj;
          B[j].TRG[1] -= sum(xj);

          yj *= invR;
          ay += yj;
          B[j].TRG[2] -= sum(yj);

          zj *= invR;
          az += zj;
          B[j].TRG[3] -= sum(zj);
        }
        for (int k=0; k<NSIMD; k++) {
          B[i+k].TRG[0] += transpose(pot, k);
          B[i+k].TRG[1] += transpose(ax, k);
          B[i+k].TRG[2] += transpose(ay, k);
          B[i+k].TRG[3] += transpose(az, k);
        }
      }
#endif
      for ( ; i<n; i++) {
        kreal_t pot = 0;
        kreal_t ax = 0;
        kreal_t ay = 0;
        kreal_t az = 0;
        for (int j=i+1; j<n; j++) {
          vec3 dX = B[j].X - B[i].X;
          real_t R2 = norm(dX) + eps2;
          if (R2 != 0) {
            real_t invR2 = 1.0 / R2;
            real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
            dX *= invR2 * invR;
            pot += invR;
            ax += dX[0];
            ay += dX[1];
            az += dX[2];
            B[j].TRG[0] += invR;
            B[j].TRG[1] -= dX[0];
            B[j].TRG[2] -= dX[1];
            B[j].TRG[3] -= dX[2];
          }
        }
        B[i].TRG[0] += pot;
        B[i].TRG[1] += ax;
        B[i].TRG[2] += ay;
        B[i].TRG[3] += az;
      }
    }

    static void P2M(C_iter C) {
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

    static void M2M(C_iter Ci, C_iter C0) {
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

    static void M2L(C_iter Ci, C_iter Cj, bool mutual) {
      complex_t Ynmi[P*P], Ynmj[P*P];
      vec3 dX = Ci->X - Cj->X - Xperiodic;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalLocal(P, rho, alpha, beta, Ynmi);
      if (mutual) evalLocal(P, rho, alpha+M_PI, beta, Ynmj);
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
              if (mutual) Lj += std::conj(Ci->M[nms]) * Cnm * Ynmj[jnkm];
            }
            for (int m=0; m<=n; m++) {
              int nms  = n * (n + 1) / 2 + m;
              int jnkm = (j + n) * (j + n) + j + n + m - k;
              real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
              Li += Cj->M[nms] * Cnm2 * Ynmi[jnkm];
              if (mutual) Lj += Ci->M[nms] * Cnm2 * Ynmj[jnkm];
            }
          }
          Ci->L[jks] += Li;
          if (mutual) Cj->L[jks] += Lj;
        }
      }
    }

    static void L2L(C_iter Ci, C_iter C0) {
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

    static void L2P(C_iter Ci) {
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
