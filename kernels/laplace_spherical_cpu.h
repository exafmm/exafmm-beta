#ifndef laplace_spherical_cpu_h
#define laplace_spherical_cpu_h
#include "laplace_p2p_cpu.h"
#include "spherical.h"

namespace exafmm {
  template<int _P>
  class LaplaceSphericalCPU : public LaplaceP2PCPU<vec<_P*(_P+1)/2,complex_t>,Spherical> {
  public:
    static const Basis basis = Spherical;                       //!< Set basis to Spherical
    static const int P = _P;                                    //!< Set order of expansion
    static const int NTERM = P*(P+1)/2;                         //!< # of terms in Laplace Spherical expansion
    typedef vec<NTERM,complex_t> vecP;                          //!< Vector type for expansion terms
    using typename LaplaceP2PCPU<vecP,Spherical>::Bodies;       //!< Vector of bodies for Laplace
    using typename LaplaceP2PCPU<vecP,Spherical>::B_iter;       //!< Iterator of body vector
    using typename LaplaceP2PCPU<vecP,Spherical>::Cells;        //!< Vector of cells for Laplace
    using typename LaplaceP2PCPU<vecP,Spherical>::C_iter;       //!< Iterator of cell vector
    using LaplaceP2PCPU<vecP,Spherical>::Xperiodic;

    static void init() {}
    static void finalize() {}

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
