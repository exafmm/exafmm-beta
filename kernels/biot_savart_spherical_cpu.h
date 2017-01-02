#ifndef biot_savart_spherical_cpu_h
#define biot_savart_spherical_cpu_h
#include "biot_savart_p2p_cpu.h"
#include "spherical.h"

namespace exafmm {
  template<int _P>
  class BiotSavartSphericalCPU : public BiotSavartP2PCPU<vec<3*_P*(_P+1)/2,complex_t>,Spherical> {
  public:
    static const Basis basis = Spherical;                       //!< Set basis to Spherical
    static const int P = _P;                                    //!< Set order of expansion
    static const int NTERM = 3*P*(P+1)/2;                       //!< # of terms in Biot-Savart Spherical expansion
    typedef vec<NTERM,complex_t> vecP;                          //!< Vector type for expansion terms
    using typename BiotSavartP2PCPU<vecP,Spherical>::Bodies;    //!< Vector of bodies for Biot-Savart
    using typename BiotSavartP2PCPU<vecP,Spherical>::B_iter;    //!< Iterator of body vector
    using typename BiotSavartP2PCPU<vecP,Spherical>::Cells;     //!< Vector of cells for Biot-Savart
    using typename BiotSavartP2PCPU<vecP,Spherical>::C_iter;    //!< Iterator of cell vector
    using BiotSavartP2PCPU<vecP,Spherical>::Xperiodic;

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
	    for (int d=0; d<3; d++) {
	      C->M[3*nms+d] += B->SRC[d] * Ynm[nm];
	    }
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
	    complex_t M[3] = {0, 0, 0};
	    for (int n=0; n<=j; n++) {
	      for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
		int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
		int nm    = n * n + n - m;
		for (int d=0; d<3; d++) {
		  M[d] += Cj->M[3*jnkms+d] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
		}
	      }
	      for (int m=k; m<=std::min(n,j+k-n); m++) {
		int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
		int nm    = n * n + n - m;
		for (int d=0; d<3; d++) {
		  M[d] += std::conj(Cj->M[3*jnkms+d]) * Ynm[nm] * real_t(oddOrEven(k+n+m));
		}
	      }
	    }
	    for (int d=0; d<3; d++) {
	      Ci->M[3*jks+d] += M[d];
	    }
	  }
	}
      }
    }

    static void M2L(C_iter Ci, C_iter Cj, bool mutual) {
      assert(mutual == false);
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
	  complex_t Li[3] = {0, 0, 0}, Lj[3] = {0, 0, 0};
	  for (int n=0; n<P-j; n++) {
	    for (int m=-n; m<0; m++) {
	      int nms  = n * (n + 1) / 2 - m;
	      int jnkm = (j + n) * (j + n) + j + n + m - k;
	      for (int d=0; d<3; d++) {
		Li[d] += std::conj(Cj->M[3*nms+d]) * Cnm * Ynmi[jnkm];
		if (mutual) Lj[d] += std::conj(Ci->M[3*nms+d]) * Cnm * Ynmj[jnkm];
	      }
	    }
	    for (int m=0; m<=n; m++) {
	      int nms  = n * (n + 1) / 2 + m;
	      int jnkm = (j + n) * (j + n) + j + n + m - k;
	      real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
	      for (int d=0; d<3; d++) {
		Li[d] += Cj->M[3*nms+d] * Cnm2 * Ynmi[jnkm];
		if (mutual) Lj[d] += Ci->M[3*nms+d] * Cnm2 * Ynmj[jnkm];
	      }
	    }
	  }
	  for (int d=0; d<3; d++) {
	    Ci->L[3*jks+d] += Li[d];
	    if (mutual) Cj->L[3*jks+d] += Lj[d];
	  }
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
	  complex_t L[3] = {0, 0, 0};
	  for (int n=j; n<P; n++) {
	    for (int m=j+k-n; m<0; m++) {
	      int jnkm = (n - j) * (n - j) + n - j + m - k;
	      int nms  = n * (n + 1) / 2 - m;
	      for (int d=0; d<3; d++) {
		L[d] += std::conj(Cj->L[3*nms+d]) * Ynm[jnkm] * real_t(oddOrEven(k));
	      }
	    }
	    for (int m=0; m<=n; m++) {
	      if( n-j >= abs(m-k) ) {
		int jnkm = (n - j) * (n - j) + n - j + m - k;
		int nms  = n * (n + 1) / 2 + m;
		for (int d=0; d<3; d++) {
		  L[d] += Cj->L[3*nms+d] * Ynm[jnkm] * real_t(oddOrEven((m-k)*(m<k)));
		}
	      }
	    }
	  }
	  for (int d=0; d<3; d++) {
	    Ci->L[3*jks+d] += L[d];
	  }
	}
      }
    }

    static void L2P(C_iter Ci) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NBODY; B++) {
	vec3 dX = B->X - Ci->X + EPS;
	vec3 spherical[3] = {0, 0, 0};
	vec3 cartesian[3] = {0, 0, 0};
	real_t r, theta, phi;
	cart2sph(dX, r, theta, phi);
	evalMultipole(P, r, theta, phi, Ynm, YnmTheta);
	for (int n=0; n<P; n++) {
	  int nm  = n * n + n;
	  int nms = n * (n + 1) / 2;
	  B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
	  for (int d=0; d<3; d++) {
	    spherical[d][0] += std::real(Ci->L[3*nms+d] * Ynm[nm]) / r * n;
	    spherical[d][1] += std::real(Ci->L[3*nms+d] * YnmTheta[nm]);
	  }
	  for( int m=1; m<=n; m++) {
	    nm  = n * n + n + m;
	    nms = n * (n + 1) / 2 + m;
	    for (int d=0; d<3; d++) {
	      spherical[d][0] += 2 * std::real(Ci->L[3*nms+d] * Ynm[nm]) / r * n;
	      spherical[d][1] += 2 * std::real(Ci->L[3*nms+d] * YnmTheta[nm]);
	      spherical[d][2] += 2 * std::real(Ci->L[3*nms+d] * Ynm[nm] * I) * m;
	    }
	  }
	}
	for (int d=0; d<3; d++) {
	  sph2cart(r, theta, phi, spherical[d], cartesian[d]);
	}
	B->TRG[0] = 1;
	B->TRG[1] += cartesian[1][2] - cartesian[2][1];
	B->TRG[2] += cartesian[2][0] - cartesian[0][2];
	B->TRG[3] += cartesian[0][1] - cartesian[1][0];
      }
    }
  };
}
#endif
