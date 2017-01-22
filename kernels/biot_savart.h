#ifndef biot_savart_h
#define biot_savart_h
#include "spherical.h"
#if EXAFMM_USE_SIMD
#include "simdvec.h"
#endif

namespace exafmm {
  class Kernel {
  public:
    static const int P = Pmax;
    static vec3      Xperiodic;
    static real_t    eps2;
    static complex_t wavek;

    static void init() {}
    static void finalize() {}

    static void normalize(Bodies) {}

    static void P2P(C_iter Ci, C_iter Cj) {
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

    static void M2L(C_iter Ci, C_iter Cj) {
      complex_t Ynmi[P*P], Ynmj[P*P];
      vec3 dX = Ci->X - Cj->X - Xperiodic;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalLocal(P, rho, alpha, beta, Ynmi);
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
	      }
	    }
	    for (int m=0; m<=n; m++) {
	      int nms  = n * (n + 1) / 2 + m;
	      int jnkm = (j + n) * (j + n) + j + n + m - k;
	      real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
	      for (int d=0; d<3; d++) {
		Li[d] += Cj->M[3*nms+d] * Cnm2 * Ynmi[jnkm];
	      }
	    }
	  }
	  for (int d=0; d<3; d++) {
	    Ci->L[3*jks+d] += Li[d];
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
