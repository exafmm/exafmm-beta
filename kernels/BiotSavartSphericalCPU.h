#include "BiotSavartP2PCPU.h"

namespace exafmm {
    inline int oddOrEven(int n) {
      return (((n) & 1) == 1) ? -1 : 1;
    }

    inline int ipow2n(int n) {
      return (n >= 0) ? 1 : oddOrEven(n);
    }

    const complex_t I(0.,1.);                                   // Imaginary unit

//! Get r,theta,phi from x,y,z
    void cart2sph(real_t & r, real_t & theta, real_t & phi, vec3 dX) {
      r = sqrt(norm(dX));                                       // r = sqrt(x^2 + y^2 + z^2)
      theta = r == 0 ? 0 : acos(dX[2] / r);                     // theta = acos(z / r)
      phi = atan2(dX[1], dX[0]);                                // phi = atan(y / x)
    }

//! Spherical to cartesian coordinates
    void sph2cart(real_t r, real_t theta, real_t phi, vec3 spherical, vec3 & cartesian) {
      cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
	+ std::cos(theta) * std::cos(phi) / r * spherical[1]
	- std::sin(phi) / r / std::sin(theta) * spherical[2];
      cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
	+ std::cos(theta) * std::sin(phi) / r * spherical[1]
	+ std::cos(phi) / r / std::sin(theta) * spherical[2];
      cartesian[2] = std::cos(theta) * spherical[0]            // z component (not z itself)
	- std::sin(theta) / r * spherical[1];
    }

//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
    void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
      real_t x = std::cos(alpha);                               // x = cos(alpha)
      real_t y = std::sin(alpha);                               // y = sin(alpha)
      real_t invY = y == 0 ? 0 : 1 / y;                         // 1 / y
      real_t fact = 1;                                          // Initialize 2 * m + 1
      real_t pn = 1;                                            // Initialize Legendre polynomial Pn
      real_t rhom = 1;                                          // Initialize rho^m
      complex_t ei = std::exp(I * beta);                        // exp(i * beta)
      complex_t eim = 1.0;                                      // Initialize exp(i * m * beta)
      for (int m=0; m<P; m++) {                                 // Loop over m in Ynm
	real_t p = pn;                                          //  Associated Legendre polynomial Pnm
	int npn = m * m + 2 * m;                                //  Index of Ynm for m > 0
	int nmn = m * m;                                        //  Index of Ynm for m < 0
	Ynm[npn] = rhom * p * eim;                              //  rho^m * Ynm for m > 0
	Ynm[nmn] = std::conj(Ynm[npn]);                         //  Use conjugate relation for m < 0
	real_t p1 = p;                                          //  Pnm-1
	p = x * (2 * m + 1) * p1;                               //  Pnm using recurrence relation
	YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim; //  theta derivative of r^n * Ynm
	rhom *= rho;                                            //  rho^m
	real_t rhon = rhom;                                     //  rho^n
	for (int n=m+1; n<P; n++) {                             //  Loop over n in Ynm
	  int npm = n * n + n + m;                              //   Index of Ynm for m > 0
	  int nmm = n * n + n - m;                              //   Index of Ynm for m < 0
	  rhon /= -(n + m);                                     //   Update factorial
	  Ynm[npm] = rhon * p * eim;                            //   rho^n * Ynm
	  Ynm[nmm] = std::conj(Ynm[npm]);                       //   Use conjugate relation for m < 0
	  real_t p2 = p1;                                       //   Pnm-2
	  p1 = p;                                               //   Pnm-1
	  p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
	  YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;// theta derivative
	  rhon *= rho;                                          //   Update rho^n
	}                                                       //  End loop over n in Ynm
	rhom /= -(2 * m + 2) * (2 * m + 1);                     //  Update factorial
	pn = -pn * fact * y;                                    //  Pn
	fact += 2;                                              //  2 * m + 1
	eim *= ei;                                              //  Update exp(i * m * beta)
      }                                                         // End loop over m in Ynm
    }

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
    void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
      real_t x = std::cos(alpha);                               // x = cos(alpha)
      real_t y = std::sin(alpha);                               // y = sin(alpha)
      real_t fact = 1;                                          // Initialize 2 * m + 1
      real_t pn = 1;                                            // Initialize Legendre polynomial Pn
      real_t invR = -1.0 / rho;                                 // - 1 / rho
      real_t rhom = -invR;                                      // Initialize rho^(-m-1)
      complex_t ei = std::exp(I * beta);                        // exp(i * beta)
      complex_t eim = 1.0;                                      // Initialize exp(i * m * beta)
      for (int m=0; m<P; m++) {                                 // Loop over m in Ynm
	real_t p = pn;                                          //  Associated Legendre polynomial Pnm
	int npn = m * m + 2 * m;                                //  Index of Ynm for m > 0
	int nmn = m * m;                                        //  Index of Ynm for m < 0
	Ynm[npn] = rhom * p * eim;                              //  rho^(-m-1) * Ynm for m > 0
	Ynm[nmn] = std::conj(Ynm[npn]);                         //  Use conjugate relation for m < 0
	real_t p1 = p;                                          //  Pnm-1
	p = x * (2 * m + 1) * p1;                               //  Pnm using recurrence relation
	rhom *= invR;                                           //  rho^(-m-1)
	real_t rhon = rhom;                                     //  rho^(-n-1)
	for (int n=m+1; n<P; n++) {                             //  Loop over n in Ynm
	  int npm = n * n + n + m;                              //   Index of Ynm for m > 0
	  int nmm = n * n + n - m;                              //   Index of Ynm for m < 0
	  Ynm[npm] = rhon * p * eim;                            //   rho^n * Ynm for m > 0
	  Ynm[nmm] = std::conj(Ynm[npm]);                       //   Use conjugate relation for m < 0
	  real_t p2 = p1;                                       //   Pnm-2
	  p1 = p;                                               //   Pnm-1
	  p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
	  rhon *= invR * (n - m + 1);                           //   rho^(-n-1)
	}                                                       //  End loop over n in Ynm
	pn = -pn * fact * y;                                    //  Pn
	fact += 2;                                              //  2 * m + 1
	eim *= ei;                                              //  Update exp(i * m * beta)
      }                                                         // End loop over m in Ynm
    }

    class BiotSavartSphericalCPU : public BiotSavartP2PCPU {
    public:
    static void setup() {}

    static void P2M(C_iter C) {
      complex_t Ynm[P*P], YnmTheta[P*P];
      for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
	vec3 dX = B->X - C->X;
	real_t rho, alpha, beta;
	cart2sph(rho, alpha, beta, dX);
	evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
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
	cart2sph(rho, alpha, beta, dX);
	evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
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
      cart2sph(rho, alpha, beta, dX);
      evalLocal(rho, alpha, beta, Ynmi);
      if (mutual) evalLocal(rho, alpha+M_PI, beta, Ynmj);
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
      cart2sph(rho, alpha, beta, dX);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
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
	cart2sph(r, theta, phi, dX);
	evalMultipole(r, theta, phi, Ynm, YnmTheta);
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
