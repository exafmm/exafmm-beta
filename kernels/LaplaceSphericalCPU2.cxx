#include "kernel.h"

#define SIGN(n) ((n >= 0) - (n < 0))
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const complex_t I(0.,1.);                                       // Imaginary unit
const real_t EPS = 1e-6;                                        // Epsilon
real_t *factorial;                                              // Factorial
real_t *prefactor;                                              // \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
real_t *Anm;                                                    // \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
complex_t *Cnm;                                                 // M2L translation matrix \f$ C_{jn}^{km} \f$

//! Precalculate M2L translation matrix
void preTraversal() {
  factorial = new real_t [P];                                   // Factorial
  prefactor = new real_t [4*P*P];                               // sqrt( (n - |m|)! / (n + |m|)! )
  Anm       = new real_t [4*P*P];                               // (-1)^n / sqrt( (n + m)! / (n - m)! )
  Cnm       = new complex_t [P*P*P*P];                          // M2L translation matrix Cjknm

  factorial[0] = 1;                                             // Initialize factorial
  for (int n=1; n<P; n++) {                                     // Loop to P
    factorial[n] = factorial[n-1] * n;                          //  n!
  }                                                             // End loop to P

  for (int n=0; n<2*P; n++) {                                   // Loop over n in Anm
    for (int m=-n; m<=n; m++) {                                 //  Loop over m in Anm
      int nm = n*n+n+m;                                         //   Index of Anm
      int nabsm = abs(m);                                       //   |m|
      real_t fnmm = EPS;                                        //   Initialize (n - m)!
      for (int i=1; i<=n-m; i++) fnmm *= i;                     //   (n - m)!
      real_t fnpm = EPS;                                        //   Initialize (n + m)!
      for (int i=1; i<=n+m; i++) fnpm *= i;                     //   (n + m)!
      real_t fnma = 1.0;                                        //   Initialize (n - |m|)!
      for (int i=1; i<=n-nabsm; i++) fnma *= i;                 //   (n - |m|)!
      real_t fnpa = 1.0;                                        //   Initialize (n + |m|)!
      for (int i=1; i<=n+nabsm; i++) fnpa *= i;                 //   (n + |m|)!
      prefactor[nm] = std::sqrt(fnma/fnpa);                     //   sqrt( (n - |m|)! / (n + |m|)! )
      Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);                //   (-1)^n / sqrt( (n + m)! / (n - m)! )
    }                                                           //  End loop over m in Anm
  }                                                             // End loop over n in Anm

  for (int j=0, jk=0, jknm=0; j<P; j++) {                       // Loop over j in Cjknm
    for (int k=-j; k<=j; k++, jk++){                            //  Loop over k in Cjknm
      for (int n=0, nm=0; n<P; n++) {                           //   Loop over n in Cjknm
	for (int m=-n; m<=n; m++, nm++, jknm++) {               //    Loop over m in Cjknm
	  const int jnkm = (j+n)*(j+n)+j+n+m-k;                 //     Index C_{j+n}^{m-k}
	  Cnm[jknm] = std::pow(I,real_t(abs(k-m)-abs(k)-abs(m)))//     Cjknm
	    * real_t(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
	}                                                       //    End loop over m in Cjknm
      }                                                         //   End loop over n in Cjknm
    }                                                           //  End loop over in k in Cjknm
  }                                                             // End loop over in j in Cjknm
}

//! Free temporary allocations
void postTraversal() {
  delete[] factorial;                                           // Free factorial
  delete[] prefactor;                                           // Free sqrt( (n - |m|)! / (n + |m|)! )
  delete[] Anm;                                                 // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
  delete[] Cnm;                                                 // Free M2L translation matrix Cjknm
}

//! Get r,theta,phi from x,y,z
void cart2sph(real_t& r, real_t& theta, real_t& phi, vec3 dX) {
  r = sqrt(norm(dX)) * 1.000001;                                // r = sqrt(x^2 + y^2 + z^2)
  theta = acos(dX[2] / r);                                      // theta = acos(z / r)
  phi = atan2(dX[1],dX[0]);                                     // phi = atan(y / x)
}

//! Spherical to cartesian coordinates
template<typename T>
void sph2cart(real_t r, real_t theta, real_t phi, T spherical, T &cartesian) {
  cartesian[0] = sin(theta) * cos(phi) * spherical[0]           // x component (not x itself)
               + cos(theta) * cos(phi) / r * spherical[1]
               - sin(phi) / r / sin(theta) * spherical[2];
  cartesian[1] = sin(theta) * sin(phi) * spherical[0]           // y component (not y itself)
               + cos(theta) * sin(phi) / r * spherical[1]
               + cos(phi) / r / sin(theta) * spherical[2];
  cartesian[2] = cos(theta) * spherical[0]                      // z component (not z itself)
               - sin(theta) / r * spherical[1];
}

//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t *Ynm, complex_t *YnmTheta) {
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = 1;                                              // Initialize rho^m
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    complex_t eim = std::exp(I * real_t(m * beta));             //  exp(i * m * beta)
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * prefactor[npn] * eim;                 //  rho^m * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
    rhom *= rho;                                                //  rho^m
    real_t rhon = rhom;                                         //  rho^n
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * prefactor[npm] * eim;               //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
      rhon *= rho;                                              //   Update rho^n
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
  }                                                             // End loop over m in Ynm
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t *Ynm) {
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
  for (int m=0; m<2*P; m++) {                                   // Loop over m in Ynm
    complex_t eim = std::exp(I * real_t(m * beta));             //  exp(i * m * beta)
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * prefactor[npn] * eim;                 //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    rhom /= rho;                                                //  rho^(-m-1)
    real_t rhon = rhom;                                         //  rho^(-n-1)
    for (int n=m+1; n<2*P; n++) {                               //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * prefactor[npm] * eim;               //   rho^n * Ynm for m > 0
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon /= rho;                                              //   rho^(-n-1)
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
  }                                                             // End loop over m in Ynm
}

void Kernel::P2M(C_iter C) const {
  complex_t Ynm[4*P*P], YnmTheta[4*P*P];
  for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {
    vec3 dX = B->X - C->X;
    real_t R = std::sqrt(norm(dX));
    if (R > C->RMAX) C->RMAX = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        int nm  = n * n + n + m;
        int nms = n * (n + 1) / 2 + m;
        C->M[nms] += B->SRC * Ynm[nm];
      }
    }
  }
#if USE_RMAX
  C->RCRIT = std::min(C->R,C->RMAX);
#else
  C->RCRIT = C->R;
#endif
}

void Kernel::M2M(C_iter Ci, C_iter C0) const {
  complex_t Ynm[4*P*P], YnmTheta[4*P*P];
  for (C_iter Cj=C0+Ci->CHILD; Cj!=C0+Ci->CHILD+Ci->NCHILD; Cj++) {
    vec3 dX = Ci->X - Cj->X;
    real_t R = std::sqrt(norm(dX)) + Cj->RCRIT;
    if (R > Ci->RMAX) Ci->RMAX = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for (int j=0; j<P; j++) {
      for (int k=0; k<=j; k++) {
	int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        complex_t M = 0;
        for (int n=0; n<=j; n++) {
          for (int m=-n; m<=std::min(k-1,n); m++) {
            if( j-n >= k-m ) {
	      int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n + m;
              M += Cj->M[jnkms] * std::pow(I,real_t(m-abs(m))) * Ynm[nm]
	        * real_t(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
	    }
          }
          for (int m=k; m<=n; m++) {
	    if( j-n >= m-k ) {
  	      int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n + m;
	      M += std::conj(Cj->M[jnkms]) * Ynm[nm]
	        * real_t(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
	    }
          }
        }
        Ci->M[jks] += M * EPS;
      }
    }
  }
#if USE_RMAX
  Ci->RCRIT = std::min(Ci->R,Ci->RMAX);
#else
  Ci->RCRIT = Ci->R;
#endif
}

void Kernel::M2L(C_iter Ci, C_iter Cj, bool mutual) const {
  complex_t Ynmi[4*P*P], Ynmj[4*P*P];
  vec3 dX = Ci->X - Cj->X - Xperiodic;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalLocal(rho,alpha,beta,Ynmi);
  if (mutual) evalLocal(rho,alpha+M_PI,beta,Ynmj);
  for (int j=0; j<P; j++) {
    for (int k=0; k<=j; k++) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex_t Li = 0, Lj = 0;
      for (int n=0; n<P; n++) {
        for (int m=-n; m<0; m++) {
	  int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 - m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          Li += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynmi[jnkm];
          if (mutual) Lj += std::conj(Ci->M[nms]) * Cnm[jknm] * Ynmj[jnkm];
        }
        for (int m=0; m<=n; m++) {
	  int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 + m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          Li += Cj->M[nms] * Cnm[jknm] * Ynmi[jnkm];
          if (mutual) Lj += Ci->M[nms] * Cnm[jknm] * Ynmj[jnkm];
        }
      }
      Ci->L[jks] += Li;
      if (mutual) Cj->L[jks] += Lj;
    }
  }
}

void Kernel::L2L(C_iter Ci, C_iter C0) const {
  complex_t Ynm[4*P*P], YnmTheta[4*P*P];
  C_iter Cj = C0 + Ci->PARENT;
  vec3 dX = Ci->X - Cj->X;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
  for (int j=0; j<P; j++) {
    for (int k=0; k<=j; k++) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for (int n=j; n<P; n++) {
        for (int m=j+k-n; m<0; m++) {
          int jnkm = (n - j) * (n - j) + n - j + m - k;
	  int nm   = n * n + n - m;
          int nms  = n * (n + 1) / 2 - m;
	  L += std::conj(Cj->L[nms]) * Ynm[jnkm]
	    * real_t(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
        }
        for (int m=0; m<=n; m++) {
          if( n-j >= abs(m-k) ) {
	    int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nm   = n * n + n + m;
            int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[nms] * std::pow(I,real_t(m-k-abs(m-k)))
	      * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
          }
        }
      }
      Ci->L[jks] += L * EPS;
    }
  }
}

void Kernel::L2P(C_iter Ci) const {
  complex_t Ynm[4*P*P], YnmTheta[4*P*P];
  for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {
    vec3 dX = B->X - Ci->X;
    vec3 spherical = 0;
    vec3 cartesian = 0;
    real_t r, theta, phi;
    cart2sph(r,theta,phi,dX);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);
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
    sph2cart(r,theta,phi,spherical,cartesian);
    B->TRG[1] += cartesian[0];
    B->TRG[2] += cartesian[1];
    B->TRG[3] += cartesian[2];
  }
}
