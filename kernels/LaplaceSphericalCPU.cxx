#include "kernel.h"
#include <iostream>

#define SIGN(n) ((n >= 0) - (n < 0))
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const complex_t I(0.,1.);

template<int m, int n>
struct Harmonics {
  static const int c0 = - n - 2 * m;
  static const int c1 = 2 * n + 2 * m - 1;
  static const int c2 = - n - 2 * m + 1;
  static const int c3 = n + 1;
  static const int c4 =  - n - m - 1;
  static const int npm = (n + m) * (n + m) + n + 2 * m;
  static const int nmm = (n + m) * (n + m) + n;
  static inline real_t Rn(const real_t &R) {
    return Harmonics<m,n-1>::Rn(R) * R / c0;
  }
  static inline real_t invRn(const real_t &invR) {
    return Harmonics<m,n-1>::invRn(invR) * invR * n;
  }
  static inline real_t Pnm(const real_t &x, const real_t &y) {
    return (x * c1 * Harmonics<m,n-1>::Pnm(x,y) + c2 * Harmonics<m,n-2>::Pnm(x,y)) / n;
  }
  static inline complex_t Em(const complex_t &ei) {
    return Harmonics<m,0>::Em(ei);
  }
  static inline void Multipole(complex_t *Ynm, complex_t *YnmTheta, const real_t &R, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m,n-1>::Multipole(Ynm,YnmTheta,R,x,y,ei);
    Ynm[npm] = Rn(R) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
    YnmTheta[npm] = Rn(R) * (c3 * Harmonics<m,n+1>::Pnm(x,y) + c4 * x * Pnm(x,y)) / y * Em(ei);
  }
  static inline void Local(complex_t *Ynm, const real_t &invR, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m,n-1>::Local(Ynm,invR,x,y,ei);
    Ynm[npm] = invRn(invR) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
  }
};

template<int m>
struct Harmonics<m,1> {
  static const int c0 = 2 * m + 1;
  static const int c1 = - c0;
  static const int c2 =  - m - 2;
  static const int npm = (m + 1) * (m + 1) + 2 * m + 1;
  static const int nmm = (m + 1) * (m + 1) + 1;
  static inline real_t Rn(const real_t &R) {
    return Harmonics<m,0>::Rn(R) * R / c1;
  }
  static inline real_t invRn(const real_t &invR) {
    return Harmonics<m,0>::invRn(invR) * invR;
  }
  static inline real_t Pnm(const real_t &x, const real_t &y) {
    return x * c0 * Harmonics<m,0>::Pnm(x,y);
  }
  static inline complex_t Em(const complex_t &ei) {
    return Harmonics<m,0>::Em(ei);
  }
  static inline void Multipole(complex_t *Ynm, complex_t *YnmTheta, const real_t &R, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m,0>::Multipole(Ynm,YnmTheta,R,x,y,ei);
    Ynm[npm] = Rn(R) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
    YnmTheta[npm] = Rn(R) * (2 * Harmonics<m,2>::Pnm(x,y) + c2 * x * Pnm(x,y)) / y * Em(ei);
  }
  static inline void Local(complex_t *Ynm, const real_t &invR, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m,0>::Local(Ynm,invR,x,y,ei);
    Ynm[npm] = invRn(invR) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
  }
};

template<int m>
struct Harmonics<m,0> {
  static const int c0 = 2 * m - 1;
  static const int c1 = - 2 * m * c0;
  static const int c2 =  - m - 1;
  static const int npm = m * m + 2 * m;
  static const int nmm = m * m;
  static inline real_t Rn(const real_t &R) {
    return Harmonics<m-1,0>::Rn(R) * R / c1;
  }
  static inline real_t invRn(const real_t &invR) {
    return Harmonics<m-1,0>::invRn(invR) * invR;
  }
  static inline real_t Pnm(const real_t &x, const real_t &y) {
    return - c0 * y * Harmonics<m-1,0>::Pnm(x,y);
  }
  static inline complex_t Em(const complex_t &ei) {
    return Harmonics<m-1,0>::Em(ei) * ei;
  }
  static inline void Multipole(complex_t *Ynm, complex_t *YnmTheta, const real_t &R, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m-1,P-m>::Multipole(Ynm,YnmTheta,R,x,y,ei);
    Ynm[npm] = Rn(R) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
    YnmTheta[npm] = Rn(R) * (Harmonics<m,1>::Pnm(x,y) + c2 * x * Pnm(x,y)) / y * Em(ei);
  }
  static inline void Local(complex_t *Ynm, const real_t &invR, const real_t &x, const real_t &y, const complex_t &ei) {
    Harmonics<m-1,P-m>::Local(Ynm,invR,x,y,ei);
    Ynm[npm] = invRn(invR) * Pnm(x,y) * Em(ei);
    Ynm[nmm] = std::conj(Ynm[npm]);
  }
};

template<>
struct Harmonics<0,0> {
  static inline real_t Rn(const real_t &) {
    return 1.0;
  }
  static inline real_t invRn(const real_t &invR) {
    return -invR;
  }
  static inline real_t Pnm(const real_t &, const real_t &) {
    return 1.0;
  }
  static inline complex_t Em(const complex_t &) {
    return 1.0;
  }
  static inline void Multipole(complex_t *Ynm, complex_t *YnmTheta, const real_t&, const real_t&, const real_t&, const complex_t&) {
    Ynm[0] = 1.0;
    YnmTheta[0] = 0.0;
  }
  static inline void Local(complex_t *Ynm, const real_t &invR, const real_t&, const real_t&, const complex_t&) {
    Ynm[0] = -invR;
  }
};

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
  complex_t ei = std::exp(I * beta);                            // exp(i * beta)
#if 1
  Harmonics<P-1,0>::Multipole(Ynm,YnmTheta,rho,x,y,ei);
#else
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = 1;                                              // Initialize rho^m
  complex_t eim = 1.0;                                          // Initialize exp(i * m * beta)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * eim;                                  //  rho^m * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * eim;    //  Theta derivative of r^n * Ynm
    rhom *= rho;                                                //  rho^m
    real_t rhon = rhom;                                         //  rho^n
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      rhon /= -(n + m);                                         //   Update factorial
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * eim;// theta derivative
      rhon *= rho;                                              //   Update rho^n
    }                                                           //  End loop over n in Ynm
    rhom /= -(2 * m + 2) * (2 * m + 1);                         //  Update factorial
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
    eim *= ei;                                                  //  Update exp(i * m * beta)
  }                                                             // End loop over m in Ynm
#endif
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t *Ynm) {
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t invR = -1.0 / rho;                                     // - 1 / rho
  complex_t ei = std::exp(I * beta);                            // exp(i * beta)
#if 1
  Harmonics<P-1,0>::Local(Ynm,invR,x,y,ei);
#else
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = -invR;                                          // Initialize rho^(-m-1)
  complex_t eim = 1.0;                                          // Initialize exp(i * m * beta)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * eim;                                  //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    rhom *= invR;                                               //  rho^(-m-1)
    real_t rhon = rhom;                                         //  rho^(-n-1)
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm for m > 0
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon *= invR * (n - m + 1);                               //   rho^(-n-1)
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
    eim *= ei;                                                  //  Update exp(i * m * beta)
  }                                                             // End loop over m in Ynm
#endif
}

void Kernel::P2M(C_iter C, real_t &Rmax) const {
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {
    vec3 dX = B->X - C->X;
    real_t R = std::sqrt(norm(dX));
    if (R > Rmax) Rmax = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
        int nm  = n * n + n - m;
        int nms = n * (n + 1) / 2 + m;
        C->M[nms] += B->SRC * Ynm[nm];
      }
    }
  }
#if USE_RMAX
  C->RCRIT = std::min(C->R,Rmax);
#else
  C->RCRIT = C->R;
#endif
}

void Kernel::M2M(C_iter Ci, real_t &Rmax) const {
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; Cj++) {
    vec3 dX = Ci->X - Cj->X;
    real_t R = std::sqrt(norm(dX)) + Cj->RCRIT;
    if (R > Rmax) Rmax = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
    for (int j=0; j<P; j++) {
      for (int k=0; k<=j; k++) {
        int jks = j * (j + 1) / 2 + k;
        complex_t M = 0;
        for (int n=0; n<=j; n++) {
          for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
            int nm    = n * n + n - m;
            M += Cj->M[jnkms] * Ynm[nm] * real_t(SIGN(m) * ODDEVEN(n));
          }
          for (int m=k; m<=std::min(n,j+k-n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            int nm    = n * n + n - m;
            M += std::conj(Cj->M[jnkms]) * Ynm[nm] * real_t(ODDEVEN(k+n+m));
          }
        }
        Ci->M[jks] += M;
      }
    }
  }
#if USE_RMAX
  Ci->RCRIT = std::min(Ci->R,Rmax);
#else
  Ci->RCRIT = Ci->R;
#endif
}

void Kernel::M2L(C_iter Ci, C_iter Cj, bool mutual) const {
  complex_t Ynmi[P*P], Ynmj[P*P];
  vec3 dX = Ci->X - Cj->X - Xperiodic;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalLocal(rho,alpha,beta,Ynmi);
  if (mutual) evalLocal(rho,alpha+M_PI,beta,Ynmj);
  for (int j=0; j<P; j++) {
    real_t Cnm = ODDEVEN(j);
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
          real_t Cnm2 = Cnm * ODDEVEN((k-m)*(k<m)+m);
          Li += Cj->M[nms] * Cnm2 * Ynmi[jnkm];
          if (mutual) Lj += Ci->M[nms] * Cnm2 * Ynmj[jnkm];
        }
      }
      Ci->L[jks] += Li;
      if (mutual) Cj->L[jks] += Lj;
    }
  }
}

void Kernel::L2L(C_iter Ci) const {
  complex_t Ynm[P*P], YnmTheta[P*P];
  C_iter Cj = Ci0 + Ci->PARENT;
  vec3 dX = Ci->X - Cj->X;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
  for (int j=0; j<P; j++) {
    for (int k=0; k<=j; k++) {
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for (int n=j; n<P; n++) {
        for (int m=j+k-n; m<0; m++) {
          int jnkm = (n - j) * (n - j) + n - j + m - k;
          int nms  = n * (n + 1) / 2 - m;
          L += std::conj(Cj->L[nms]) * Ynm[jnkm] * real_t(ODDEVEN(k));
        }
        for (int m=0; m<=n; m++) {
          if( n-j >= abs(m-k) ) {
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[nms] * Ynm[jnkm] * real_t(ODDEVEN((m-k)*(m<k)));
          }
        }
      }
      Ci->L[jks] += L;
    }
  }
}

void Kernel::L2P(C_iter Ci) const {
  complex_t Ynm[P*P], YnmTheta[P*P];
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
