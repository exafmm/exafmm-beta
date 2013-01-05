#include "kernel.h"
#include <iostream>

#define SIGN(n) ((n >= 0) - (n < 0))
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const complex_t I(0.,1.);

template<int j, int k, int n, int m>
struct M2LTemplate {
  static const int nms = n * (n + 1) / 2 + abs(m);
  static const int jnkm = (j + n) * (j + n) + j + n + m - k;
  static const int Cnm = ODDEVEN(j);
  static const int Cnm2 = Cnm * ODDEVEN((k-m)*(k<m)+m);
  static const int mm1 = (m-1) < -n ? 0 : m-1;
  static inline void negative(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k,n,mm1>::negative(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    if (m >= -n) {
      Li2 += std::conj(Mj[nms]) * Ynmi[jnkm] * real_t(Cnm);
      if (mutual) Lj2 += std::conj(Mi[nms]) * Ynmj[jnkm] * real_t(Cnm);
    }
  }
  static inline void positive(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k,n,m-1>::positive(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm2);
    if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm2);
  }
};

template<int j, int k, int n>
struct M2LTemplate<j,k,n,0> {
  static const int nms = n * (n + 1) / 2;
  static const int jnkm = (j + n) * (j + n) + j + n - k;
  static const int Cnm = ODDEVEN(j);
  static const int Cnm2 = Cnm * ODDEVEN(k*(k<0));
  static inline void negative(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k,n-1,n-1>::positive(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
  }
  static inline void positive(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k,n,-1>::negative(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm2);
    if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm2);
  }
};

template<int j, int k>
struct M2LTemplate<j,k,0,0> {
  static const int nms = 0;
  static const int jks = j * (j + 1) / 2 + k - 1;
  static const int jnkm = j * j + j - k;
  static const int Cnm = ODDEVEN(j);
  static const int Cnm2 = Cnm * ODDEVEN(k*(k<0));
  static inline void negative(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k-1,P-j-1,P-j-1>::positive(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li[jks] = Li2;
    Lj[jks] = Lj2;
    Li2 = Lj2 = 0;
  }
  static inline void positive(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,k,0,-1>::negative(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm2);
    if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm2);
  }
};

template<int j>
struct M2LTemplate<j,0,0,0> {
  static const int nms = 0;
  static const int jks = j * (j - 1) / 2 + j - 1;
  static const int jnkm = j * j + j;
  static const int Cnm = ODDEVEN(j);
  static inline void negative(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j-1,j-1,P-j,P-j>::positive(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li[jks] = Li2;
    Lj[jks] = Lj2;
    Li2 = Lj2 = 0;
  }
  static inline void positive(const vecM &Mi, const vecM &Mj, vecL &Li, vecL &Lj, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    M2LTemplate<j,0,0,-1>::negative(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
    Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm);
    if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm);
  }
};

template<>
struct M2LTemplate<0,0,0,0> {
  static const int nms = 0;
  static const int jnkm = 0;
  static const int Cnm = 1;
  static inline void negative(const vecM&, const vecM&, vecL&, vecL&, complex_t, complex_t, complex_t*, complex_t*, bool) {}
  static inline void positive(const vecM &Mi, const vecM &Mj, vecL&, vecL&, complex_t &Li2, complex_t &Lj2,
                              complex_t *Ynmi, complex_t *Ynmj, bool mutual) {
    Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm);
    if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm);
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
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t *Ynm) {
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t invR = -1.0 / rho;                                     // - 1 / rho
  complex_t ei = std::exp(I * beta);                            // exp(i * beta)
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
  vecM Mi = Ci->M;
  vecM Mj = Cj->M;
  complex_t Li2 = 0, Lj2 = 0;
#if 1
  vecL Li = I * real_t(0), Lj = I * real_t(0);
  M2LTemplate<P-1,P-1,0,0>::positive(Mi,Mj,Li,Lj,Li2,Lj2,Ynmi,Ynmj,mutual);
  int jks = (P + 2) * (P - 1) / 2;
  Li[jks] = Li2;
  Lj[jks] = Lj2;
  Ci->L += Li;
  Cj->L += Lj;
#else
  for (int j=0; j<P; j++) {
    int Cnm = ODDEVEN(j);
    for (int k=0; k<=j; k++) {
      int jks = j * (j + 1) / 2 + k;
      Li2 = Lj2 = 0;
      for (int n=0; n<P-j; n++) {
        for (int m=-n; m<0; m++) {
          int nms  = n * (n + 1) / 2 - m;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          Li2 += std::conj(Mj[nms]) * Ynmi[jnkm] * real_t(Cnm);
          if (mutual) Lj2 += std::conj(Mi[nms]) * Ynmj[jnkm] * real_t(Cnm);
        }
        for (int m=0; m<=n; m++) {
          int nms  = n * (n + 1) / 2 + m;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          int Cnm2 = Cnm * ODDEVEN((k-m)*(k<m)+m);
          Li2 += Mj[nms] * Ynmi[jnkm] * real_t(Cnm2);
          if (mutual) Lj2 += Mi[nms] * Ynmj[jnkm] * real_t(Cnm2);
        }
      }
      Ci->L[jks] += Li2;
      if (mutual) Cj->L[jks] += Lj2;
    }
  }
#endif
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
