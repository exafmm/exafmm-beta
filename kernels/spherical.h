#ifndef spherical_h
#define spherical_h
#include "types.h"

namespace exafmm {
  inline int oddOrEven(int n) {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  inline int ipow2n(int n) {
    return (n >= 0) ? 1 : oddOrEven(n);
  }

  //! Get r,theta,phi from x,y,z
  void cart2sph(vec3 dX, real_t & r, real_t & theta, real_t & phi) {
    r = sqrt(norm(dX));                                         // r = sqrt(x^2 + y^2 + z^2)
    theta = r == 0 ? 0 : acos(dX[2] / r);                       // theta = acos(z / r)
    phi = atan2(dX[1], dX[0]);                                  // phi = atan(y / x)
  }

  //! Spherical to cartesian coordinates
  void sph2cart(real_t r, real_t theta, real_t phi, vec3 spherical, vec3 & cartesian) {
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
      + std::cos(theta) * std::cos(phi) / r * spherical[1]
      - std::sin(phi) / r / std::sin(theta) * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
      + std::cos(theta) * std::sin(phi) / r * spherical[1]
      + std::cos(phi) / r / std::sin(theta) * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0]               // z component (not z itself)
      - std::sin(theta) / r * spherical[1];
  }

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(int P, real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t invY = y == 0 ? 0 : 1 / y;                           // 1 / y
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t rhom = 1;                                            // Initialize rho^m
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim; //  theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real_t rhon = rhom;                                       //  rho^n
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        rhon /= -(n + m);                                       //   Update factorial
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      rhom /= -(2 * m + 2) * (2 * m + 1);                       //  Update factorial
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(int P, real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t invR = -1.0 / rho;                                   // - 1 / rho
    real_t rhom = -invR;                                        // Initialize rho^(-m-1)
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      rhom *= invR;                                             //  rho^(-m-1)
      real_t rhon = rhom;                                       //  rho^(-n-1)
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        rhon *= invR * (n - m + 1);                             //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }
}
#endif
