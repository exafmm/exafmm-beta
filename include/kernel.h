/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef kernel_h
#define kernel_h
#include "sort.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;                                          //!< P^2
const int  P4 = P2 * P2;                                        //!< P^4

//! Unified CPU/GPU kernel class
class KernelBase : public Sort {
protected:
  vect                 X0;                                      //!< Center of root cell
  real                 R0;                                      //!< Radius of root cell
  C_iter               Ci0;                                     //!< icells.begin()
  C_iter               Cj0;                                     //!< jcells.begin()

  int                  ATOMS;                                   //!< Number of atom types in Van der Waals
  std::vector<real>    RSCALE;                                  //!< Scaling parameter for Van der Waals
  std::vector<real>    GSCALE;                                  //!< Scaling parameter for Van der Waals
  real                 KSIZE;                                   //!< Number of waves in Ewald summation
  real                 ALPHA;                                   //!< Scaling parameter for Ewald summation
  real                 SIGMA;                                   //!< Scaling parameter for Ewald summation

  std::vector<int>     keysHost;                                //!< Offsets for rangeHost
  std::vector<int>     rangeHost;                               //!< Offsets for sourceHost
  std::vector<gpureal> sourceHost;                              //!< Sources on host
  std::vector<gpureal> targetHost;                              //!< Targets on host
  std::vector<gpureal> constHost;                               //!< Constants on host
  Map                  sourceBegin;                             //!< Define map for offset of source cells
  Map                  sourceSize;                              //!< Define map for size of source cells
  Map                  targetBegin;                             //!< Define map for offset of target cells
  size_t               keysDevcSize;                            //!< Size of offsets for rangeDevc
  size_t               rangeDevcSize;                           //!< Size of offsets for sourceDevc
  size_t               sourceDevcSize;                          //!< Size of sources on device
  size_t               targetDevcSize;                          //!< Size of targets on device
  int                 *keysDevc;                                //!< Offsets for rangeDevc
  int                 *rangeDevc;                               //!< Offsets for sourceDevc
  gpureal             *sourceDevc;                              //!< Sources on device
  gpureal             *targetDevc;                              //!< Targets on device

  real *factorial;                                              //!< Factorial
  real *prefactor;                                              //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  real *Anm;                                                    //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  complex *Cnm;                                                 //!< M2L translation matrix \f$ C_{jn}^{km} \f$
public:
  real NP2P;                                                    //!< Number of P2P kernel calls
  real NM2P;                                                    //!< Number of M2P kernel calls
  real NM2L;                                                    //!< Number of M2L kernel calls

protected:
//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1;                                              // Initialize rho^m
    for( int m=0; m!=P; ++m ) {                                 // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real rhon = rhom;                                         //  rho^n
      for( int n=m+1; n!=P; ++n ) {                             //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const {
    const complex I(0.,1.);                                     // Imaginary unit
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
    for( int m=0; m!=2*P; ++m ) {                               // Loop over m in Ynm
      complex eim = std::exp(I * real(m * beta));               //  exp(i * m * beta)
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * prefactor[npn] * eim;               //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
      rhom /= rho;                                              //  rho^(-m-1)
      real rhon = rhom;                                         //  rho^(-n-1)
      for( int n=m+1; n!=2*P; ++n ) {                           //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eim;             //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
        rhon /= rho;                                            //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
    }                                                           // End loop over m in Ynm
  }

public:
//! Constructor
  KernelBase() : X0(0), R0(0), keysDevcSize(0), rangeDevcSize(0),
                 sourceDevcSize(0), targetDevcSize(0),
                 NP2P(0), NM2P(0), NM2L(0) {}
//! Destructor
  ~KernelBase() {}

//! Set center of root cell
  void setX0(vect x0) {X0 = x0;}
//! Set radius of root cell
  void setR0(real r0) {R0 = r0;}

//! Get center of root cell
  vect getX0() const {return X0;}
//! Get radius of root cell
  real getR0() const {return R0;}

//! Set center and size of root cell
  void setDomain(Bodies &bodies, vect x0=0, real r0=M_PI) {
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over bodies
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = (xmax[d] + xmin[d]) / 2;                          // Calculate center of domain
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( X0[0]-R0 < x0[0]-r0 || x0[0]+r0 < X0[0]+R0            //  Check for outliers in x direction
       || X0[1]-R0 < x0[1]-r0 || x0[1]+r0 < X0[1]+R0            //  Check for outliers in y direction
       || X0[2]-R0 < x0[2]-r0 || x0[2]+r0 < X0[2]+R0 ) {        //  Check for outliers in z direction
        std::cout << "Error: Particles located outside periodic domain : " << std::endl;// Print error message
        std::cout << X0-R0 << std::endl;
        std::cout << X0+R0 << std::endl;
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is r0
    } else {
      R0 += 1e-5;                                               // Add some leeway to root radius
    }                                                           // Endif for periodic boundary condition
  }

//! Precalculate M2L translation matrix
  void preCalculation() {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new real  [P];                                  // Factorial
    prefactor = new real  [4*P2];                               // sqrt( (n - |m|)! / (n + |m|)! )
    Anm       = new real  [4*P2];                               // (-1)^n / sqrt( (n + m)! / (n - m)! )
    Cnm       = new complex [P4];                               // M2L translation matrix Cjknm

    factorial[0] = 1;                                           // Initialize factorial
    for( int n=1; n!=P; ++n ) {                                 // Loop to P
      factorial[n] = factorial[n-1] * n;                        //  n!
    }                                                           // End loop to P

    for( int n=0; n!=2*P; ++n ) {                               // Loop over n in Anm
      for( int m=-n; m<=n; ++m ) {                              //  Loop over m in Anm
        int nm = n*n+n+m;                                       //   Index of Anm
        int nabsm = abs(m);                                     //   |m|
        real fnmm = EPS;                                        //   Initialize (n - m)!
        for( int i=1; i<=n-m; ++i ) fnmm *= i;                  //   (n - m)!
        real fnpm = EPS;                                        //   Initialize (n + m)!
        for( int i=1; i<=n+m; ++i ) fnpm *= i;                  //   (n + m)!
        real fnma = 1.0;                                        //   Initialize (n - |m|)!
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;              //   (n - |m|)!
        real fnpa = 1.0;                                        //   Initialize (n + |m|)!
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;              //   (n + |m|)!
        prefactor[nm] = std::sqrt(fnma/fnpa);                   //   sqrt( (n - |m|)! / (n + |m|)! )
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);              //   (-1)^n / sqrt( (n + m)! / (n - m)! )
      }                                                         //  End loop over m in Anm
    }                                                           // End loop over n in Anm

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {                   // Loop over j in Cjknm
      for( int k=-j; k<=j; ++k, ++jk ){                         //  Loop over k in Cjknm
        for( int n=0, nm=0; n!=P; ++n ) {                       //   Loop over n in Cjknm
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {            //    Loop over m in Cjknm
            const int jnkm = (j+n)*(j+n)+j+n+m-k;               //     Index C_{j+n}^{m-k}
            Cnm[jknm] = std::pow(I,real(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                      * real(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * EPS;
          }                                                     //    End loop over m in Cjknm
        }                                                       //   End loop over n in Cjknm
      }                                                         //  End loop over in k in Cjknm
    }                                                           // End loop over in j in Cjknm
  }

//! Free temporary allocations
  void postCalculation() {
    delete[] factorial;                                         // Free factorial
    delete[] prefactor;                                         // Free sqrt( (n - |m|)! / (n + |m|)! )
    delete[] Anm;                                               // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
    delete[] Cnm;                                               // Free M2L translation matrix Cjknm
  }

//! Set paramters for Van der Waals
  void setVanDerWaals(int atoms, double *rscale, double *gscale) {
    assert(atoms <= 16);                                        // Change GPU constant memory alloc if needed
    THETA = .1;                                                 // Force opening angle to be small
    ATOMS = atoms;                                              // Set number of atom types
    RSCALE.resize(ATOMS*ATOMS);                                 // Resize rscale vector
    GSCALE.resize(ATOMS*ATOMS);                                 // Resize gscale vector
    for( int i=0; i!=ATOMS*ATOMS; ++i ) {                       // Loop over scale vector
      RSCALE[i] = rscale[i];                                    //  Set rscale vector
      GSCALE[i] = gscale[i];                                    //  Set gscale vector
    }                                                           // End loop over scale vector
  }

//! Set paramters for Ewald summation
  void setEwald(real ksize, real alpha, real sigma) {
    KSIZE = ksize;                                              // Set number of waves
    ALPHA = alpha;                                              // Set scaling parameter
    SIGMA = sigma;                                              // Set scaling parameter
  }

};

template<Equation equation>
class Kernel : public KernelBase {
public:
  void initialize();                                            //!< Initialize kernels
  void P2M(C_iter Ci) const;                                    //!< Evaluate P2M kernel on CPU
  void M2M(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2M kernel on CPU
  void M2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2L kernel on CPU
  void M2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2P kernel on CPU
  void P2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate P2P kernel on CPU
  void L2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate L2L kernel on CPU
  void L2P(C_iter Ci) const;                                    //!< Evaluate L2P kernel on CPU
  void EwaldReal(C_iter Ci, C_iter Cj) const;                   //!< Evaluate Ewald real part on CPU
  void EwaldWave(Bodies &bodies) const;                         //!< Evaluate Ewald wave part on CPU
  void P2M();                                                   //!< Evaluate P2M kernel on GPU
  void M2M();                                                   //!< Evaluate M2M kernel on GPU
  void M2L();                                                   //!< Evaluate M2L kernel on GPU
  void M2P();                                                   //!< Evaluate M2P kernel on GPU
  void P2P();                                                   //!< Evalaute P2P kernel on GPU
  void L2L();                                                   //!< Evaluate L2L kernel on GPU
  void L2P();                                                   //!< Evaluate L2P kernel on GPU
  void EwaldReal();                                             //!< Evaluate Ewald real part on GPU
  void EwaldWave();                                             //!< Evalaute Ewald wave part on GPU
  void finalize();                                              //!< Finalize kernels

  void allocate();                                              //!< Allocate GPU variables
  void hostToDevice();                                          //!< Copy from host to device
  void deviceToHost();
};

#endif
