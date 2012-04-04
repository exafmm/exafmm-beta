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

//! Unified CPU/GPU kernel class
class KernelBase : public Sort {
protected:
  C_iter               Ci0;                                     //!< Begin iterator for target cells
  C_iter               Cj0;                                     //!< Begin iterator for source cells

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
  vect X0;                                                      //!< Center of root cell
  real R0;                                                      //!< Radius of root cell

private:
  real getBmax(vect const&X, C_iter C) const {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

protected:
  void setCenter(C_iter C) const {
/*
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += std::abs(B->SRC);
      X += B->X * std::abs(B->SRC);
    }
    for( C_iter c=Cj0+C->CHILD; c!=Cj0+C->CHILD+C->NCHILD; ++c ) {
      m += std::abs(c->M[0]);
      X += c->X * std::abs(c->M[0]);
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
*/
  }

public:
//! Constructor
  KernelBase() : Ci0(), Cj0(), ATOMS(), RSCALE(), GSCALE(), KSIZE(), ALPHA(), SIGMA(),
                 keysHost(), rangeHost(), sourceHost(), targetHost(), constHost(),
                 sourceBegin(), sourceSize(), targetBegin(),
                 keysDevcSize(0), rangeDevcSize(0),
                 sourceDevcSize(0), targetDevcSize(0),
                 keysDevc(), rangeDevc(), sourceDevc(), targetDevc(),
                 factorial(), prefactor(), Anm(), Cnm(),
                 X0(0), R0(0) {}
//! Destructor
  ~KernelBase() {}
//! Copy constructor
  KernelBase(const KernelBase&) : Sort(), Ci0(), Cj0(), ATOMS(), RSCALE(), GSCALE(), KSIZE(), ALPHA(), SIGMA(),
                 keysHost(), rangeHost(), sourceHost(), targetHost(), constHost(),
                 sourceBegin(), sourceSize(), targetBegin(),
                 keysDevcSize(0), rangeDevcSize(0),
                 sourceDevcSize(0), targetDevcSize(0),
                 keysDevc(), rangeDevc(), sourceDevc(), targetDevc(),
                 factorial(), prefactor(), Anm(), Cnm(),
                 X0(0), R0(0) {}
//! Overload assignment
  KernelBase &operator=(const KernelBase) {return *this;}

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
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( xmin[0] < x0[0]-r0 || x0[0]+r0 < xmax[0]              //  Check for outliers in x direction
       || xmin[1] < x0[1]-r0 || x0[1]+r0 < xmax[1]              //  Check for outliers in y direction
       || xmin[2] < x0[2]-r0 || x0[2]+r0 < xmax[2] ) {          //  Check for outliers in z direction
        std::cout << "Error: Particles located outside periodic domain : " << std::endl;// Print error message
        std::cout << xmin << std::endl;
        std::cout << xmax << std::endl;
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is r0
    } else {
      for( int d=0; d!=3; ++d ) {                               // Loop over each dimension
        X0[d] = (xmax[d] + xmin[d]) / 2;                        // Calculate center of domain
        X0[d] = int(X0[d]+.5);                                  //  Shift center to nearest integer
        R0 = std::max(xmax[d] - X0[d], R0);                     //  Calculate max distance from center
        R0 = std::max(X0[d] - xmin[d], R0);                     //  Calculate max distance from center
      }                                                         // End loop over each dimension
      R0 *= 1.000001;                                           // Add some leeway to root radius
    }                                                           // Endif for periodic boundary condition
  }

//! Precalculate M2L translation matrix
  void preCalculation() {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new real  [P];                                  // Factorial
    prefactor = new real  [4*P*P];                              // sqrt( (n - |m|)! / (n + |m|)! )
    Anm       = new real  [4*P*P];                              // (-1)^n / sqrt( (n + m)! / (n - m)! )
    Cnm       = new complex [P*P*P*P];                          // M2L translation matrix Cjknm

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
  void evalMultipole(real rho, real alpha, real beta, complex *Ynm) const;//!< Set Spherical multipole expansion
  void evalMultipoleTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const;//!< With dY/dtheta
  void evalLocal(real rho, real alpha, real beta, complex *Ynm) const;//!< Set Spherical local expansion
  void evalLocalTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const;//!< With dY/dtheta
  void P2M(C_iter Ci) const;                                    //!< Evaluate P2M kernel on CPU
  void M2M(C_iter Ci) const;                                    //!< Evaluate M2M kernel on CPU
  void M2L(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2L kernel on CPU
  void M2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate M2P kernel on CPU
  void P2P(C_iter Ci, C_iter Cj) const;                         //!< Evaluate P2P kernel on CPU
  void L2L(C_iter Ci) const;                                    //!< Evaluate L2L kernel on CPU
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
