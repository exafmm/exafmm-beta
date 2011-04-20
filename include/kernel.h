#ifndef kernel_h
#define kernel_h
#define KERNEL
#include "logger.h"
#undef KERNEL
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;                                          // P^2
const int  P4 = P2 * P2;                                        // P^4
const real EPS = 1e-6;                                          // Single precision epsilon

class Kernel : public Logger {                                  // Unified CPU/GPU kernel class
protected:
  B_iter BI0;                                                   // Target bodies begin iterator
  B_iter BIN;                                                   // Target bodies end iterator
  B_iter BJ0;                                                   // Source bodies begin iterator
  B_iter BJN;                                                   // Source bodies end iterator
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator
  vect   X0;                                                    // Center of root cell
  real   R0;                                                    // Radius of root cell
  vect   Xperiodic;                                             // Coordinate offset of periodic image

  std::vector<int>     keysHost;                                // Offsets for rangeHost
  std::vector<int>     rangeHost;                               // Offsets for sourceHost
  std::vector<gpureal> constHost;                               // Constants on host
  std::vector<gpureal> sourceHost;                              // Sources on host
  std::vector<gpureal> targetHost;                              // Targets on host
  Map                  sourceBegin;                             // Define map for offset of source cells
  Map                  sourceSize;                              // Define map for size of source cells
  Map                  targetBegin;                             // Define map for offset of target cells

  double *prefactor, *Anm;                                      // Auxiliary variables for spherical harmonics
  complex *Ynm, *YnmTheta, *Cnm;                                // Auxiliary variables for spherical harmonics

private:
  void cart2sph(real& r, real& theta, real& phi, vect dist) {
    r = std::sqrt(norm(dist))+EPS;
    theta = std::acos(dist[2] / r);
    if( std::abs(dist[0]) + std::abs(dist[1]) < EPS ) {
      phi = 0;
    } else if( std::abs(dist[0]) < EPS ) {
      phi = dist[1] / std::abs(dist[1]) * M_PI * 0.5;
    } else if( dist[0] > 0 ) {
      phi = std::atan(dist[1] / dist[0]);
    } else {
      phi = std::atan(dist[1] / dist[0]) + M_PI;
    }
  }

  template<typename T>
  void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) {
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]
                 + cos(theta) * cos(phi) / r * spherical[1]
                 - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]
                 + cos(theta) * sin(phi) / r * spherical[1]
                 + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]
                 - sin(theta) / r * spherical[1];
  }

  void evalMultipole(real rho, real alpha, real beta) {
    const complex I(0.,1.);                                     // Imaginary unit
    double x = std::cos(alpha);
    double y = std::sin(alpha);
    double s = std::sqrt(1 - x * x);
    double fact = 1;
    double pn = 1;
    double rhom = 1;
    for( int m=0; m!=P; ++m ) {
      complex eim = std::exp(I * double(m * beta));
      double p = pn;
      int npn = m * m + 2 * m;
      int nmn = m * m;
      Ynm[npn] = rhom * p * prefactor[npn] * eim;
      Ynm[nmn] = std::conj(Ynm[npn]);
      double p1 = p;
      p = x * (2 * m + 1) * p;
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;
      rhom *= rho;
      double rhon = rhom;
      for( int n=m+1; n!=P; ++n ) {
        int npm = n * n + n + m;
        int nmm = n * n + n - m;
        Ynm[npm] = rhon * p * prefactor[npm] * eim;
        Ynm[nmm] = std::conj(Ynm[npm]);
        double p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;
        rhon *= rho;
      }
      pn = -pn * fact * s;
      fact += 2;
    }
  }

  void evalLocal(real rho, real alpha, real beta) {
    const complex I(0.,1.);                                     // Imaginary unit
    double x = std::cos(alpha);
    double y = std::sin(alpha);
    double s = std::sqrt(1 - x * x);
    double fact = 1;
    double pn = 1;
    double rhom = 1.0 / rho;
    for( int m=0; m!=2*P; ++m ) {
      complex eim = std::exp(I * double(m * beta));
      double p = pn;
      int npn = m * m + 2 * m;
      int nmn = m * m;
      Ynm[npn] = rhom * p * prefactor[npn] * eim;
      Ynm[nmn] = std::conj(Ynm[npn]);
      double p1 = p;
      p = x * (2 * m + 1) * p;
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;
      rhom /= rho;
      double rhon = rhom;
      for( int n=m+1; n!=2*P; ++n ) {
        int npm = n * n + n + m;
        int nmm = n * n + n - m;
        Ynm[npm] = rhon * p * prefactor[npm] * eim;
        Ynm[nmm] = std::conj(Ynm[npm]);
        double p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;
        rhon /= rho;
      }
      pn = -pn * fact * s;
      fact += 2;
    }
  }

protected:
  int getLevel(bigint index) {                                  // Get level from cell index
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

public:
  Kernel() : X0(0), R0(0) {}                                    // Constructor
  ~Kernel() {}                                                  // Destructor

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain(Bodies &bodies) {                              // Set center and size of root cell
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    X0 = 0;                                                     // Initialize center and size of root cell
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->X;                                               //  Sum positions
    }                                                           // End loop over bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 += 1e-5;                                                 // Add some leeway to root radius
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      X0 = 0;                                                   //  Center is [0, 0, 0]
      R0 = M_PI;                                                //  Radius is M_PI
    }                                                           // Endif for periodic boundary condition
  }

  void preCalculation() {
    const complex I(0.,1.);                                     // Imaginary unit
    prefactor = new double  [4*P2];
    Anm       = new double  [4*P2];
    Ynm       = new complex [4*P2];
    YnmTheta  = new complex [4*P2];
    Cnm       = new complex [P4];

    for( int n=0; n!=2*P; ++n ) {
      for( int m=-n; m<=n; ++m ) {
        int nm = n*n+n+m;
        int nabsm = abs(m);
        double fnmm = 1.0;
        for( int i=1; i<=n-m; ++i ) fnmm *= i;
        double fnpm = 1.0;
        for( int i=1; i<=n+m; ++i ) fnpm *= i;
        double fnma = 1.0;
        for( int i=1; i<=n-nabsm; ++i ) fnma *= i;
        double fnpa = 1.0;
        for( int i=1; i<=n+nabsm; ++i ) fnpa *= i;
        prefactor[nm] = std::sqrt(fnma/fnpa);
        Anm[nm] = ODDEVEN(n)/std::sqrt(fnmm*fnpm);
      }
    }

    for( int j=0, jk=0, jknm=0; j!=P; ++j ) {
      for( int k=-j; k<=j; ++k, ++jk ){
        for( int n=0, nm=0; n!=P; ++n ) {
          for( int m=-n; m<=n; ++m, ++nm, ++jknm ) {
            const int jnkm = (j+n)*(j+n)+j+n+m-k;
            Cnm[jknm] = std::pow(I,double(abs(k-m)-abs(k)-abs(m)))*(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]);
          }
        }
      }
    }
  }

  void postCalculation() {
    delete[] prefactor;
    delete[] Anm;
    delete[] Ynm;
    delete[] YnmTheta;
    delete[] Cnm;
  }

  void LaplaceInit();
  void LaplaceP2M();
  void LaplaceM2M();
  void LaplaceM2M_CPU();
  void LaplaceM2L();
  void LaplaceM2P();
  void LaplaceP2P();
  void LaplaceP2P_CPU();
  void LaplaceL2L();
  void LaplaceL2P();
  void LaplaceFinal();

  void BiotSavartInit();
  void BiotSavartP2M();
  void BiotSavartM2M();
  void BiotSavartM2M_CPU();
  void BiotSavartM2L();
  void BiotSavartM2P();
  void BiotSavartP2P();
  void BiotSavartP2P_CPU();
  void BiotSavartL2L();
  void BiotSavartL2P();
  void BiotSavartFinal();

  void StretchingInit();
  void StretchingP2M();
  void StretchingM2M();
  void StretchingM2M_CPU();
  void StretchingM2L();
  void StretchingM2P();
  void StretchingP2P();
  void StretchingP2P_CPU();
  void StretchingL2L();
  void StretchingL2P();
  void StretchingFinal();

  void GaussianInit();
  void GaussianP2M();
  void GaussianM2M();
  void GaussianM2M_CPU();
  void GaussianM2L();
  void GaussianM2P();
  void GaussianP2P();
  void GaussianP2P_CPU();
  void GaussianL2L();
  void GaussianL2P();
  void GaussianFinal();
};

#endif
