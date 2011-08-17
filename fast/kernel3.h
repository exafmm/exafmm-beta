#ifndef kernel_h
#define kernel_h
#include <types.h>
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;
const int  P4 = P2 * P2;
const real EPS = 1e-6; 

class Kernel {
private:
  real DMAX;
  double *factorial, *prefactor, *Anm;
  complex *Ynm, *YnmTheta, *Cnm;

protected:
  vect   X0;
  real   R0;
  B_iter B0, BN;
  C_iter C0, CN;

private:
  real getBmax(vect const&X, C_iter C) {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  void cart2sph(real& r, real& theta, real& phi, vect dist) const {
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
  void sph2cart(real r, real theta, real phi, T spherical, T &cartesian) const {
    cartesian[0] = sin(theta) * cos(phi) * spherical[0]
                 + cos(theta) * cos(phi) / r * spherical[1]
                 - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0]
                 + cos(theta) * sin(phi) / r * spherical[1]
                 + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]
                 - sin(theta) / r * spherical[1];
  }

  void evalMultipole(real rho, real alpha, real beta) const {
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

  void evalLocal(real rho, real alpha, real beta) const {
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

public:
  Kernel() : X0(0), R0(0) {
    const complex I(0.,1.);                                     // Imaginary unit
    factorial = new double  [P];
    prefactor = new double  [4*P2];
    Anm       = new double  [4*P2];
    Ynm       = new complex [4*P2];
    YnmTheta  = new complex [4*P2];
    Cnm       = new complex [P4];

    factorial[0] = 1;
    for( int n=1; n!=P; ++n ) {
      factorial[n] = factorial[n-1] * n;
    }

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

  ~Kernel() {
    delete[] factorial;
    delete[] prefactor;
    delete[] Anm;
    delete[] Ynm;
    delete[] YnmTheta;
    delete[] Cnm;
  }

  void setCenter(C_iter C) {
    DMAX = 0;
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += B->SRC[0];
      X += B->X * B->SRC[0];
    }
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      m += c->M[0].real();
      X += c->X * c->M[0].real();
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  void P2M(C_iter C) {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vect dist = C->X - B->X;
      real R = std::sqrt(norm(dist));
      if( R > DMAX ) DMAX = R;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta);
      for( int n=0; n!=P; ++n ) {
        for( int m=0; m<=n; ++m ) {
          const int nm  = n * n + n + m;
          const int nms = n * (n + 1) / 2 + m;
          C->M[nms] += double(B->SRC[0]) * Ynm[nm];
        }
      }
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2M(C_iter CI) {
    for( C_iter CJ=C0+CI->CHILD; CJ!=C0+CI->CHILD+CI->NCHILD; ++CJ ) {
      vect dist = CI->X - CJ->X;
      real R = std::sqrt(norm(dist)) + CJ->RCRIT;
      if( R > DMAX ) DMAX = R;
      const complex I(0.,1.);                                       // Imaginary unit
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta);
      for( int j=0; j!=P; ++j ) {
        for( int k=0; k<=j; ++k ) {
          const int jk = j * j + j + k;
          const int jks = j * (j + 1) / 2 + k;
          complex M = 0;
          for( int n=0; n<=j; ++n ) {
            for( int m=-n; m<=std::min(k-1,n); ++m ) {
              if( j-n >= k-m ) {
                const int jnkm  = (j - n) * (j - n) + j - n + k - m;
                const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                const int nm    = n * n + n + m;
                M += CJ->M[jnkms] * std::pow(I,double(m-abs(m))) * Ynm[nm]
                   * double(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
              }
            }
            for( int m=k; m<=n; ++m ) {
              if( j-n >= m-k ) {
                const int jnkm  = (j - n) * (j - n) + j - n + k - m;
                const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                const int nm    = n * n + n + m;
                M += std::conj(CJ->M[jnkms]) * Ynm[nm]
                   * double(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
              }
            }
          }
          CI->M[jks] += M;
        }
      }
    }
    CI->RCRIT = std::min(CI->R,DMAX);
  }

  void P2P(C_iter CI, C_iter CJ, bool mutual=true) const {
    for( B_iter BI=CI->LEAF; BI!=CI->LEAF+CI->NDLEAF; ++BI ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=CJ->LEAF; BJ!=CJ->LEAF+CJ->NDLEAF; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0 * mutual;
        BJ->TRG[1] += dR[0] * mutual;
        BJ->TRG[2] += dR[1] * mutual;
        BJ->TRG[3] += dR[2] * mutual;
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void P2P(C_iter C) const {
    unsigned NJ = C->NDLEAF;
    for( B_iter BI=C->LEAF; BI!=C->LEAF+C->NDLEAF; ++BI, --NJ ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=BI+1; BJ!=BI+NJ; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0;
        BJ->TRG[1] += dR[0];
        BJ->TRG[2] += dR[1];
        BJ->TRG[3] += dR[2];
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void M2L(C_iter CI, C_iter CJ, bool mutual=true) const {
    vect dist = CI->X - CJ->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalLocal(rho,alpha,beta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j * j + j + k;
        const int jks = j * (j + 1) / 2 + k;
        complex L = 0;
        for( int n=0; n!=P; ++n ) {
          for( int m=-n; m<0; ++m ) {
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 - m;
            const int jknm = jk * P2 + nm;
            const int jnkm = (j + n) * (j + n) + j + n + m - k;
            L += std::conj(CJ->M[nms]) * Cnm[jknm] * Ynm[jnkm];
          }
          for( int m=0; m<=n; ++m ) {
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
            const int jknm = jk * P2 + nm;
            const int jnkm = (j + n) * (j + n) + j + n + m - k;
            L += CJ->M[nms] * Cnm[jknm] * Ynm[jnkm];
          }
        }
        CI->L[jks] += L;
      }
    }
    rho = mutual;
  }

  void L2L(C_iter CI) const {
    C_iter CJ = C0 + CI->PARENT;
    vect dist = CI->X - CJ->X;
    const complex I(0.,1.);                                       // Imaginary unit
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,beta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j * j + j + k;
        const int jks = j * (j + 1) / 2 + k;
        complex L = 0;
        for( int n=j; n!=P; ++n ) {
          for( int m=j+k-n; m<0; ++m ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n - m;
            const int nms  = n * (n + 1) / 2 - m;
            L += std::conj(CJ->L[nms]) * Ynm[jnkm]
               * double(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
          for( int m=0; m<=n; ++m ) {
            if( n-j >= abs(m-k) ) {
              const int jnkm = (n - j) * (n - j) + n - j + m - k;
              const int nm   = n * n + n + m;
              const int nms  = n * (n + 1) / 2 + m;
              L += CJ->L[nms] * std::pow(I,double(m-k-abs(m-k)))
                 * Ynm[jnkm] * double(Anm[jnkm] * Anm[jk] / Anm[nm]);
            }
          }
        }
        CI->L[jks] += L;
      }
    }
  }

  void L2P(C_iter CI) const {
    const complex I(0.,1.);                                       // Imaginary unit
    for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NCLEAF; ++B ) {
      vect dist = B->X - CI->X;
      B->TRG /= B->SRC[0];
      vect spherical = 0;
      vect cartesian = 0;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipole(r,theta,phi);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        B->TRG[0] -= (CI->L[nms] * Ynm[nm]).real();
        spherical[0] += (CI->L[nms] * Ynm[nm]).real() / r * n;
        spherical[1] += (CI->L[nms] * YnmTheta[nm]).real();
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          B->TRG[0] -= 2 * (CI->L[nms] * Ynm[nm]).real();
          spherical[0] += 2 * (CI->L[nms] * Ynm[nm]).real() / r * n;
          spherical[1] += 2 * (CI->L[nms] * YnmTheta[nm]).real();
          spherical[2] += 2 * (CI->L[nms] * Ynm[nm] * I).real() * m;
        }
      }
      sph2cart(r,theta,phi,spherical,cartesian);
      B->TRG[1] += cartesian[0];
      B->TRG[2] += cartesian[1];
      B->TRG[3] += cartesian[2];
    }
  }
};

#endif
