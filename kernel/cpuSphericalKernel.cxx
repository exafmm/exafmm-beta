#include "kernel.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const int  P2 = P * P;
const int  P4 = P2 * P2;
const real EPS = 1e-6;
double *prefactor, *Anm;
complex *Ynm, *YnmTheta, *Cnm, I(0.0,1.0);

void Kernel::initialize() {
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

void evalMultipole(real rho, real alpha, real beta) {
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

void Kernel::P2M() {
  for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {
    vect dist = B->pos - CJ->X;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        const int nm  = n * n + n + m;
        const int nms = n * (n + 1) / 2 + m;
        CJ->M[nms] += double(B->scal)*Ynm[nm];
      }
    }
  }
}

void Kernel::M2M() {
  vect dist = CI->X - CJ->X;
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
            M += CJ->M[jnkms]*std::pow(I,double(m-abs(m)))*Ynm[nm]*double(ODDEVEN(n)*Anm[nm]*Anm[jnkm]/Anm[jk]);
          }
        }
        for( int m=k; m<=n; ++m ) {
          if( j-n >= m-k ) {
            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            const int nm    = n * n + n + m;
            M += std::conj(CJ->M[jnkms])*Ynm[nm]*double(ODDEVEN(k+n+m)*Anm[nm]*Anm[jnkm]/Anm[jk]);
          }
        }
      }
      CI->M[jks] += M;
    }
  }
}

void Kernel::M2L() {
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
          L += std::conj(CJ->M[nms])*Cnm[jknm]*Ynm[jnkm];
        }
        for( int m=0; m<=n; ++m ) {
          const int nm   = n * n + n + m;
          const int nms  = n * (n + 1) / 2 + m;
          const int jknm = jk * P2 + nm;
          const int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += CJ->M[nms]*Cnm[jknm]*Ynm[jnkm];
        }
      }
      CI->L[jks] += L;
    }
  }
}

void Kernel::M2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CJ->X;
    vect acc = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->pot += (CJ->M[nms]*Ynm[nm]).real();
      acc[0] -= (CJ->M[nms]*Ynm[nm]).real()/r*(n+1);
      acc[1] += (CJ->M[nms]*YnmTheta[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->pot += 2*(CJ->M[nms]*Ynm[nm]).real();
        acc[0] -= 2*(CJ->M[nms]*Ynm[nm]).real()/r*(n+1);
        acc[1] += 2*(CJ->M[nms]*YnmTheta[nm]).real();
        acc[2] += 2*(CJ->M[nms]*Ynm[nm]*I).real()*m;
      }
    }
    B->acc[0] += sin(theta)*cos(phi)*acc[0]+cos(theta)*cos(phi)/r*acc[1]-sin(phi)/r/sin(theta)*acc[2];
    B->acc[1] += sin(theta)*sin(phi)*acc[0]+cos(theta)*sin(phi)/r*acc[1]+cos(phi)/r/sin(theta)*acc[2];
    B->acc[2] += cos(theta)*acc[0]-sin(theta)/r*acc[1];
  }
}

void Kernel::P2P() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->pos - BJ->pos;
      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR3 = BJ->scal * invR * invR * invR;
      BI->pot += BJ->scal * invR;
      BI->acc -= dist * invR3;
    }
  }
}

void Kernel::L2L() {
  vect dist = CI->X - CJ->X;
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
          L += std::conj(CJ->L[nms])*Ynm[jnkm]*double(ODDEVEN(k)*Anm[jnkm]*Anm[jk]/Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            const int jnkm = (n - j) * (n - j) + n - j + m - k;
            const int nm   = n * n + n + m;
            const int nms  = n * (n + 1) / 2 + m;
            L += CJ->L[nms]*std::pow(I,double(m-k-abs(m-k)))*Ynm[jnkm]*double(Anm[jnkm]*Anm[jk]/Anm[nm]);
          }
        }
      }
      CI->L[jks] += L;
    }
  }
}

void Kernel::L2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->pos - CI->X;
    vect acc = 0;
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi);
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->pot += (CI->L[nms]*Ynm[nm]).real();
      acc[0] += (CI->L[nms]*Ynm[nm]).real()/r*n;
      acc[1] += (CI->L[nms]*YnmTheta[nm]).real();
      for( int m=1; m<=n; ++m ) {
        nm  = n * n + n + m;
        nms = n * (n + 1) / 2 + m;
        B->pot += 2*(CI->L[nms]*Ynm[nm]).real();
        acc[0] += 2*(CI->L[nms]*Ynm[nm]).real()/r*n;
        acc[1] += 2*(CI->L[nms]*YnmTheta[nm]).real();
        acc[2] += 2*(CI->L[nms]*Ynm[nm]*I).real()*m;
      }
    }
    B->acc[0] += sin(theta)*cos(phi)*acc[0]+cos(theta)*cos(phi)/r*acc[1]-sin(phi)/r/sin(theta)*acc[2];
    B->acc[1] += sin(theta)*sin(phi)*acc[0]+cos(theta)*sin(phi)/r*acc[1]+cos(phi)/r/sin(theta)*acc[2];
    B->acc[2] += cos(theta)*acc[0]-sin(theta)/r*acc[1];
  }
}

void Kernel::finalize() {
  delete[] prefactor;
  delete[] Anm;
  delete[] Ynm;
  delete[] YnmTheta;
  delete[] Cnm;
}
