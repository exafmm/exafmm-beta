#include "kernel.h"
#include "spherical.h"

void Kernel::initialize() {
  precalculate();
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

void Kernel::M2L() {
  vect dist = CI->X - CJ->X - Xperiodic;
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
    vect dist = B->pos - CJ->X - Xperiodic;
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
      vect dist = BI->pos - BJ->pos - Xperiodic;
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
  postcalculate();
}
