#ifndef spherical_h
#define spherical_h

const int  P2 = P * P;
const int  P4 = P2 * P2;
const real EPS = 1e-6;
static double *prefactor, *Anm;
static complex *Ynm, *YnmTheta, *Cnm, I(0.0,1.0);

namespace {
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
}

#endif
