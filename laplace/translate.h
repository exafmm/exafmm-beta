#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)
#define IPOW2N(n) ((n >= 0) ? 1 : ODDEVEN(n))
const complex_t I(0.0,1.0);

void cart2sph(real_t & r, real_t & theta, real_t & phi, vec3 dX) {
  r = sqrt(norm(dX));
  theta = r == 0 ? 0 : acos(dX[2] / r);
  phi = atan2(dX[1], dX[0]);
}

template<typename T>
void sph2cart(real_t r, real_t theta, real_t phi, T spherical, T & cartesian) {
  cartesian[0] = sin(theta) * cos(phi) * spherical[0]
    + cos(theta) * cos(phi) / r * spherical[1]
    - sin(phi) / r / sin(theta) * spherical[2];
  cartesian[1] = sin(theta) * sin(phi) * spherical[0]
    + cos(theta) * sin(phi) / r * spherical[1]
    + cos(phi) / r / sin(theta) * spherical[2];
  cartesian[2] = cos(theta) * spherical[0]
    - sin(theta) / r * spherical[1];
}

void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
  real_t x = std::cos(alpha);
  real_t y = std::sin(alpha);
  real_t fact = 1;
  real_t pn = 1;
  real_t rhom = 1;
  complex_t ei = std::exp(I * beta);
  complex_t eim = 1.0;
  for (int m=0; m<P; m++) {
    real_t p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    Ynm[npn] = rhom * p * eim;
    Ynm[nmn] = std::conj(Ynm[npn]);
    real_t p1 = p;
    p = x * (2 * m + 1) * p1;
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * eim;
    rhom *= rho;
    real_t rhon = rhom;
    for (int n=m+1; n<P; n++) {
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      rhon /= -(n + m);
      Ynm[npm] = rhon * p * eim;
      Ynm[nmm] = std::conj(Ynm[npm]);
      real_t p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * eim;
      rhon *= rho;
    }
    rhom /= -(2 * m + 2) * (2 * m + 1);
    pn = -pn * fact * y;
    fact += 2;
    eim *= ei;
  }
}

void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
  real_t x = std::cos(alpha);
  real_t y = std::sin(alpha);
  real_t fact = 1;
  real_t pn = 1;
  real_t invR = -1.0 / rho;
  real_t rhom = -invR;
  complex_t ei = std::exp(I * beta);
  complex_t eim = 1.0;
  for (int m=0; m<P; m++) {
    real_t p = pn;
    int npn = m * m + 2 * m;
    int nmn = m * m;
    Ynm[npn] = rhom * p * eim;
    Ynm[nmn] = std::conj(Ynm[npn]);
    real_t p1 = p;
    p = x * (2 * m + 1) * p1;
    rhom *= invR;
    real_t rhon = rhom;
    for (int n=m+1; n<P; n++) {
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      Ynm[npm] = rhon * p * eim;
      Ynm[nmm] = std::conj(Ynm[npm]);
      real_t p2 = p1;
      p1 = p;
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
      rhon *= invR * (n - m + 1);
    }
    pn = -pn * fact * y;
    fact += 2;
    eim *= ei;
  }
}
