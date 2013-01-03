#include "kernel.h"

#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

const real_t EPS = 1e-6;
const real_t SCALING = 1e-6;

//! Get r,theta,phi from x,y,z
void cart2sph(real_t& r, real_t& theta, real_t& phi, vec3 dX) {
  r = sqrt(norm(dX)) * (1 + EPS);                               // r = sqrt(x^2 + y^2 + z^2)
  if( r < EPS ) {                                               // If r == 0
    theta = 0;                                                  //  theta can be anything so we set it to 0
  } else {                                                      // If r != 0
    theta = acos(dX[2] / r);                                    //  theta = acos(z / r)
  }                                                             // End if for r == 0
  if( fabs(dX[0]) + fabs(dX[1]) < EPS ) {                       // If |x| < eps & |y| < eps
    phi = 0;                                                    //  phi can be anything so we set it to 0
  } else if( fabs(dX[0]) < EPS ) {                              // If |x| < eps
    phi = dX[1] / fabs(dX[1]) * M_PI * 0.5;                     //  phi = sign(y) * pi / 2
  } else if( dX[0] > 0 ) {                                      // If x > 0
    phi = atan(dX[1] / dX[0]);                                  //  phi = atan(y / x)
  } else {                                                      // If x < 0
    phi = atan(dX[1] / dX[0]) + M_PI;                           //  phi = atan(y / x) + pi
  }                                                             // End if for x,y cases
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
void evalMultipole(real_t rho, real_t alpha, real_t beta, real_t *prefactor, complex_t *Ynm, complex_t *YnmTheta) {
  const complex_t I(0.,1.);                                     // Imaginary unit
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = 1;                                              // Initialize rho^m
  for( int m=0; m!=P; ++m ) {                                   // Loop over m in Ynm
    complex_t eim = std::exp(I * real_t(m * beta));             //  exp(i * m * beta)
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * prefactor[npn] * eim;                 //  rho^m * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
    rhom *= rho;                                                //  rho^m
    real_t rhon = rhom;                                         //  rho^n
    for( int n=m+1; n!=P; ++n ) {                               //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * prefactor[npm] * eim;               //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
      rhon *= rho;                                              //   Update rho^n
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
  }                                                             // End loop over m in Ynm
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
void evalLocal(real_t rho, real_t alpha, real_t beta, real_t *prefactor, complex_t *Ynm, complex_t *YnmTheta) {
  const complex_t I(0.,1.);                                     // Imaginary unit
  real_t x = std::cos(alpha);                                   // x = cos(alpha)
  real_t y = std::sin(alpha);                                   // y = sin(alpha)
  real_t fact = 1;                                              // Initialize 2 * m + 1
  real_t pn = 1;                                                // Initialize Legendre polynomial Pn
  real_t rhom = 1.0 / rho;                                      // Initialize rho^(-m-1)
  for( int m=0; m!=P; ++m ) {                                   // Loop over m in Ynm
    complex_t eim = std::exp(I * real_t(m * beta));             //  exp(i * m * beta)
    real_t p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * prefactor[npn] * eim;                 //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_t p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * prefactor[npn] * eim;// theta derivative of r^n * Ynm
    rhom /= rho;                                                //  rho^(-m-1)
    real_t rhon = rhom;                                         //  rho^(-n-1)
    for( int n=m+1; n!=P; ++n ) {                               //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * prefactor[npm] * eim;               //   rho^n * Ynm for m > 0
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_t p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * prefactor[npm] * eim;// theta derivative
      rhon /= rho;                                              //   rho^(-n-1)
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
  }                                                             // End loop over m in Ynm
}

#if __SSE__
inline float vecSum4(__m128 reg) {
  float mem[4];
  _mm_store_ps(mem, reg);
  return mem[0] + mem[1] + mem[2] + mem[3];
}
#endif

#if __AVX__
inline float vecSum8(__m256 reg) {
  float mem[8];
  _mm256_store_ps(mem, reg);
  return mem[0] + mem[1] + mem[2] + mem[3] + mem[4] + mem[5] + mem[6] + mem[7];
}
#endif

void Kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) const {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  int i = 0;
#if __AVX__
  for ( ; i<=ni-8; i+=8) {
    __m256 pot = _mm256_setzero_ps();
    __m256 ax = _mm256_setzero_ps();
    __m256 ay = _mm256_setzero_ps();
    __m256 az = _mm256_setzero_ps();

    __m256 xi = _mm256_setr_ps(Bi[i].X[0],Bi[i+1].X[0],Bi[i+2].X[0],Bi[i+3].X[0],
      Bi[i+4].X[0],Bi[i+5].X[0],Bi[i+6].X[0],Bi[i+7].X[0]) - _mm256_set1_ps(Xperiodic[0]);
    __m256 yi = _mm256_setr_ps(Bi[i].X[1],Bi[i+1].X[1],Bi[i+2].X[1],Bi[i+3].X[1],
      Bi[i+4].X[1],Bi[i+5].X[1],Bi[i+6].X[1],Bi[i+7].X[1]) - _mm256_set1_ps(Xperiodic[1]);
    __m256 zi = _mm256_setr_ps(Bi[i].X[2],Bi[i+1].X[2],Bi[i+2].X[2],Bi[i+3].X[2],
      Bi[i+4].X[2],Bi[i+5].X[2],Bi[i+6].X[2],Bi[i+7].X[2]) - _mm256_set1_ps(Xperiodic[2]);
    __m256 mi = _mm256_setr_ps(Bi[i].SRC,Bi[i+1].SRC,Bi[i+2].SRC,Bi[i+3].SRC,
      Bi[i+4].SRC,Bi[i+5].SRC,Bi[i+6].SRC,Bi[i+7].SRC);
    __m256 R2 = _mm256_set1_ps(EPS2);

    __m256 x2 = _mm256_set1_ps(Bj[0].X[0]);
    x2 = _mm256_sub_ps(x2, xi);
    __m256 y2 = _mm256_set1_ps(Bj[0].X[1]);
    y2 = _mm256_sub_ps(y2, yi);
    __m256 z2 = _mm256_set1_ps(Bj[0].X[2]);
    z2 = _mm256_sub_ps(z2, zi);
    __m256 mj = _mm256_set1_ps(Bj[0].SRC);

    __m256 xj = x2;
    x2 = _mm256_mul_ps(x2, x2);
    R2 = _mm256_add_ps(R2, x2);
    __m256 yj = y2;
    y2 = _mm256_mul_ps(y2, y2);
    R2 = _mm256_add_ps(R2, y2);
    __m256 zj = z2;
    z2 = _mm256_mul_ps(z2, z2);
    R2 = _mm256_add_ps(R2, z2);

    x2 = _mm256_set1_ps(Bj[1].X[0]);
    y2 = _mm256_set1_ps(Bj[1].X[1]);
    z2 = _mm256_set1_ps(Bj[1].X[2]);
    for (int j=0; j<nj; j++) {
      __m256 invR = _mm256_rsqrt_ps(R2);
      __m256 mask = _mm256_cmp_ps(R2, _mm256_setzero_ps(), _CMP_GT_OQ);
      invR = _mm256_and_ps(invR, mask);
      R2 = _mm256_set1_ps(EPS2);
      x2 = _mm256_sub_ps(x2, xi);
      y2 = _mm256_sub_ps(y2, yi);
      z2 = _mm256_sub_ps(z2, zi);

      mj = _mm256_mul_ps(mj, invR);
      mj = _mm256_mul_ps(mj, mi);
      pot = _mm256_add_ps(pot, mj);
      if (mutual) Bj[j].TRG[0] += vecSum8(mj);
      invR = _mm256_mul_ps(invR, invR);
      invR = _mm256_mul_ps(invR, mj);
      mj = _mm256_set1_ps(Bj[j+1].SRC);

      xj = _mm256_mul_ps(xj, invR);
      ax = _mm256_add_ps(ax, xj);
      if (mutual) Bj[j].TRG[1] -= vecSum8(xj);
      xj = x2;
      x2 = _mm256_mul_ps(x2, x2);
      R2 = _mm256_add_ps(R2, x2);
      x2 = _mm256_set1_ps(Bj[j+2].X[0]);

      yj = _mm256_mul_ps(yj, invR);
      ay = _mm256_add_ps(ay, yj);
      if (mutual) Bj[j].TRG[2] -= vecSum8(yj);
      yj = y2;
      y2 = _mm256_mul_ps(y2, y2);
      R2 = _mm256_add_ps(R2, y2);
      y2 = _mm256_set1_ps(Bj[j+2].X[1]);

      zj = _mm256_mul_ps(zj, invR);
      az = _mm256_add_ps(az, zj);
      if (mutual) Bj[j].TRG[3] -= vecSum8(zj);
      zj = z2;
      z2 = _mm256_mul_ps(z2, z2);
      R2 = _mm256_add_ps(R2, z2);
      z2 = _mm256_set1_ps(Bj[j+2].X[2]);
    }
    for (int k=0; k<8; k++) {
      Bi[i+k].TRG[0] += ((float*)&pot)[k];
      Bi[i+k].TRG[1] += ((float*)&ax)[k];
      Bi[i+k].TRG[2] += ((float*)&ay)[k];
      Bi[i+k].TRG[3] += ((float*)&az)[k];
    }
  }
#endif // __AVX__

#if __SSE__
  for ( ; i<=ni-4; i+=4) {
    __m128 pot = _mm_setzero_ps();
    __m128 ax = _mm_setzero_ps();
    __m128 ay = _mm_setzero_ps();
    __m128 az = _mm_setzero_ps();

    __m128 xi = _mm_setr_ps(Bi[i].X[0], Bi[i+1].X[0], Bi[i+2].X[0], Bi[i+3].X[0]) - _mm_load1_ps(&Xperiodic[0]);
    __m128 yi = _mm_setr_ps(Bi[i].X[1], Bi[i+1].X[1], Bi[i+2].X[1], Bi[i+3].X[1]) - _mm_load1_ps(&Xperiodic[1]);
    __m128 zi = _mm_setr_ps(Bi[i].X[2], Bi[i+1].X[2], Bi[i+2].X[2], Bi[i+3].X[2]) - _mm_load1_ps(&Xperiodic[2]);
    __m128 mi = _mm_setr_ps(Bi[i].SRC,  Bi[i+1].SRC,  Bi[i+2].SRC,  Bi[i+3].SRC);
    __m128 R2 = _mm_set1_ps(EPS2);

    __m128 x2 = _mm_load1_ps(&Bj[0].X[0]);
    x2 = _mm_sub_ps(x2, xi);
    __m128 y2 = _mm_load1_ps(&Bj[0].X[1]);
    y2 = _mm_sub_ps(y2, yi);
    __m128 z2 = _mm_load1_ps(&Bj[0].X[2]);
    z2 = _mm_sub_ps(z2, zi);
    __m128 mj = _mm_load1_ps(&Bj[0].SRC);

    __m128 xj = x2;
    x2 = _mm_mul_ps(x2, x2);
    R2 = _mm_add_ps(R2, x2);
    __m128 yj = y2;
    y2 = _mm_mul_ps(y2, y2);
    R2 = _mm_add_ps(R2, y2);
    __m128 zj = z2;
    z2 = _mm_mul_ps(z2, z2);
    R2 = _mm_add_ps(R2, z2);

    x2 = _mm_load_ps(&Bj[1].X[0]);
    y2 = x2;
    z2 = x2;
    for (int j=0; j<nj; j++) {
      __m128 invR = _mm_rsqrt_ps(R2);
      __m128 mask = _mm_cmpgt_ps(R2, _mm_setzero_ps());
      invR = _mm_and_ps(invR, mask);
      R2 = _mm_set1_ps(EPS2);
      x2 = _mm_shuffle_ps(x2, x2, _MM_SHUFFLE(0,0,0,0));
      x2 = _mm_sub_ps(x2, xi);
      y2 = _mm_shuffle_ps(y2, y2, _MM_SHUFFLE(1,1,1,1));
      y2 = _mm_sub_ps(y2, yi);
      z2 = _mm_shuffle_ps(z2, z2, _MM_SHUFFLE(2,2,2,2));
      z2 = _mm_sub_ps(z2, zi);

      mj = _mm_mul_ps(mj, invR);
      mj = _mm_mul_ps(mj, mi);
      pot = _mm_add_ps(pot, mj);
      if (mutual) Bj[j].TRG[0] += vecSum4(mj);
      invR = _mm_mul_ps(invR, invR);
      invR = _mm_mul_ps(invR, mj);
      mj = _mm_load_ps(&Bj[j+1].X[0]);
      mj = _mm_shuffle_ps(mj, mj, _MM_SHUFFLE(3,3,3,3));

      xj = _mm_mul_ps(xj, invR);
      ax = _mm_add_ps(ax, xj);
      if (mutual) Bj[j].TRG[1] -= vecSum4(xj);
      xj = x2;
      x2 = _mm_mul_ps(x2, x2);
      R2 = _mm_add_ps(R2, x2);
      x2 = _mm_load_ps(&Bj[j+2].X[0]);

      yj = _mm_mul_ps(yj, invR);
      ay = _mm_add_ps(ay, yj);
      if (mutual) Bj[j].TRG[2] -= vecSum4(yj);
      yj = y2;
      y2 = _mm_mul_ps(y2, y2);
      R2 = _mm_add_ps(R2, y2);
      y2 = x2;

      zj = _mm_mul_ps(zj, invR);
      az = _mm_add_ps(az, zj);
      if (mutual) Bj[j].TRG[3] -= vecSum4(zj);
      zj = z2;
      z2 = _mm_mul_ps(z2, z2);
      R2 = _mm_add_ps(R2, z2);
      z2 = x2;
    }
    for (int k=0; k<4; k++) {
      Bi[i+k].TRG[0] += ((float*)&pot)[k];
      Bi[i+k].TRG[1] += ((float*)&ax)[k];
      Bi[i+k].TRG[2] += ((float*)&ay)[k];
      Bi[i+k].TRG[3] += ((float*)&az)[k];
    }
  }
#endif // __SSE__

  for ( ; i<ni; i++) {
    real_t pot = 0;
    vec3 acc = 0;
    for (int j=0; j<nj; j++) {
      vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0f / R2;
        real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
        if (mutual) {
          Bj[j].TRG[0] += invR;
          Bj[j].TRG[1] += dX[0];
          Bj[j].TRG[2] += dX[1];
          Bj[j].TRG[3] += dX[2];
        }
      }
    }
    Bi[i].TRG[0] += pot;
    Bi[i].TRG[1] -= acc[0];
    Bi[i].TRG[2] -= acc[1];
    Bi[i].TRG[3] -= acc[2];
  }
}

void Kernel::P2P(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  int i = 0;
#if __AVX__
  for ( ; i<=n-8; i+=8) {
    __m256 pot = _mm256_setzero_ps();
    __m256 ax = _mm256_setzero_ps();
    __m256 ay = _mm256_setzero_ps();
    __m256 az = _mm256_setzero_ps();

    __m256 xi = _mm256_setr_ps(B[i].X[0],B[i+1].X[0],B[i+2].X[0],B[i+3].X[0],
      B[i+4].X[0],B[i+5].X[0],B[i+6].X[0],B[i+7].X[0]) - _mm256_set1_ps(Xperiodic[0]);
    __m256 yi = _mm256_setr_ps(B[i].X[1],B[i+1].X[1],B[i+2].X[1],B[i+3].X[1],
      B[i+4].X[1],B[i+5].X[1],B[i+6].X[1],B[i+7].X[1]) - _mm256_set1_ps(Xperiodic[1]);
    __m256 zi = _mm256_setr_ps(B[i].X[2],B[i+1].X[2],B[i+2].X[2],B[i+3].X[2],
      B[i+4].X[2],B[i+5].X[2],B[i+6].X[2],B[i+7].X[2]) - _mm256_set1_ps(Xperiodic[2]);
    __m256 mi = _mm256_setr_ps(B[i].SRC,B[i+1].SRC,B[i+2].SRC,B[i+3].SRC,
      B[i+4].SRC,B[i+5].SRC,B[i+6].SRC,B[i+7].SRC);
    __m256 R2 = _mm256_set1_ps(EPS2);

    __m256 x2 = _mm256_set1_ps(B[i+1].X[0]);
    x2 = _mm256_sub_ps(x2, xi);
    __m256 y2 = _mm256_set1_ps(B[i+1].X[1]);
    y2 = _mm256_sub_ps(y2, yi);
    __m256 z2 = _mm256_set1_ps(B[i+1].X[2]);
    z2 = _mm256_sub_ps(z2, zi);
    __m256 mj = _mm256_set1_ps(B[i+1].SRC);

    __m256 xj = x2;
    x2 = _mm256_mul_ps(x2, x2);
    R2 = _mm256_add_ps(R2, x2);
    __m256 yj = y2;
    y2 = _mm256_mul_ps(y2, y2);
    R2 = _mm256_add_ps(R2, y2);
    __m256 zj = z2;
    z2 = _mm256_mul_ps(z2, z2);
    R2 = _mm256_add_ps(R2, z2);

    x2 = _mm256_set1_ps(B[i+2].X[0]);
    y2 = _mm256_set1_ps(B[i+2].X[1]);
    z2 = _mm256_set1_ps(B[i+2].X[2]);
    for (int j=i+1; j<n; j++) {
      __m256 invR = _mm256_rsqrt_ps(R2);
      __m256 mask = _mm256_cmp_ps(_mm256_setr_ps(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7),
        _mm256_set1_ps(j), _CMP_LT_OQ);
      mask = _mm256_and_ps(mask, _mm256_cmp_ps(R2, _mm256_setzero_ps(), _CMP_GT_OQ));
      invR = _mm256_and_ps(invR, mask);
      R2 = _mm256_set1_ps(EPS2);
      x2 = _mm256_sub_ps(x2, xi);
      y2 = _mm256_sub_ps(y2, yi);
      z2 = _mm256_sub_ps(z2, zi);

      mj = _mm256_mul_ps(mj, invR);
      mj = _mm256_mul_ps(mj, mi);
      pot = _mm256_add_ps(pot, mj);
      B[j].TRG[0] += vecSum8(mj);
      invR = _mm256_mul_ps(invR, invR);
      invR = _mm256_mul_ps(invR, mj);
      mj = _mm256_set1_ps(B[j+1].SRC);

      xj = _mm256_mul_ps(xj, invR);
      ax = _mm256_add_ps(ax, xj);
      B[j].TRG[1] -= vecSum8(xj);
      xj = x2;
      x2 = _mm256_mul_ps(x2, x2);
      R2 = _mm256_add_ps(R2, x2);
      x2 = _mm256_set1_ps(B[j+2].X[0]);

      yj = _mm256_mul_ps(yj, invR);
      ay = _mm256_add_ps(ay, yj);
      B[j].TRG[2] -= vecSum8(yj);
      yj = y2;
      y2 = _mm256_mul_ps(y2, y2);
      R2 = _mm256_add_ps(R2, y2);
      y2 = _mm256_set1_ps(B[j+2].X[1]);

      zj = _mm256_mul_ps(zj, invR);
      az = _mm256_add_ps(az, zj);
      B[j].TRG[3] -= vecSum8(zj);
      zj = z2;
      z2 = _mm256_mul_ps(z2, z2);
      R2 = _mm256_add_ps(R2, z2);
      z2 = _mm256_set1_ps(B[j+2].X[2]);
    }
    for (int k=0; k<8; k++) {
      B[i+k].TRG[0] += ((float*)&pot)[k];
      B[i+k].TRG[1] += ((float*)&ax)[k];
      B[i+k].TRG[2] += ((float*)&ay)[k];
      B[i+k].TRG[3] += ((float*)&az)[k];
    }
  }
#endif // __AVX__

#if __SSE__
  for ( ; i<=n-4; i+=4) {
    __m128 pot = _mm_setzero_ps();
    __m128 ax = _mm_setzero_ps();
    __m128 ay = _mm_setzero_ps();
    __m128 az = _mm_setzero_ps();

    __m128 xi = _mm_setr_ps(B[i].X[0], B[i+1].X[0], B[i+2].X[0], B[i+3].X[0]) - _mm_load1_ps(&Xperiodic[0]);
    __m128 yi = _mm_setr_ps(B[i].X[1], B[i+1].X[1], B[i+2].X[1], B[i+3].X[1]) - _mm_load1_ps(&Xperiodic[1]);
    __m128 zi = _mm_setr_ps(B[i].X[2], B[i+1].X[2], B[i+2].X[2], B[i+3].X[2]) - _mm_load1_ps(&Xperiodic[2]);
    __m128 mi = _mm_setr_ps(B[i].SRC,  B[i+1].SRC,  B[i+2].SRC,  B[i+3].SRC);
    __m128 R2 = _mm_set1_ps(EPS2);

    __m128 x2 = _mm_load1_ps(&B[i+1].X[0]);
    x2 = _mm_sub_ps(x2, xi);
    __m128 y2 = _mm_load1_ps(&B[i+1].X[1]);
    y2 = _mm_sub_ps(y2, yi);
    __m128 z2 = _mm_load1_ps(&B[i+1].X[2]);
    z2 = _mm_sub_ps(z2, zi);
    __m128 mj = _mm_load1_ps(&B[i+1].SRC);

    __m128 xj = x2;
    x2 = _mm_mul_ps(x2, x2);
    R2 = _mm_add_ps(R2, x2);
    __m128 yj = y2;
    y2 = _mm_mul_ps(y2, y2);
    R2 = _mm_add_ps(R2, y2);
    __m128 zj = z2;
    z2 = _mm_mul_ps(z2, z2);
    R2 = _mm_add_ps(R2, z2);

    x2 = _mm_load_ps(&B[i+2].X[0]);
    y2 = x2;
    z2 = x2;
    for (int j=i+1; j<n; j++) {
      __m128 invR = _mm_rsqrt_ps(R2);
      __m128 mask = _mm_cmplt_ps(_mm_setr_ps(i, i+1, i+2, i+3), _mm_set1_ps(j));
      mask = _mm_and_ps(mask, _mm_cmpgt_ps(R2, _mm_setzero_ps()));
      invR = _mm_and_ps(invR, mask);
      R2 = _mm_set1_ps(EPS2);
      x2 = _mm_shuffle_ps(x2, x2, _MM_SHUFFLE(0,0,0,0));
      x2 = _mm_sub_ps(x2, xi);
      y2 = _mm_shuffle_ps(y2, y2, _MM_SHUFFLE(1,1,1,1));
      y2 = _mm_sub_ps(y2, yi);
      z2 = _mm_shuffle_ps(z2, z2, _MM_SHUFFLE(2,2,2,2));
      z2 = _mm_sub_ps(z2, zi);

      mj = _mm_mul_ps(mj, invR);
      mj = _mm_mul_ps(mj, mi);
      pot = _mm_add_ps(pot, mj);
      B[j].TRG[0] += vecSum4(mj);
      invR = _mm_mul_ps(invR, invR);
      invR = _mm_mul_ps(invR, mj);
      mj = _mm_load_ps(&B[j+1].X[0]);
      mj = _mm_shuffle_ps(mj, mj, _MM_SHUFFLE(3,3,3,3));

      xj = _mm_mul_ps(xj, invR);
      ax = _mm_add_ps(ax, xj);
      B[j].TRG[1] -= vecSum4(xj);
      xj = x2;
      x2 = _mm_mul_ps(x2, x2);
      R2 = _mm_add_ps(R2, x2);
      x2 = _mm_load_ps(&B[j+2].X[0]);

      yj = _mm_mul_ps(yj, invR);
      ay = _mm_add_ps(ay, yj);
      B[j].TRG[2] -= vecSum4(yj);
      yj = y2;
      y2 = _mm_mul_ps(y2, y2);
      R2 = _mm_add_ps(R2, y2);
      y2 = x2;

      zj = _mm_mul_ps(zj, invR);
      az = _mm_add_ps(az, zj);
      B[j].TRG[3] -= vecSum4(zj);
      zj = z2;
      z2 = _mm_mul_ps(z2, z2);
      R2 = _mm_add_ps(R2, z2);
      z2 = x2;
    }
    for (int k=0; k<4; k++) {
      B[i+k].TRG[0] += ((float*)&pot)[k];
      B[i+k].TRG[1] += ((float*)&ax)[k];
      B[i+k].TRG[2] += ((float*)&ay)[k];
      B[i+k].TRG[3] += ((float*)&az)[k];
    }
  }
#endif // __SSE__

  for ( ; i<n; i++) {
    real_t pot = 0;
    vec3 acc = 0;
    for (int j=i+1; j<n; j++) {
      vec3 dX = B[i].X - B[j].X;
      real_t R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
        B[j].TRG[0] += invR;
        B[j].TRG[1] += dX[0];
        B[j].TRG[2] += dX[1];
        B[j].TRG[3] += dX[2];
      }
    }
    B[i].TRG[0] += pot;
    B[i].TRG[1] -= acc[0];
    B[i].TRG[2] -= acc[1];
    B[i].TRG[3] -= acc[2];
  }
}

void Kernel::P2M(C_iter C, real_t &Rmax) const {
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {
    vec3 dX = B->X - C->X;
    real_t R = std::sqrt(norm(dX));
    if (R > Rmax) Rmax = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,-beta,prefactor,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      for( int m=0; m<=n; ++m ) {
        int nm  = n * n + n + m;
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
  const complex_t I(0.,1.);
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; Cj++) {
    vec3 dX = Ci->X - Cj->X;
    real_t R = std::sqrt(norm(dX)) + Cj->RCRIT;
    if (R > Rmax) Rmax = R;
    real_t rho, alpha, beta;
    cart2sph(rho,alpha,beta,dX);
    evalMultipole(rho,alpha,-beta,prefactor,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        complex_t M = 0;
        for( int n=0; n<=j; ++n ) {
          for( int m=-n; m<=std::min(k-1,n); ++m ) {
            if( j-n >= k-m ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n + m;
              M += Cj->M[jnkms] * std::pow(I,real_t(m-abs(m))) * Ynm[nm]
                 * real_t(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          for( int m=k; m<=n; ++m ) {
            if( j-n >= m-k ) {
              int jnkm  = (j - n) * (j - n) + j - n + k - m;
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n + m;
              M += std::conj(Cj->M[jnkms]) * Ynm[nm]
                 * real_t(ODDEVEN(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        Ci->M[jks] += M * SCALING;
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
  complex_t Ynm[P*P], YnmTheta[P*P];
  vec3 dX = Ci->X - Cj->X - Xperiodic;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalLocal(rho,alpha,beta,prefactor,Ynm,YnmTheta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for( int n=0; n!=P-j; ++n ) {
        for( int m=-n; m<0; ++m ) {
          int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 - m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynm[jnkm];
        }
        for( int m=0; m<=n; ++m ) {
          int nm   = n * n + n + m;
          int nms  = n * (n + 1) / 2 + m;
          int jknm = jk * P * P + nm;
          int jnkm = (j + n) * (j + n) + j + n + m - k;
          L += Cj->M[nms] * Cnm[jknm] * Ynm[jnkm];
        }
      }
      Ci->L[jks] += L;
    }
  }
}

void Kernel::L2L(C_iter Ci) const {
  const complex_t I(0.,1.);
  complex_t Ynm[P*P], YnmTheta[P*P];
  C_iter Cj = Ci0 + Ci->PARENT;
  vec3 dX = Ci->X - Cj->X;
  real_t rho, alpha, beta;
  cart2sph(rho,alpha,beta,dX);
  evalMultipole(rho,alpha,beta,prefactor,Ynm,YnmTheta);
  for( int j=0; j!=P; ++j ) {
    for( int k=0; k<=j; ++k ) {
      int jk = j * j + j + k;
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for( int n=j; n!=P; ++n ) {
        for( int m=j+k-n; m<0; ++m ) {
          int jnkm = (n - j) * (n - j) + n - j + m - k;
          int nm   = n * n + n - m;
          int nms  = n * (n + 1) / 2 - m;
          L += std::conj(Cj->L[nms]) * Ynm[jnkm]
             * real_t(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
        }
        for( int m=0; m<=n; ++m ) {
          if( n-j >= abs(m-k) ) {
            int jnkm = (n - j) * (n - j) + n - j + m - k;
            int nm   = n * n + n + m;
            int nms  = n * (n + 1) / 2 + m;
            L += Cj->L[nms] * std::pow(I,real_t(m-k-abs(m-k)))
               * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
          }
        }
      }
      Ci->L[jks] += L * SCALING;
    }
  }
}

void Kernel::L2P(C_iter Ci) const {
  const complex_t I(0.,1.);
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {
    vec3 dX = B->X - Ci->X;
    vec3 spherical = 0;
    vec3 cartesian = 0;
    real_t r, theta, phi;
    cart2sph(r,theta,phi,dX);
    evalMultipole(r,theta,phi,prefactor,Ynm,YnmTheta);
    B->TRG /= B->SRC;
    for( int n=0; n!=P; ++n ) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
      spherical[0] += std::real(Ci->L[nms] * Ynm[nm]) / r * n;
      spherical[1] += std::real(Ci->L[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
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

void Kernel::preCalculation() {
  const complex_t I(0.,1.);                                   // Imaginary unit
  factorial = new real_t    [P];                              // Factorial
  prefactor = new real_t    [P*P];                            // sqrt( (n - |m|)! / (n + |m|)! )
  Anm       = new real_t    [P*P];                            // (-1)^n / sqrt( (n + m)! / (n - m)! )
  Cnm       = new complex_t [P*P*P*P];                        // M2L translation matrix Cjknm

  factorial[0] = 1;                                           // Initialize factorial
  for (int n=1; n<P; n++) {                                   // Loop to P
    factorial[n] = factorial[n-1] * n;                        //  n!
  }                                                           // End loop to P

  for( int n=0; n!=P; ++n ) {                                 // Loop over n in Anm
    for( int m=-n; m<=n; ++m ) {                              //  Loop over m in Anm
      int nm = n*n+n+m;                                       //   Index of Anm
      int nabsm = abs(m);                                     //   |m|
      real_t fnmm = SCALING;                                  //   Initialize (n - m)!
      for( int i=1; i<=n-m; ++i ) fnmm *= i;                  //   (n - m)!
      real_t fnpm = SCALING;                                  //   Initialize (n + m)!
      for( int i=1; i<=n+m; ++i ) fnpm *= i;                  //   (n + m)!
      real_t fnma = 1.0;                                      //   Initialize (n - |m|)!
      for( int i=1; i<=n-nabsm; ++i ) fnma *= i;              //   (n - |m|)!
      real_t fnpa = 1.0;                                      //   Initialize (n + |m|)!
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
          Cnm[jknm] = std::pow(I,real_t(abs(k-m)-abs(k)-abs(m)))//     Cjknm
                    * real_t(ODDEVEN(j)*Anm[nm]*Anm[jk]/Anm[jnkm]) * SCALING;
        }                                                     //    End loop over m in Cjknm
      }                                                       //   End loop over n in Cjknm
    }                                                         //  End loop over in k in Cjknm
  }                                                           // End loop over in j in Cjknm
}

void Kernel::postCalculation() {
  delete[] factorial;                                         // Free factorial
  delete[] prefactor;                                         // Free sqrt( (n - |m|)! / (n + |m|)! )
  delete[] Anm;                                               // Free (-1)^n / sqrt( (n + m)! / (n - m)! )
  delete[] Cnm;                                               // Free M2L translation matrix Cjknm
}
