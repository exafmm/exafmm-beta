#include "kernel.h"
#include "simdvec.h"

real_t kernel::eps2;
vec3 kernel::Xperiodic;

void kernel::P2P(C_iter Ci, C_iter Cj, bool mutual) {
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NBODY;
  int nj = Cj->NBODY;
  int i = 0;
#if USE_SIMD
  for ( ; i<=ni-NSIMD; i+=NSIMD) {
    simdvec zero = 0.0;
    ksimdvec pot = zero;
    ksimdvec ax = zero;
    ksimdvec ay = zero;
    ksimdvec az = zero;

    simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(Bi,i);
    simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(Bi,i);
    simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(Bi,i);
    simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(Bi,i);
    simdvec R2 = eps2;

    simdvec xj = Xperiodic[0];
    xi -= xj;
    simdvec yj = Xperiodic[1];
    yi -= yj;
    simdvec zj = Xperiodic[2];
    zi -= zj;

    simdvec x2 = Bj[0].X[0];
    x2 -= xi;
    simdvec y2 = Bj[0].X[1];
    y2 -= yi;
    simdvec z2 = Bj[0].X[2];
    z2 -= zi;
    simdvec mj = Bj[0].SRC;

    xj = x2;
    R2 += x2 * x2;
    yj = y2;
    R2 += y2 * y2;
    zj = z2;
    R2 += z2 * z2;
    simdvec invR;

    x2 = Bj[1].X[0];
    y2 = Bj[1].X[1];
    z2 = Bj[1].X[2];
    for (int j=0; j<nj-2; j++) {
      invR = rsqrt(R2); 
      invR &= R2 > zero;
      R2 = eps2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      if (mutual) Bj[j].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = Bj[j+1].SRC;

      xj *= invR;
      ax += xj;
      if (mutual) Bj[j].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;
      x2 = Bj[j+2].X[0];

      yj *= invR;
      ay += yj;
      if (mutual) Bj[j].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;
      y2 = Bj[j+2].X[1];

      zj *= invR;
      az += zj;
      if (mutual) Bj[j].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
      z2 = Bj[j+2].X[2];
    }
    if ( nj > 1 ) {
      invR = rsqrt(R2);
      invR &= R2 > zero;
      R2 = eps2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      if (mutual) Bj[nj-2].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = Bj[nj-1].SRC;

      xj *= invR;
      ax += xj;
      if (mutual) Bj[nj-2].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;

      yj *= invR;
      ay += yj;
      if (mutual) Bj[nj-2].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;

      zj *= invR;
      az += zj;
      if (mutual) Bj[nj-2].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
    }
    invR = rsqrt(R2);
    invR &= R2 > zero;
    mj *= invR * mi;
    pot += mj;
    if (mutual) Bj[nj-1].TRG[0] += sum(mj);
    invR = invR * invR * mj;

    xj *= invR;
    ax += xj;
    if (mutual) Bj[nj-1].TRG[1] -= sum(xj);
    yj *= invR;
    ay += yj;
    if (mutual) Bj[nj-1].TRG[2] -= sum(yj);
    zj *= invR;
    az += zj;
    if (mutual) Bj[nj-1].TRG[3] -= sum(zj);
    for (int k=0; k<NSIMD; k++) {
      Bi[i+k].TRG[0] += transpose(pot,k);
      Bi[i+k].TRG[1] += transpose(ax,k);
      Bi[i+k].TRG[2] += transpose(ay,k);
      Bi[i+k].TRG[3] += transpose(az,k);
    }
  }
#endif
  for ( ; i<ni; i++) {
    kreal_t pot = 0; 
    kreal_t ax = 0;
    kreal_t ay = 0;
    kreal_t az = 0;
    for (int j=0; j<nj; j++) {
      vec3 dX = Bi[i].X - Bj[j].X - Xperiodic;
      real_t R2 = norm(dX) + eps2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = Bi[i].SRC * Bj[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        ax += dX[0];
        ay += dX[1];
        az += dX[2];
        if (mutual) {
          Bj[j].TRG[0] += invR;
          Bj[j].TRG[1] += dX[0];
          Bj[j].TRG[2] += dX[1];
          Bj[j].TRG[3] += dX[2];
        }
      }
    }
    Bi[i].TRG[0] += pot;
    Bi[i].TRG[1] -= ax;
    Bi[i].TRG[2] -= ay;
    Bi[i].TRG[3] -= az;
  }
}

void kernel::P2P(C_iter C) {
  B_iter B = C->BODY;
  int n = C->NBODY;
  int i = 0;
#if USE_SIMD
  for ( ; i<=n-NSIMD; i+=NSIMD) {
    simdvec zero = 0;
    ksimdvec pot = zero;
    ksimdvec ax = zero;
    ksimdvec ay = zero;
    ksimdvec az = zero;

    simdvec index = SIMD<simdvec,0,NSIMD>::setIndex(i);
    simdvec xi = SIMD<simdvec,0,NSIMD>::setBody(B,i);
    simdvec yi = SIMD<simdvec,1,NSIMD>::setBody(B,i);
    simdvec zi = SIMD<simdvec,2,NSIMD>::setBody(B,i);
    simdvec mi = SIMD<simdvec,3,NSIMD>::setBody(B,i);
    simdvec R2 = eps2;

    simdvec x2 = B[i+1].X[0];
    x2 -= xi;
    simdvec y2 = B[i+1].X[1];
    y2 -= yi;
    simdvec z2 = B[i+1].X[2];
    z2 -= zi;
    simdvec mj = B[i+1].SRC;

    simdvec xj = x2;
    R2 += x2 * x2;
    simdvec yj = y2;
    R2 += y2 * y2;
    simdvec zj = z2;
    R2 += z2 * z2;
    simdvec invR;

    x2 = B[i+2].X[0];
    y2 = B[i+2].X[1];
    z2 = B[i+2].X[2];
    for (int j=i+1; j<n-2; j++) {
      invR = rsqrt(R2);
      invR &= index < j;
      invR &= R2 > zero;
      R2 = eps2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      B[j].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = B[j+1].SRC;

      xj *= invR;
      ax += xj;
      B[j].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;
      x2 = B[j+2].X[0];

      yj *= invR;
      ay += yj;
      B[j].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;
      y2 = B[j+2].X[1];

      zj *= invR;
      az += zj;
      B[j].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
      z2 = B[j+2].X[2];
    }
    if ( n-2 > i ) {
      invR = rsqrt(R2);
      invR &= index < n-2;
      invR &= R2 > zero;
      R2 = eps2;
      x2 -= xi;
      y2 -= yi;
      z2 -= zi;

      mj *= invR * mi;
      pot += mj;
      B[n-2].TRG[0] += sum(mj);
      invR = invR * invR * mj;
      mj = B[n-1].SRC;

      xj *= invR;
      ax += xj;
      B[n-2].TRG[1] -= sum(xj);
      xj = x2;
      R2 += x2 * x2;

      yj *= invR;
      ay += yj;
      B[n-2].TRG[2] -= sum(yj);
      yj = y2;
      R2 += y2 * y2;

      zj *= invR;
      az += zj;
      B[n-2].TRG[3] -= sum(zj);
      zj = z2;
      R2 += z2 * z2;
    }
    invR = rsqrt(R2);
    invR &= index < n-1;
    invR &= R2 > zero;
    mj *= invR * mi;
    pot += mj;
    B[n-1].TRG[0] += sum(mj);
    invR = invR * invR * mj;

    xj *= invR;
    ax += xj;
    B[n-1].TRG[1] -= sum(xj);
    yj *= invR;
    ay += yj;
    B[n-1].TRG[2] -= sum(yj);
    zj *= invR;
    az += zj;
    B[n-1].TRG[3] -= sum(zj);
    for (int k=0; k<NSIMD; k++) {
      B[i+k].TRG[0] += transpose(pot,k);
      B[i+k].TRG[1] += transpose(ax,k);
      B[i+k].TRG[2] += transpose(ay,k);
      B[i+k].TRG[3] += transpose(az,k);
    }
  }
#endif
  for ( ; i<n; i++) {
    kreal_t pot = 0;
    kreal_t ax = 0;
    kreal_t ay = 0;
    kreal_t az = 0;
    for (int j=i+1; j<n; j++) {
      vec3 dX = B[i].X - B[j].X;
      real_t R2 = norm(dX) + eps2;
      if (R2 != 0) {
        real_t invR2 = 1.0 / R2;
        real_t invR = B[i].SRC * B[j].SRC * sqrt(invR2);
        dX *= invR2 * invR;
        pot += invR;
        ax += dX[0];
        ay += dX[1];
        az += dX[2];
        B[j].TRG[0] += invR;
        B[j].TRG[1] += dX[0];
        B[j].TRG[2] += dX[1];
        B[j].TRG[3] += dX[2];
      }
    }
    B[i].TRG[0] += pot;
    B[i].TRG[1] -= ax;
    B[i].TRG[2] -= ay;
    B[i].TRG[3] -= az;
  }
}
