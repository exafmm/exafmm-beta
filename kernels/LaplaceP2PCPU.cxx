#include "kernel.h"

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
