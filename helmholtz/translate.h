void getAnm(real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  Anm1[0][0] = 1;
  Anm2[0][0] = 1;
  for (int m=0; m<=P; m++) {
    if(m>0) Anm1[m][m] = sqrt((2 * m - 1.0) / (2 * m));
    if(m<P) Anm1[m+1][m] = sqrt(2 * m + 1.0);
    for (int n=m+2; n<=P; n++) {
      Anm1[n][m] = 2 * n - 1;
      Anm2[n][m] = sqrt((n + m - 1.0) * (n - m - 1.0));
      Anm1[n][m] /= sqrt(real_t(n - m) * (n + m));
      Anm2[n][m] /= sqrt(real_t(n - m) * (n + m));
    }
  }
}

void cart2sph(vec3 dX, real_t & r, real_t & theta, real_t & phi) {
  r = sqrt(norm(dX));
  theta = r == 0 ? 0 : acos(dX[2] / r);
  phi = atan2(dX[1], dX[0]);
}

void get_Ynm(int nterms, real_t x, real_t Ynm[P+1][P+1], real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t y = -sqrt((1 - x) * (1 + x));
  Ynm[0][0] = 1;
  for (int m=0; m<=nterms; m++) {
    if (m > 0) Ynm[m][m] = Ynm[m-1][m-1] * y * Anm1[m][m];
    if (m < nterms) Ynm[m+1][m] = x * Ynm[m][m] * Anm1[m+1][m];
    for (int n=m+2; n<=nterms; n++) {
      Ynm[n][m] = Anm1[n][m] * x * Ynm[n-1][m] - Anm2[n][m] * Ynm[n-2][m];
    }
  }
  for (int n=0; n<=nterms; n++) {
    for (int m=0; m<=n; m++) {
      Ynm[n][m] *= sqrt(2 * n + 1.0);
    }
  }
}

void get_jn(int nterms, complex_t z, real_t scale, complex_t * jn, int ifder, complex_t * jnd) {
  int iscale[nterms+2];
  if (abs(z) < eps) {
    jn[0] = 1;
    for (int i=1; i<=nterms; i++) {
      jn[i] = 0;
    }
    if (ifder) {
      for (int i=0; i<=nterms; i++) {
	jnd[i] = 0;
      }
      jnd[1] = 1.0 / (3 * scale);
    }
    return;
  }
  complex_t zinv = 1.0 / z;
  jn[nterms-1] = 0;
  jn[nterms] = 1;
  real_t coef = 2 * nterms + 1;
  complex_t ztmp = coef * zinv;
  jn[nterms+1] = ztmp;
  int ntop = nterms + 1;
  for (int i=0; i<=ntop; i++) {
    iscale[i] = 0;
  }
  jn[ntop] = 0;
  jn[ntop-1] = 1;
  for (int i=ntop-1; i>0; i--) {
    coef = 2 * i + 1;
    ztmp = coef * zinv * jn[i] - jn[i+1];
    jn[i-1] = ztmp;
    if (abs(ztmp) > 1.0/eps) {
      jn[i] *= eps;
      jn[i-1] *= eps;
      iscale[i] = 1;
    }
  }
  real_t scalinv = 1.0 / scale;
  coef = 1;
  for (int i=0; i<ntop; i++) {
    coef *= scalinv;
    if(iscale[i] == 1) coef *= eps;
    jn[i+1] *= coef;
  }
  complex_t fj0 = sin(z) * zinv;
  complex_t fj1 = fj0 * zinv - cos(z) * zinv;
  if (abs(fj1) > abs(fj0)) {
    ztmp = fj1 / (jn[1] * scale);
  } else {
    ztmp = fj0 / jn[0];
  }
  for (int i=0; i<=nterms; i++) {
    jn[i] *= ztmp;
  }
  if (ifder) {
    jn[nterms+1] *= ztmp;
    jnd[0] = -jn[1] * scale;
    for (int i=1; i<=nterms; i++) {
      coef = i / (2 * i + 1.0);
      jnd[i] = coef * scalinv * jn[i-1] - (1 - coef) * scale * jn[i+1];
    }
  }
}
