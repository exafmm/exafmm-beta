void getAnm() {
  Anm1[0] = 1;
  Anm2[0] = 1;
  for (int m=0; m<=P; m++) {
    int ms = m * (m + 1) / 2 + m;
    int mps = (m + 1) * (m + 2) / 2 + m;
    if(m>0) Anm1[ms] = sqrt((2 * m - 1.0) / (2 * m));
    if(m<P) Anm1[mps] = sqrt(2 * m + 1.0);
    for (int n=m+2; n<=P; n++) {
      int nms = n * (n + 1) / 2 + m;
      Anm1[nms] = 2 * n - 1;
      Anm2[nms] = sqrt((n + m - 1.0) * (n - m - 1.0));
      Anm1[nms] /= sqrt(real_t(n - m) * (n + m));
      Anm2[nms] /= sqrt(real_t(n - m) * (n + m));
    }
  }
}

void cart2sph(vec3 dX, real_t & r, real_t & theta, real_t & phi) {
  r = sqrt(norm(dX));
  theta = r == 0 ? 0 : acos(dX[2] / r);
  phi = atan2(dX[1], dX[0]);
}

void rotate(real_t theta, int nterms, complex_t Mnm[(P+1)*(P+1)],
	    complex_t Mrot[(P+1)*(P+1)]) {
  real_t Rnm1[P+1][2*P+1];
  real_t Rnm2[P+1][2*P+1];
  real_t sqrtCnm[2*P+1][2];
  for (int m=0; m<=2*nterms; m++) {
    sqrtCnm[m][0] = sqrt(m+0.0);
  }
  sqrtCnm[0][1] = 0;
  sqrtCnm[1][1] = 0;
  for (int m=2; m<=2*nterms; m++) {
    sqrtCnm[m][1] = sqrt(m * (m - 1) / 2.0);
  }
  real_t ctheta = cos(theta);
  if (fabs(ctheta) < eps) ctheta = 0;
  real_t stheta = sin(-theta);
  if (fabs(stheta) < eps) stheta = 0;
  real_t hsthta = stheta / sqrt(2.0);
  real_t cthtap = sqrt(2.0) * cos(theta * .5) * cos(theta * .5);
  real_t cthtan =-sqrt(2.0) * sin(theta * .5) * sin(theta * .5);
  Rnm1[0][P] = 1;
  Mrot[0] = Mnm[0] * Rnm1[0][P];
  for (int n=1; n<=nterms; n++) {
    for (int m=-n; m<0; m++) {
      Rnm2[0][P+m] = -sqrtCnm[n-m][1] * Rnm1[0][P+m+1];
      if (m > (1 - n)) {
	Rnm2[0][P+m] += sqrtCnm[n+m][1] * Rnm1[0][P+m-1];
      }
      Rnm2[0][P+m] *= hsthta; 
      if (m > -n) {
	Rnm2[0][P+m] += Rnm1[0][P+m] * ctheta * sqrtCnm[n+m][0] * sqrtCnm[n-m][0];
      }
      Rnm2[0][P+m] /= n;
    }
    Rnm2[0][P] = Rnm1[0][P] * ctheta;
    if (n > 1) {
      Rnm2[0][P] += hsthta * sqrtCnm[n][1] * (2 * Rnm1[0][P-1]) / n;
    }
    for (int m=1; m<=n; m++) {
      Rnm2[0][P+m] = Rnm2[0][P-m];
      if (m % 2 == 0) {
	Rnm2[m][P] = Rnm2[0][P+m];
      } else {
	Rnm2[m][P] =-Rnm2[0][P+m];
      }
    }
    for (int mp=1; mp<=n; mp++) {
      real_t scale = 1 / (sqrt(2.0) * sqrtCnm[n+mp][1]);
      for (int m=mp; m<=n; m++) {
	Rnm2[mp][P+m] = Rnm1[mp-1][P+m-1] * cthtap * sqrtCnm[n+m][1];
	Rnm2[mp][P-m] = Rnm1[mp-1][P-m+1] * cthtan * sqrtCnm[n+m][1];
	if (m < (n - 1)) {
	  Rnm2[mp][P+m] -= Rnm1[mp-1][P+m+1] * cthtan * sqrtCnm[n-m][1];
	  Rnm2[mp][P-m] -= Rnm1[mp-1][P-m-1] * cthtap * sqrtCnm[n-m][1];
	}
	if (m < n) {
	  real_t d = stheta * sqrtCnm[n+m][0] * sqrtCnm[n-m][0];
	  Rnm2[mp][P+m] += Rnm1[mp-1][P+m] * d;
	  Rnm2[mp][P-m] += Rnm1[mp-1][P-m] * d;
	}
	Rnm2[mp][P+m] *= scale;
	Rnm2[mp][P-m] *= scale;
	if (m > mp) {
	  if ((mp+m) % 2 == 0) {
	    Rnm2[m][P+mp] = Rnm2[mp][P+m];
	    Rnm2[m][P-mp] = Rnm2[mp][P-m];
	  } else {
	    Rnm2[m][P+mp] =-Rnm2[mp][P+m];
	    Rnm2[m][P-mp] =-Rnm2[mp][P-m];
	  }
	}
      }
    }
    for (int m=-n; m<=n; m++) {
      int nn = n * n + n;
      int nm = n * n + n + m;
      Mrot[nm] = Mnm[nn] * Rnm2[0][P+m];
      for (int mp=1; mp<=n; mp++) {
	nm = n * n + n + m;
	int npm = n * n + n + mp;
	int nmm = n * n + n - mp;
	Mrot[nm] += Mnm[npm] * Rnm2[mp][P+m] + Mnm[nmm] * Rnm2[mp][P-m];
      }
    }
    for (int m=-n; m<=n; m++) {
      for (int mp=0; mp<=n; mp++) {
	Rnm1[mp][P+m] = Rnm2[mp][P+m];
      }
    }    
  }
}

void get_Ynm(int nterms, real_t x, real_t Ynm[(P+1)*(P+2)/2]) {
  real_t y = -sqrt((1 - x) * (1 + x));
  Ynm[0] = 1;
  for (int m=0; m<=nterms; m++) {
    int ms = m * (m + 1) / 2 + m;
    int mms = m * (m - 1) / 2 + m - 1;
    int mps = (m + 1) * (m + 2) / 2 + m;
    if (m > 0) Ynm[ms] = Ynm[mms] * y * Anm1[ms];
    if (m < nterms) Ynm[mps] = x * Ynm[ms] * Anm1[mps];
    for (int n=m+2; n<=nterms; n++) {
      int nms = n * (n + 1) / 2 + m;
      int nm1 = n * (n - 1) / 2 + m;
      int nm2 = (n - 1) * (n - 2) / 2 + m;
      Ynm[nms] = Anm1[nms] * x * Ynm[nm1] - Anm2[nms] * Ynm[nm2];
    }
  }
  for (int n=0; n<=nterms; n++) {
    for (int m=0; m<=n; m++) {
      int nms = n * (n + 1) / 2 + m;
      Ynm[nms] *= sqrt(2 * n + 1.0);
    }
  }
}

void get_Ynmd(int nterms, real_t x, real_t Ynm[(P+1)*(P+2)/2], real_t Ynmd[(P+1)*(P+2)/2]) {
  real_t y = -sqrt((1 - x) * (1 + x));
  real_t y2 = y * y;
  Ynm[0] = 1;
  Ynmd[0] = 0;
  Ynm[1] = x * Ynm[0] * Anm1[1];
  Ynmd[1] = (x * Ynmd[0] + Ynm[0]) * Anm1[1];
  for (int n=2; n<=nterms; n++) {
    int ns = n * (n + 1) / 2;
    int nm1 = n * (n - 1) / 2;
    int nm2 = (n - 1) * (n - 2) / 2;
    Ynm[ns] = Anm1[ns] * x * Ynm[nm1] - Anm2[ns] * Ynm[nm2];
    Ynmd[ns] = Anm1[ns] * (x * Ynmd[nm1] + Ynm[nm1]) - Anm2[ns] * Ynmd[nm2];
  }
  for (int m=1; m<=nterms; m++) {
    int ms = m * (m + 1) / 2 + m;
    int mms = m * (m - 1) / 2 + m - 1;
    int mps = (m + 1) * (m + 2) / 2 + m;
    if (m == 1) Ynm[ms] = -Ynm[mms] * Anm1[ms];
    if (m > 1) Ynm[ms] = Ynm[mms] * y * Anm1[ms];
    if (m > 0) Ynmd[ms] = -Ynm[ms] * m * x;
    if (m < nterms) Ynm[mps] = x * Ynm[ms] * Anm1[mps];
    if (m < nterms) Ynmd[mps] = (x * Ynmd[ms] + y2 * Ynm[ms]) * Anm1[mps];
    for (int n=m+2; n<=nterms; n++) {
      int nms = n * (n + 1) / 2 + m;
      int nm1 = n * (n - 1) / 2 + m;
      int nm2 = (n - 1) * (n - 2) / 2 + m;
      Ynm[nms] = Anm1[nms] * x * Ynm[nm1] - Anm2[nms] * Ynm[nm2];
      Ynmd[nms] = Anm1[nms] * (x * Ynmd[nm1] + y2 * Ynm[nm1]) - Anm2[nms] * Ynmd[nm2];
    }
  }
  for (int n=0; n<=nterms; n++) {
    for (int m=0; m<=n; m++) {
      int nms = n * (n + 1) / 2 + m;
      Ynm[nms] *= sqrt(2 * n + 1.0);
      Ynmd[nms] *= sqrt(2 * n + 1.0);
    }
  }
}

void get_hn(int nterms, complex_t z, real_t scale, complex_t * hn) {
  if (abs(z) < eps) {
    for (int i=0; i<=nterms; i++) {
      hn[i] = 0;
    }
    return;
  }
  complex_t zi = I * z;
  complex_t zinv = scale / z;
  hn[0] = exp(zi) / zi;
  hn[1] = hn[0] * (zinv - I * scale);
  real_t scale2 = scale * scale;
  for (int i=2; i<=nterms; i++) {
    hn[i] = zinv * (2 * i - 1.0) * hn[i-1] - scale2 * hn[i-2];
  }
}

void get_hnd(int nterms, complex_t z, real_t scale, complex_t * hn, complex_t * hnd) {
  if (abs(z) < eps) {
    for (int i=0; i<=nterms; i++) {
      hn[i] = 0;
      hnd[i] = 0;
    }
    return;
  }
  complex_t zi = I * z;
  complex_t zinv = 1.0 / z;
  hn[0] = exp(zi) / zi;
  hn[1] = hn[0] * (zinv - I) * scale;
  hnd[0] = -hn[1] / scale;
  hnd[1] = -zinv * 2.0 * hn[1] + scale * hn[0];
  for (int i=2; i<=nterms; i++) {
    hn[i] = (zinv * (2 * i - 1.0) * hn[i-1] - scale * hn[i-2]) * scale;
    hnd[i] = -zinv * (i + 1.0) * hn[i] + scale * hn[i-1];
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
  for (int i=1; i<=ntop; i++) {
    coef *= scalinv;
    if(iscale[i-1] == 1) coef *= eps;
    jn[i] *= coef;
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

void polynomial(real_t x, int n, real_t & pol, real_t & der, real_t & sum) {
  sum = 0.5 + x * x * 1.5;
  real_t pk = 1;
  real_t pkp1 = x;
  if (n < 2) {
    der = 0;
    sum = 0.5;
    if (n == 0) return;
    der = 1;
    sum += x * x * 1.5;
    return;
  }
  for (int k=1; k<n; k++) {
    real_t pkm1 = pk;
    pk = pkp1;
    pkp1 = ((2 * k + 1) * x * pk - k * pkm1) / (k + 1);
    sum += pkp1 * pkp1 * (k + 1.5);
  }
  pol = pkp1;
  der = n * (x * pkp1 - pk) / (x * x - 1);
}

void legendre(int nquad, real_t xquad[2*P], real_t wquad[2*P]) {
  real_t pol = 0, der, sum;
  real_t h = M_PI / (2 * nquad);
  for (int i=1; i<=nquad; i++) {
    xquad[nquad-i] = cos((2 * i - 1) * h);
  }
  xquad[nquad/2] = 0;
  for (int i=0; i<nquad/2; i++) {
    real_t xk = xquad[i];
    int ifout = 0;
    for (int k=0; k<10; k++) {
      polynomial(xk,nquad,pol,der,sum);
      real_t delta = -pol / der;
      xk += delta;
      if (fabs(delta) < eps) ifout++;
      if (ifout == 3) break;
    }
    xquad[i] = xk;
    xquad[nquad-i-1] = -xk;
  }
  for (int i=0; i<(nquad+1)/2; i++) {
    polynomial(xquad[i],nquad,pol,der,sum);
    wquad[i] = 1 / sum;
    wquad[nquad-i-1] = wquad[i];
  }
}
