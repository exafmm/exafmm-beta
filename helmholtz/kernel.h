void P2P(int * icell, complex_t * pi, cvec3 * Fi, int * jcell, vec3 * Xj, complex_t * qj, complex_t wavek) {
  for (int i=icell[7]; i<icell[7]+icell[8]; i++) {
    complex_t p = 0.0;
    complex_t F[3] = {0.0,0.0,0.0};
    for (int j=jcell[7]; j<jcell[7]+jcell[8]; j++) {
      vec3 dX = Xj[i] - Xj[j];
      real_t R2 = norm(dX);
      if (R2 != 0) {
	real_t R = sqrt(R2);
	complex_t coef1 = qj[j] * exp(I * wavek * R) / R;
	complex_t coef2 = (1.0 - I * wavek * R) * coef1 / R2;
	p += coef1;
	F[0] += coef2 * dX[0];
	F[1] += coef2 * dX[1];
	F[2] += coef2 * dX[2];
      }
    }
    pi[i] += p;
    Fi[i][0] += F[0];
    Fi[i][1] += F[1];
    Fi[i][2] += F[2];
  }
}

void P2M(complex_t wavek, real_t scale, vec3 * Xj, complex_t * qj, int nj, vec3 Xi, complex_t Mi[P+1][2*P+1],
	 real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1];
  complex_t ephi[P+1], jn[P+2], jnd[P+2], Mnm[(P+1)*(P+1)];
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Mnm[nm] = 0;
    }
  }
  for (int i=0; i<nj; i++) {
    vec3 dX = Xj[i] - Xi;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    real_t ctheta = cos(theta);
    ephi[1] = exp(I * phi);
    for (int n=2; n<=P; n++) {
      ephi[n] = ephi[n-1] * ephi[1];
    }
    get_Ynm(P, ctheta, Ynm, Anm1, Anm2);
    complex_t z = wavek * r;
    get_jn(P, z, scale, jn, 0, jnd);
    for (int n=0; n<=P; n++) {
      jn[n] *= qj[i];
    }
    for (int n=0; n<=P; n++) {
      int nm = n * n + n;
      Mnm[nm] += Ynm[n][0] * jn[n];
      for (int m=1; m<=n; m++) {
	int npm = n * n + n + m;
	int nmm = n * n + n - m;
	complex_t Ynmjn = Ynm[n][m] * jn[n];
	Mnm[npm] += Ynmjn * conj(ephi[m]);
	Mnm[nmm] += Ynmjn * ephi[m];
      }
    }
  }
  for (int n=0; n<=P; n++) {
    int nm = n * n + n;
    Mi[n][P] += Mnm[nm] * I * wavek;
    for (int m=1; m<=n; m++) {
      int npm = n * n + n + m;
      int nmm = n * n + n - m;
      Mi[n][P+m] += Mnm[npm] * I * wavek;
      Mi[n][P-m] += Mnm[nmm] * I * wavek;
    }
  }
}

void M2M(complex_t wavek, real_t scalej, vec3 Xj, complex_t Mj[P+1][2*P+1],
	 real_t scalei, vec3 Xi, complex_t Mi[P+1][2*P+1],
	 real_t radius, real_t xquad[2*P], real_t wquad[2*P], int nquad,
	 real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1];
  complex_t Mnm[P+1][2*P+1];
  complex_t Mrot[P+1][2*P+1];
  complex_t phitemp[nquad][2*P+1];
  complex_t hn[P+1];
  complex_t ephi[2*P+1];
  vec3 dX = Xi - Xj;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(I * phi);
  ephi[P] = 1;
  ephi[P-1] = conj(ephi[P+1]);
  for (int n=2; n<=P; n++) {
    ephi[P+n] = ephi[P+n-1] * ephi[P+1];
    ephi[P-n] = conj(ephi[P+n]);
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = Mj[n][P+m] * ephi[P+m];
    }
  }
  rotate(theta, P, Mnm, Mrot);
  for (int l=0; l<nquad; l++) {
    for (int m=-P; m<=P; m++) {
      phitemp[l][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    real_t ctheta = xquad[l];
    real_t stheta = sqrt(1 - ctheta * ctheta);
    real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
    rj = sqrt(rj);
    real_t cthetaj = (r + radius * ctheta) / rj;
    complex_t z = wavek * rj;
    get_Ynm(P, cthetaj, Ynm, Anm1, Anm2);
    get_hn(P, z, scalej, hn);
    for (int m=-P; m<=P; m++) {
      int mabs = abs(m);
      for (int n=mabs; n<=P; n++) {
	phitemp[l][P+m] += Mrot[n][P+m] * hn[n] * Ynm[n][mabs];
      }
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    get_Ynm(P, xquad[l], Ynm, Anm1, Anm2);
    for (int m=-P; m<=P; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[l][P+m] * wquad[l] * .5;
      for (int n=mabs; n<=P; n++) {
	Mnm[n][P+m] += z * Ynm[n][mabs];
      }
    }
  }
  complex_t z = wavek * radius;
  get_hn(P, z, scalei, hn);
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] /= hn[n];
    }
  }
  rotate(-theta, P, Mnm, Mrot);
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = ephi[P-m] * Mrot[n][P+m];
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mi[n][P+m] += Mnm[n][P+m];
    }
  }
}

void M2L(complex_t wavek, real_t scalej, vec3 Xj, complex_t Mj[P+1][2*P+1],
	 real_t scalei, vec3 Xi, complex_t Li[P+1][2*P+1], int Popt, real_t radius,
	 real_t xquad[2*P], real_t wquad[2*P], int nquad,
	 real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1], Ynmd[P+1][P+1];
  complex_t phitemp[nquad][2*Popt+1], phitempn[nquad][2*Popt+1];
  complex_t hn[P+1], hnd[P+1], jn[P+2], jnd[P+2], ephi[2*P+1];
  complex_t Mnm[P+1][2*P+1], Mrot[P+1][2*P+1];
  complex_t Lnm[P+1][2*P+1], Lrot[P+1][2*P+1], Lnmd[P+1][2*P+1];
  vec3 dX = Xi - Xj;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(I * phi);
  ephi[P] = 1;
  ephi[P-1] = conj(ephi[P+1]);
  for (int n=2; n<=P; n++) {
    ephi[P+n] = ephi[P+n-1] * ephi[P+1];
    ephi[P-n] = conj(ephi[P+n]);
  }
  for (int n=0; n<=Popt; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = Mj[n][P+m] * ephi[P+m];
      Lnm[n][P+m] = 0;
    }
  }
  rotate(theta, Popt, Mnm, Mrot);
  for (int l=0; l<nquad; l++) {
    for (int m=-Popt; m<=Popt; m++) {
      phitemp[l][Popt+m] = 0;
      phitempn[l][Popt+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    real_t ctheta = xquad[l];
    real_t stheta = sqrt(1 - ctheta * ctheta);
    real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
    rj = sqrt(rj);
    real_t cthetaj = (r + radius * ctheta) / rj;
    real_t sthetaj = sqrt(1 - cthetaj * cthetaj);
    real_t rn = sthetaj * stheta + cthetaj * ctheta;
    real_t thetan = (cthetaj * stheta - ctheta * sthetaj) / rj;
    complex_t z = wavek * rj;
    get_Ynmd(Popt, cthetaj, Ynm, Ynmd, Anm1, Anm2);
    get_hnd(Popt, z, scalej, hn, hnd);
    for (int n=0; n<=Popt; n++) {
      hnd[n] *= wavek;
    }
    for (int n=1; n<=Popt; n++) {
      for (int m=1; m<=n; m++) {
	Ynm[n][m] *= sthetaj;
      }
    }
    phitemp[l][Popt] = Mrot[0][P] * hn[0];
    phitempn[l][Popt] = Mrot[0][P] * hnd[0] * rn;
    for (int n=1; n<=Popt; n++) {
      phitemp[l][Popt] += Mrot[n][P] * hn[n] * Ynm[n][0];
      complex_t ut1 = hnd[n] * rn;
      complex_t ut2 = hn[n] * thetan;
      complex_t ut3 = ut1 * Ynm[n][0] - ut2 * Ynmd[n][0] * sthetaj;
      phitempn[l][Popt] += ut3 * Mrot[n][P];
      for (int m=1; m<=n; m++) {
	z = hn[n] * Ynm[n][m];
	phitemp[l][Popt+m] += Mrot[n][P+m] * z;
	phitemp[l][Popt-m] += Mrot[n][P-m] * z;
	ut3 = ut1 * Ynm[n][m] - ut2 * Ynmd[n][m];
	phitempn[l][Popt+m] += ut3 * Mrot[n][P+m];
	phitempn[l][Popt-m] += ut3 * Mrot[n][P-m];
      }
    }
  }
  for (int n=0; n<=Popt; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = 0;
      Lnmd[n][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    real_t cthetaj = xquad[l];
    get_Ynm(Popt, cthetaj, Ynm, Anm1, Anm2);
    for (int m=-Popt; m<=Popt; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[l][Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<=Popt; n++) {
	Lnm[n][P+m] += z * Ynm[n][mabs];
      }
      z = phitempn[l][Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<=Popt; n++) {
	Lnmd[n][P+m] += z * Ynm[n][mabs];
      }
    }
  }
  complex_t z = wavek * radius;
  get_jn(Popt, z, scalei, jn, 1, jnd);
  for (int n=0; n<=Popt; n++) {
    for (int m=-n; m<=n; m++) {
      complex_t zh = jn[n];
      complex_t zhn = jnd[n] * wavek;
      complex_t z = zh * zh + zhn * zhn;
      Lnm[n][P+m] = (zh * Lnm[n][P+m] + zhn * Lnmd[n][P+m]) / z;
    }
  }
  rotate(-theta, Popt, Lnm, Lrot);
  for (int n=0; n<=Popt; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = ephi[P-m] * Lrot[n][P+m];
    }
  }
  for (int n=0; n<=Popt; n++) {
    for (int m=-n; m<=n; m++) {
      Li[n][P+m] += Lnm[n][P+m];
    }
  }
}

void L2L(complex_t wavek, real_t scalej, vec3 Xj, complex_t Lj[P+1][2*P+1],
         real_t scalei, vec3 Xi, complex_t Li[P+1][2*P+1],
         real_t radius, real_t xquad[2*P], real_t wquad[2*P], int nquad,
         real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1], Ynmd[P+1][P+1];
  complex_t Lnm[P+1][2*P+1];
  complex_t Lnmd[P+1][2*P+1];
  complex_t Lrot[P+1][2*P+1];
  complex_t phitemp[nquad][2*P+1], phitempn[nquad][2*P+1];
  complex_t jn[P+2], jnd[P+2], ephi[2*P+1];
  vec3 dX = Xi - Xj;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(I * phi);
  ephi[P] = 1;
  ephi[P-1] = conj(ephi[P+1]);
  for (int n=2; n<=P; n++) {
    ephi[P+n] = ephi[P+n-1] * ephi[P+1];
    ephi[P-n] = conj(ephi[P+n]);
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = Lj[n][P+m] * ephi[P+m];
    }
  }
  rotate(theta, P, Lnm, Lrot);
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    for (int m=-P; m<=P; m++) {
      phitemp[l][P+m] = 0;
      phitempn[l][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    real_t ctheta = xquad[l];
    real_t stheta = sqrt(1 - ctheta * ctheta);
    real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
    rj = sqrt(rj);
    real_t cthetaj = (r + radius * ctheta) / rj;
    real_t sthetaj = sqrt(1 - cthetaj * cthetaj);
    real_t rn = sthetaj * stheta + cthetaj * ctheta;
    real_t thetan = (cthetaj * stheta - ctheta * sthetaj) / rj;
    complex_t z = wavek * rj;
    get_Ynmd(P, cthetaj, Ynm, Ynmd, Anm1, Anm2);
    get_jn(P, z, scalej, jn, 1, jnd);
    for (int n=0; n<=P; n++) {
      jnd[n] *= wavek;
    }
    for (int n=1; n<=P; n++) {
      for (int m=1; m<=n; m++) {
	Ynm[n][m] *= sthetaj;
      }
    }
    phitemp[l][P] = Lrot[0][P] * jn[0];
    phitempn[l][P] = Lrot[0][P] * jnd[0] * rn;
    for (int n=1; n<=P; n++) {
      phitemp[l][P] += Lrot[n][P] * jn[n] * Ynm[n][0];
      complex_t ut1 = jnd[n] * rn;
      complex_t ut2 = jn[n] * thetan;
      complex_t ut3 = ut1 * Ynm[n][0] - ut2 * Ynmd[n][0] * sthetaj;
      phitempn[l][P] += ut3 * Lrot[n][P];
      for (int m=1; m<=n; m++) {
	z = jn[n] * Ynm[n][m];
	phitemp[l][P+m] += Lrot[n][P+m] * z;
	phitemp[l][P-m] += Lrot[n][P-m] * z;
	ut3 = ut1 * Ynm[n][m] - ut2 * Ynmd[n][m];
	phitempn[l][P+m] += ut3 * Lrot[n][P+m];
	phitempn[l][P-m] += ut3 * Lrot[n][P-m];
      }
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = 0;
      Lnmd[n][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    real_t cthetaj = xquad[l];
    get_Ynm(P, cthetaj, Ynm, Anm1, Anm2);
    for (int m=-P; m<=P; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[l][P+m] * wquad[l] * .5;
      for (int n=mabs; n<=P; n++) {
        Lnm[n][P+m] += z * Ynm[n][mabs];
      }
      z = phitempn[l][P+m] * wquad[l] * .5;
      for (int n=mabs; n<=P; n++) {
        Lnmd[n][P+m] += z * Ynm[n][mabs];
      }
    }
  }
  complex_t z = wavek * radius;
  get_jn(P, z, scalei, jn, 1, jnd);
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      complex_t zh = jn[n];
      complex_t zhn = jnd[n] * wavek;
      complex_t z = zh * zh + zhn * zhn;
      Lnm[n][P+m] = (zh * Lnm[n][P+m] + zhn * Lnmd[n][P+m]) / z;
    }
  }
  rotate(-theta, P, Lnm, Lrot);
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Lnm[n][P+m] = ephi[P-m] * Lrot[n][P+m];
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Li[n][P+m] += Lnm[n][P+m];
    }
  }
}

void L2P(complex_t wavek, real_t scalej, vec3 Xj, complex_t Lj[P+1][2*P+1],
	 vec3 * Xi, int ni, complex_t * pi, cvec3 * Fi,
         real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1], Ynmd[P+1][P+1];
  complex_t ephi[P+1], jn[P+2], jnd[P+2];
  for (int i=0; i<ni; i++) {
    vec3 dX = Xi[i] - Xj;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    real_t ctheta = cos(theta);
    real_t stheta = sin(theta);
    real_t cphi = cos(phi);
    real_t sphi = sin(phi);
    ephi[1] = exp(I * phi);
    for (int n=2; n<=P; n++) {
      ephi[n] = ephi[n-1] * ephi[1];
    }
    real_t rx = stheta * cphi;
    real_t thetax = ctheta * cphi;
    real_t phix = -sphi;
    real_t ry = stheta * sphi;
    real_t thetay = ctheta * sphi;
    real_t phiy = cphi;
    real_t rz = ctheta;
    real_t thetaz = -stheta;
    real_t phiz = 0;
    get_Ynmd(P, ctheta, Ynm, Ynmd, Anm1, Anm2);
    complex_t z = wavek * r;
    get_jn(P, z, scalej, jn, 1, jnd);
    pi[i] += Lj[0][P] * jn[0];
    for (int n=0; n<=P; n++) {
      jnd[n] *= wavek;
    }
    complex_t ur = Lj[0][P] * jnd[0];
    complex_t utheta = 0;
    complex_t uphi = 0;
    for (int n=1; n<=P; n++) {
      pi[i] += Lj[n][P] * jn[n] * Ynm[n][0];
      ur += jnd[n] * Ynm[n][0] * Lj[n][P];
      complex_t jnuse = jn[n+1] * scalej + jn[n-1] / scalej;
      jnuse = wavek * jnuse / (2 * n + 1.0);
      utheta -= Lj[n][P] * jnuse * Ynmd[n][0] * stheta;
      for (int m=1; m<=n; m++) {
	complex_t ztmp1 = jn[n] * Ynm[n][m] * stheta;
	complex_t ztmp2 = Lj[n][P+m] * ephi[m];
	complex_t ztmp3 = Lj[n][P-m] * conj(ephi[m]);
	complex_t ztmpsum = ztmp2 + ztmp3;
	pi[i] += ztmp1 * ztmpsum;
	ur += jnd[n] * Ynm[n][m] * stheta * ztmpsum;
	utheta -= ztmpsum * jnuse * Ynmd[n][m];
	ztmpsum = real_t(m) * I * (ztmp2 - ztmp3);
	uphi += jnuse * Ynm[n][m] * ztmpsum;
      }
    }
    complex_t ux = ur * rx + utheta * thetax + uphi * phix;
    complex_t uy = ur * ry + utheta * thetay + uphi * phiy;
    complex_t uz = ur * rz + utheta * thetaz + uphi * phiz;
    Fi[i][0] -= ux;
    Fi[i][1] -= uy;
    Fi[i][2] -= uz;
  }
}
