void P2M(complex_t wavek, real_t scale, vec3 * Xj, complex_t * qj, int nj, vec3 Xi, complex_t Mi[P+1][2*P+1],
	 real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t Ynm[P+1][P+1];
  complex_t imag(0.0,1.0);
  complex_t ephi[P+1], jn[P+1], jnd[P+1], Mnm[P+1][2*P+1];
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = 0;
    }
  }
  for (int i=0; i<nj; i++) {
    vec3 dX = Xj[i] - Xi;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    real_t ctheta = cos(theta);
    real_t stheta = sin(theta);
    ephi[1] = exp(imag * phi);
    for (int n=2; n<=P; n++) {
      ephi[n] = ephi[n-1] * ephi[1];
    }
    get_Ynm(P, ctheta, Ynm, Anm1, Anm2);
    complex_t z = wavek * r;
    get_jn(P, z, scale, jn, 0, jnd);
    for (int n=0; n<=P; n++) {
      jn[n] *= qj[i];
    }
    Mnm[0][P] += jn[0];
    for (int n=1; n<=P; n++) {
      Mnm[n][P] += Ynm[n][0] * jn[n];
      for (int m=1; m<=n; m++) {
	complex_t Ynmjn = Ynm[n][m] * jn[n];
	Mnm[n][P+m] += Ynmjn * conj(ephi[m]);
	Mnm[n][P-m] += Ynmjn * ephi[m];
      }
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mi[n][P+m] += Mnm[n][P+m] * imag * wavek;
    }
  }
}

void M2M(complex_t wavek, real_t scalej, vec3 Xj, complex_t Mj[P+1][2*P+1],
	 real_t scalei, vec3 Xi, complex_t Mi[P+1][2*P+1],
	 real_t radius, real_t xquad[2*P], real_t wquad[2*P], int nquad,
	 real_t Anm1[P+1][P+1], real_t Anm2[P+1][P+1]) {
  real_t ynm[P+1][P+1];
  complex_t imag(0.0,1.0);
  complex_t Mnm[P+1][2*P+1];
  complex_t Mrot[P+1][2*P+1];
  complex_t phitemp[nquad][2*P+1];
  complex_t hn[P+1];
  complex_t ephi[2*P+1];
  vec3 dX = Xi - Xj;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(imag * phi);
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
    get_Ynm(P, cthetaj, ynm, Anm1, Anm2);
    get_hn(P, z, scalej, hn);
    for (int m=-P; m<=P; m++) {
      int mabs = abs(m);
      for (int n=mabs; n<=P; n++) {
	phitemp[l][P+m] += Mrot[n][P+m] * hn[n] * ynm[n][mabs];
      }
    }
  }
  for (int n=0; n<=P; n++) {
    for (int m=-n; m<=n; m++) {
      Mnm[n][P+m] = 0;
    }
  }
  for (int l=0; l<nquad; l++) {
    get_Ynm(P, xquad[l], ynm, Anm1, Anm2);
    for (int m=-P; m<=P; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[l][P+m] * wquad[l] * .5;
      for (int n=mabs; n<=P; n++) {
	Mnm[n][P+m] += z * ynm[n][mabs];
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
  real_t ynm[P+1][P+1], ynmd[P+1][P+1];
  complex_t imag(0.0,1.0);
  complex_t phitemp[nquad][2*Popt+1], phitempn[nquad][2*Popt+1];
  complex_t hn[P+1], hnd[P+1], jn[P+2], jnd[P+2], ephi[2*P+1];
  complex_t Mnm[P+1][2*P+1], Mrot[P+1][2*P+1];
  complex_t Lnm[P+1][2*P+1], Lrot[P+1][2*P+1], Lnmd[P+1][2*P+1];
  vec3 dX = Xi - Xj;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(imag * phi);
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
    get_Ynmd(Popt, cthetaj, ynm, ynmd, Anm1, Anm2);
    get_hnd(Popt, z, scalej, hn, hnd);
    for (int n=0; n<=Popt; n++) {
      hnd[n] *= wavek;
    }
    for (int n=1; n<=Popt; n++) {
      for (int m=1; m<=n; m++) {
	ynm[n][m] *= sthetaj;
      }
    }
    phitemp[l][Popt] = Mrot[0][P] * hn[0];
    phitempn[l][Popt] = Mrot[0][P] * hnd[0] * rn;
    for (int n=1; n<=Popt; n++) {
      phitemp[l][Popt] += Mrot[n][P] * hn[n] * ynm[n][0];
      complex_t ut1 = hnd[n] * rn;
      complex_t ut2 = hn[n] * thetan;
      complex_t ut3 = ut1 * ynm[n][0] - ut2 * ynmd[n][0] * sthetaj;
      phitempn[l][Popt] += ut3 * Mrot[n][P];
      for (int m=1; m<=n; m++) {
	z = hn[n] * ynm[n][m];
	phitemp[l][Popt+m] += Mrot[n][P+m] * z;
	phitemp[l][Popt-m] += Mrot[n][P-m] * z;
	ut3 = ut1 * ynm[n][m] - ut2 * ynmd[n][m];
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
    get_Ynm(Popt, cthetaj, ynm, Anm1, Anm2);
    for (int m=-Popt; m<=Popt; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[l][Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<=Popt; n++) {
	Lnm[n][P+m] += z * ynm[n][mabs];
      }
      z = phitempn[l][Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<=Popt; n++) {
	Lnmd[n][P+m] += z * ynm[n][mabs];
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
