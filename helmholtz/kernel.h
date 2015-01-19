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
  complex_t fhs[P+1];
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
}
