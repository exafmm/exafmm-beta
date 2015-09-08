void P2P(C_iter Ci, C_iter Cj) {
  for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {
    complex_t p = 0.0;
    complex_t F[3] = {0.0,0.0,0.0};
    for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {
      vec3 dX = Bi->X - Bj->X;
      real_t R2 = norm(dX);
      if (R2 != 0) {
	real_t R = sqrt(R2);
	complex_t coef1 = Bj->SRC * exp(I * wavek * R) / R;
	complex_t coef2 = (1.0 - I * wavek * R) * coef1 / R2;
	p += coef1;
	F[0] += coef2 * dX[0];
	F[1] += coef2 * dX[1];
	F[2] += coef2 * dX[2];
      }
    }
    Bi->TRG[0] += p;
    Bi->TRG[1] += F[0];
    Bi->TRG[2] += F[1];
    Bi->TRG[3] += F[2];
  }
}

void P2M(C_iter C) {
  real_t Ynm[P*(P+1)/2];
  complex_t ephi[P], jn[P+1], jnd[P+1];
  vecP Mnm = complex_t(0,0);
  real_t scale = C->R;
  real_t kscale = scale * abs(wavek);
  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
    vec3 dX = B->X - C->X;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    real_t ctheta = cos(theta);
    ephi[1] = exp(I * phi);
    for (int n=2; n<P; n++) {
      ephi[n] = ephi[n-1] * ephi[1];
    }
    get_Ynm(P, ctheta, Ynm);
    complex_t z = wavek * r;
    get_jn(P, z, kscale, jn, 0, jnd);
    for (int n=0; n<P; n++) {
      jn[n] *= B->SRC;
    }
    for (int n=0; n<P; n++) {
      int nm = n * n + n;
      int nms = n * (n + 1) / 2;
      Mnm[nm] += Ynm[nms] * jn[n];
      for (int m=1; m<=n; m++) {
	nms = n * (n + 1) / 2 + m;
	int npm = n * n + n + m;
	int nmm = n * n + n - m;
	complex_t Ynmjn = Ynm[nms] * jn[n];
	Mnm[npm] += Ynmjn * conj(ephi[m]);
	Mnm[nmm] += Ynmjn * ephi[m];
      }
    }
  }
  C->M += Mnm * I * wavek;
}

void M2M(C_iter Ci, C_iter C0) {
  real_t Ynm[P*(P+1)/2];
  complex_t phitemp[2*P], hn[P], ephi[2*P];
  vecP Mnm = complex_t(0,0);
  vecP Mrot = complex_t(0,0);
  real_t kscalei = Ci->R * abs(wavek);
  for (C_iter Cj=C0+Ci->ICHILD; Cj!=C0+Ci->ICHILD+Ci->NCHILD; Cj++) {
    real_t kscalej = Cj->R * abs(wavek);
    real_t radius = Cj->R * sqrt(3.0);
    vec3 dX = Ci->X - Cj->X;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    ephi[P+1] = exp(I * phi);
    ephi[P] = 1;
    ephi[P-1] = conj(ephi[P+1]);
    for (int n=2; n<P; n++) {
      ephi[P+n] = ephi[P+n-1] * ephi[P+1];
      ephi[P-n] = conj(ephi[P+n]);
    }
    for (int n=0; n<P; n++) {
      for (int m=-n; m<=n; m++) {
	int nm = n * n + n + m;
	Mnm[nm] = Cj->M[nm] * ephi[P+m];
      }
    }
    rotate(theta, P, Mnm, Mrot);
    for (int n=0; n<P; n++) {
      for (int m=-n; m<=n; m++) {
	int nm = n * n + n + m;
	Mnm[nm] = 0;
      }
    }
    for (int l=0; l<nquad; l++) {
      real_t ctheta = xquad[l];
      real_t stheta = sqrt(1 - ctheta * ctheta);
      real_t rj = (r + radius * ctheta) * (r + radius * ctheta) + (radius * stheta) * (radius * stheta);
      rj = sqrt(rj);
      real_t cthetaj = (r + radius * ctheta) / rj;
      complex_t z = wavek * rj;
      get_Ynm(P, cthetaj, Ynm);
      get_hn(P, z, kscalej, hn);
      for (int m=-P+1; m<P; m++) {
	int mabs = abs(m);
	phitemp[P+m] = 0;
	for (int n=mabs; n<P; n++) {
	  int nm = n * n + n + m;
	  int nms = n * (n + 1) / 2 + mabs;
	  phitemp[P+m] += Mrot[nm] * hn[n] * Ynm[nms];
	}
      }
      get_Ynm(P, xquad[l], Ynm);
      for (int m=-P+1; m<P; m++) {
	int mabs = abs(m);
	complex_t z = phitemp[P+m] * wquad[l] * .5;
	for (int n=mabs; n<P; n++) {
	  int nm = n * n + n + m;
	  int nms = n * (n + 1) / 2 + mabs;
	  Mnm[nm] += z * Ynm[nms];
	}
      }
    }
    complex_t z = wavek * radius;
    get_hn(P, z, kscalei, hn);
    for (int n=0; n<P; n++) {
      for (int m=-n; m<=n; m++) {
	int nm = n * n + n + m;
	Mnm[nm] /= hn[n];
      }
    }
    rotate(-theta, P, Mnm, Mrot);
    for (int n=0; n<P; n++) {
      for (int m=-n; m<=n; m++) {
	int nm = n * n + n + m;
	Mnm[nm] = ephi[P-m] * Mrot[nm];
      }
    }
    Ci->M += Mnm;
  }
}

void M2L(C_iter Ci, C_iter Cj, vec3 Xperiodic, bool mutual) {
  assert(mutual == false);
  real_t Ynm[P*(P+1)/2], Ynmd[P*(P+1)/2];
  complex_t phitemp[2*P], phitempn[2*P];
  complex_t hn[P], hnd[P], jn[P+1], jnd[P+1], ephi[2*P];
  vecP Lnm = complex_t(0,0);
  vecP Lnmd = complex_t(0,0);
  vecP Mnm, Mrot, Lrot;
  real_t kscalej = Cj->R * abs(wavek);
  real_t kscalei = Ci->R * abs(wavek);
  real_t radius = Cj->R * sqrt(3.0) * .5;
  vec3 dX = Ci->X - Cj->X - Xperiodic;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  dX /= Cj->R;
  if (fabs(dX[0]) > EPS) dX[0] = fabs(dX[0]) - .5;
  if (fabs(dX[1]) > EPS) dX[1] = fabs(dX[1]) - .5;
  if (fabs(dX[2]) > EPS) dX[2] = fabs(dX[2]) - .5;
  real_t rr = sqrt(norm(dX));
  real_t coef1 = P * 1.65 - 15.5;
  real_t coef2 = P * 0.25 + 3.0;
  int Popt = coef1 / (rr * rr) + coef2;  
  ephi[P+1] = exp(I * phi);
  ephi[P] = 1;
  ephi[P-1] = conj(ephi[P+1]);
  for (int n=2; n<P; n++) {
    ephi[P+n] = ephi[P+n-1] * ephi[P+1];
    ephi[P-n] = conj(ephi[P+n]);
  }
  for (int n=0; n<Popt; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Mnm[nm] = Cj->M[nm] * ephi[P+m];
    }
  }
  rotate(theta, Popt, Mnm, Mrot);
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
    get_Ynmd(Popt, cthetaj, Ynm, Ynmd);
    get_hnd(Popt, z, kscalej, hn, hnd);
    for (int n=0; n<Popt; n++) {
      hnd[n] *= wavek;
    }
    for (int n=1; n<Popt; n++) {
      for (int m=1; m<=n; m++) {
	int nms = n * (n + 1) / 2 + m;
	Ynm[nms] *= sthetaj;
      }
    }
    for (int m=-Popt+1; m<Popt; m++) {
      phitemp[Popt+m] = 0;
      phitempn[Popt+m] = 0;
    }
    phitemp[Popt] = Mrot[0] * hn[0];
    phitempn[Popt] = Mrot[0] * hnd[0] * rn;
    for (int n=1; n<Popt; n++) {
      int nm = n * n + n;
      int nms = n * (n + 1) / 2;
      phitemp[Popt] += Mrot[nm] * hn[n] * Ynm[nms];
      complex_t ut1 = hnd[n] * rn;
      complex_t ut2 = hn[n] * thetan;
      complex_t ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms] * sthetaj;
      phitempn[Popt] += ut3 * Mrot[nm];
      for (int m=1; m<=n; m++) {
	nms = n * (n + 1) / 2 + m;
	int npm = n * n + n + m;
	int nmm = n * n + n - m;
	z = hn[n] * Ynm[nms];
	phitemp[Popt+m] += Mrot[npm] * z;
	phitemp[Popt-m] += Mrot[nmm] * z;
	ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms];
	phitempn[Popt+m] += ut3 * Mrot[npm];
	phitempn[Popt-m] += ut3 * Mrot[nmm];
      }
    }
    get_Ynm(Popt, xquad[l], Ynm);
    for (int m=-Popt+1; m<Popt; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<Popt; n++) {
	int nm = n * n + n + m;
	int nms = n * (n + 1) / 2 + mabs;
	Lnm[nm] += z * Ynm[nms];
      }
      z = phitempn[Popt+m] * wquad[l] * .5;
      for (int n=mabs; n<Popt; n++) {
	int nm = n * n + n + m;
	int nms = n * (n + 1) / 2 + mabs;
	Lnmd[nm] += z * Ynm[nms];
      }
    }
  }
  complex_t z = wavek * radius;
  get_jn(Popt, z, kscalei, jn, 1, jnd);
  for (int n=0; n<Popt; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      complex_t zh = jn[n];
      complex_t zhn = jnd[n] * wavek;
      complex_t z = zh * zh + zhn * zhn;
      Lnm[nm] = (zh * Lnm[nm] + zhn * Lnmd[nm]) / z;
    }
  }
  rotate(-theta, Popt, Lnm, Lrot);
  for (int n=0; n<Popt; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Lnm[nm] = ephi[P-m] * Lrot[nm];
    }
  }
  Ci->L += Lnm;
}

void L2L(C_iter Ci, C_iter C0) {
  real_t Ynm[P*(P+1)/2], Ynmd[P*(P+1)/2];
  complex_t phitemp[2*P], phitempn[2*P];
  complex_t jn[P+1], jnd[P+1], ephi[2*P];
  vecP Lnm, Lnmd, Lrot;
  real_t kscalei = Ci->R * abs(wavek);
  C_iter Cj = C0 + Ci->IPARENT;
  real_t kscalej = Cj->R * abs(wavek);
  real_t radius = Cj->R * sqrt(3.0) * .5;
  vec3 dX = Ci->X - Cj->X;
  real_t r, theta, phi;
  cart2sph(dX, r, theta, phi);
  ephi[P+1] = exp(I * phi);
  ephi[P] = 1;
  ephi[P-1] = conj(ephi[P+1]);
  for (int n=2; n<P; n++) {
    ephi[P+n] = ephi[P+n-1] * ephi[P+1];
    ephi[P-n] = conj(ephi[P+n]);
  }
  for (int n=0; n<P; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Lnm[nm] = Cj->L[nm] * ephi[P+m];
    }
  }
  rotate(theta, P, Lnm, Lrot);
  for (int n=0; n<P; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Lnm[nm] = 0;
      Lnmd[nm] = 0;
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
    get_Ynmd(P, cthetaj, Ynm, Ynmd);
    get_jn(P, z, kscalej, jn, 1, jnd);
    for (int n=0; n<P; n++) {
      jnd[n] *= wavek;
    }
    for (int n=1; n<P; n++) {
      for (int m=1; m<=n; m++) {
	int nms = n * (n + 1) / 2 + m;
	Ynm[nms] *= sthetaj;
      }
    }
    for (int m=-P+1; m<P; m++) {
      phitemp[P+m] = 0;
      phitempn[P+m] = 0;
    }
    phitemp[P] = Lrot[0] * jn[0];
    phitempn[P] = Lrot[0] * jnd[0] * rn;
    for (int n=1; n<P; n++) {
      int nm = n * n + n;
      int nms = n * (n + 1) / 2;
      phitemp[P] += Lrot[nm] * jn[n] * Ynm[nms];
      complex_t ut1 = jnd[n] * rn;
      complex_t ut2 = jn[n] * thetan;
      complex_t ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms] * sthetaj;
      phitempn[P] += ut3 * Lrot[nm];
      for (int m=1; m<=n; m++) {
	nms = n * (n + 1) / 2 + m;
	int npm = n * n + n + m;
	int nmm = n * n + n - m;
	z = jn[n] * Ynm[nms];
	phitemp[P+m] += Lrot[npm] * z;
	phitemp[P-m] += Lrot[nmm] * z;
	ut3 = ut1 * Ynm[nms] - ut2 * Ynmd[nms];
	phitempn[P+m] += ut3 * Lrot[npm];
	phitempn[P-m] += ut3 * Lrot[nmm];
      }
    }
    get_Ynm(P, xquad[l], Ynm);
    for (int m=-P+1; m<P; m++) {
      int mabs = abs(m);
      complex_t z = phitemp[P+m] * wquad[l] * .5;
      for (int n=mabs; n<P; n++) {
	int nm = n * n + n + m;
	int nms = n * (n + 1) / 2 + mabs;
        Lnm[nm] += z * Ynm[nms];
      }
      z = phitempn[P+m] * wquad[l] * .5;
      for (int n=mabs; n<P; n++) {
	int nm = n * n + n + m;
	int nms = n * (n + 1) / 2 + mabs;
        Lnmd[nm] += z * Ynm[nms];
      }
    }
  }
  complex_t z = wavek * radius;
  get_jn(P, z, kscalei, jn, 1, jnd);
  for (int n=0; n<P; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      complex_t zh = jn[n];
      complex_t zhn = jnd[n] * wavek;
      complex_t z = zh * zh + zhn * zhn;
      Lnm[nm] = (zh * Lnm[nm] + zhn * Lnmd[nm]) / z;
    }
  }
  rotate(-theta, P, Lnm, Lrot);
  for (int n=0; n<P; n++) {
    for (int m=-n; m<=n; m++) {
      int nm = n * n + n + m;
      Lnm[nm] = ephi[P-m] * Lrot[nm];
    }
  }
  Ci->L += Lnm;
}

void L2P(C_iter C) {
  real_t Ynm[P*(P+1)/2], Ynmd[P*(P+1)/2];
  complex_t ephi[P], jn[P+1], jnd[P+1];
  real_t kscale = C->R * abs(wavek);
  for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
    vecP Lj = C->L;
    cvec4 TRG = complex_t(0,0);
    vec3 dX = B->X - C->X;
    real_t r, theta, phi;
    cart2sph(dX, r, theta, phi);
    real_t ctheta = cos(theta);
    real_t stheta = sin(theta);
    real_t cphi = cos(phi);
    real_t sphi = sin(phi);
    ephi[1] = exp(I * phi);
    for (int n=2; n<P; n++) {
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
    get_Ynmd(P, ctheta, Ynm, Ynmd);
    complex_t z = wavek * r;
    get_jn(P, z, kscale, jn, 1, jnd);
    TRG[0] += Lj[0] * jn[0];
    for (int n=0; n<P; n++) {
      jnd[n] *= wavek;
    }
    complex_t ur = Lj[0] * jnd[0];
    complex_t utheta = 0;
    complex_t uphi = 0;
    for (int n=1; n<P; n++) {
      int nm = n * n + n;
      int nms = n * (n + 1) / 2;
      TRG[0] += Lj[nm] * jn[n] * Ynm[nms];
      ur += jnd[n] * Ynm[nms] * Lj[nm];
      complex_t jnuse = jn[n+1] * kscale + jn[n-1] / kscale;
      jnuse = wavek * jnuse / (2 * n + 1.0);
      utheta -= Lj[nm] * jnuse * Ynmd[nms] * stheta;
      for (int m=1; m<=n; m++) {
	int npm = n * n + n + m;
	int nmm = n * n + n - m;
	int nms = n * (n + 1) / 2 + m;
	complex_t ztmp1 = jn[n] * Ynm[nms] * stheta;
	complex_t ztmp2 = Lj[npm] * ephi[m];
	complex_t ztmp3 = Lj[nmm] * conj(ephi[m]);
	complex_t ztmpsum = ztmp2 + ztmp3;
	TRG[0] += ztmp1 * ztmpsum;
	ur += jnd[n] * Ynm[nms] * stheta * ztmpsum;
	utheta -= ztmpsum * jnuse * Ynmd[nms];
	ztmpsum = real_t(m) * I * (ztmp2 - ztmp3);
	uphi += jnuse * Ynm[nms] * ztmpsum;
      }
    }
    complex_t ux = ur * rx + utheta * thetax + uphi * phix;
    complex_t uy = ur * ry + utheta * thetay + uphi * phiy;
    complex_t uz = ur * rz + utheta * thetaz + uphi * phiz;
    TRG[1] -= ux;
    TRG[2] -= uy;
    TRG[3] -= uz;
    B->TRG += TRG;
  }
}
