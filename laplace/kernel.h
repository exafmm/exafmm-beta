void P2P(int * icell, real_t * pi, vec3 * Fi, int * jcell, vec3 * Xj, real_t * qj) {
  for (int i=icell[7]; i<icell[7]+icell[8]; i++) {
    for (int j=jcell[7]; j<jcell[7]+jcell[8]; j++) {
      vec3 dX = Xj[i] - Xj[j];
      real_t R2 = norm(dX);
      real_t invR2 = 1.0 / R2;
      if( R2 == 0 ) invR2 = 0;
      real_t invR = qj[j] * sqrt(invR2);
      real_t invR3 = invR2 * invR;
      pi[i] += invR;
      Fi[i][0] -= dX[0] * invR3;
      Fi[i][1] -= dX[1] * invR3;
      Fi[i][2] -= dX[2] * invR3;
    }
  }
}

void P2M(vec3 * Xj, real_t * qj, int nj, vec3 Xi, complex_t Mi[NTERM]) {
  complex_t Ynm[P*P], YnmTheta[P*P], Mnm[NTERM];
  for (int n=0; n<NTERM; n++) {
    Mnm[n] = 0;
  }
  for (int i=0; i<nj; i++) {
    vec3 dX = Xj[i] - Xi;
    real_t rho, alpha, beta;
    cart2sph(rho, alpha, beta, dX);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int n=0; n<P; n++) {
      for (int m=0; m<=n; m++) {
	int nm  = n * n + n - m;
	int nms = n * (n + 1) / 2 + m;
	Mnm[nms] += qj[i] * Ynm[nm];
      }
    }
  }
  for (int n=0; n<NTERM; n++) {
    Mi[n] += Mnm[n];
  }
}

void M2M(vec3 Xj, complex_t Mj[NTERM], vec3 Xi, complex_t Mi[NTERM]) {
  complex_t Ynm[P*P], YnmTheta[P*P];
  vec3 dX = Xi - Xj;
  real_t rho, alpha, beta;
  cart2sph(rho, alpha, beta, dX);
  evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
  for (int j=0; j<P; j++) {
    for (int k=0; k<=j; k++) {
      int jks = j * (j + 1) / 2 + k;
      complex_t M = 0;
      for (int n=0; n<=j; n++) {
	for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
	  int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
	  int nm    = n * n + n - m;
	  M += Mj[jnkms] * Ynm[nm] * real_t(IPOW2N(m) * ODDEVEN(n));
	}
	for (int m=k; m<=std::min(n,j+k-n); m++) {
	  int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
	  int nm    = n * n + n - m;
	  M += std::conj(Mj[jnkms]) * Ynm[nm] * real_t(ODDEVEN(k+n+m));
	}
      }
      Mi[jks] += M;
    }
  }
}

void M2L(vec3 Xj, complex_t Mj[NTERM], vec3 Xi, complex_t Li[NTERM]) {
  complex_t Ynm[P*P];
  vec3 dX = Xi - Xj;
  real_t rho, alpha, beta;
  cart2sph(rho, alpha, beta, dX);
  evalLocal(rho, alpha, beta, Ynm);
  for (int j=0; j<P; j++) {
    real_t Cnm = ODDEVEN(j);
    for (int k=0; k<=j; k++) {
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for (int n=0; n<P-j; n++) {
	for (int m=-n; m<0; m++) {
	  int nms  = n * (n + 1) / 2 - m;
	  int jnkm = (j + n) * (j + n) + j + n + m - k;
	  L += std::conj(Mj[nms]) * Cnm * Ynm[jnkm];
	}
	for (int m=0; m<=n; m++) {
	  int nms  = n * (n + 1) / 2 + m;
	  int jnkm = (j + n) * (j + n) + j + n + m - k;
	  real_t Cnm2 = Cnm * ODDEVEN((k-m)*(k<m)+m);
	  L += Mj[nms] * Cnm2 * Ynm[jnkm];
	}
      }
      Li[jks] += L;
    }
  }
}

void L2L(vec3 Xj, complex_t Lj[NTERM], vec3 Xi, complex_t Li[NTERM]) {
  complex_t Ynm[P*P], YnmTheta[P*P];
  vec3 dX = Xi - Xj;
  real_t rho, alpha, beta;
  cart2sph(rho, alpha, beta, dX);
  evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
  for (int j=0; j<P; j++) {
    for (int k=0; k<=j; k++) {
      int jks = j * (j + 1) / 2 + k;
      complex_t L = 0;
      for (int n=j; n<P; n++) {
	for (int m=j+k-n; m<0; m++) {
	  int jnkm = (n - j) * (n - j) + n - j + m - k;
	  int nms  = n * (n + 1) / 2 - m;
	  L += std::conj(Lj[nms]) * Ynm[jnkm] * real_t(ODDEVEN(k));
	}
	for (int m=0; m<=n; m++) {
	  if( n-j >= abs(m-k) ) {
	    int jnkm = (n - j) * (n - j) + n - j + m - k;
	    int nms  = n * (n + 1) / 2 + m;
	    L += Lj[nms] * Ynm[jnkm] * real_t(ODDEVEN((m-k)*(m<k)));
	  }
	}
      }
      Li[jks] += L;
    }
  }
}

void L2P(vec3 Xj, complex_t Lj[NTERM], vec3 * Xi, int ni, real_t * pi, vec3 * Fi) {
  complex_t Ynm[P*P], YnmTheta[P*P];
  for (int i=0; i<ni; i++) {
    vec3 dX = Xi[i] - Xj;
    vec3 spherical = 0;
    vec3 cartesian = 0;
    real_t r, theta, phi;
    cart2sph(r, theta, phi, dX);
    evalMultipole(r, theta, phi, Ynm, YnmTheta);
    for (int n=0; n<P; n++) {
      int nm  = n * n + n;
      int nms = n * (n + 1) / 2;
      pi[i] += std::real(Lj[nms] * Ynm[nm]);
      spherical[0] += std::real(Lj[nms] * Ynm[nm]) / r * n;
      spherical[1] += std::real(Lj[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; m++) {
	nm  = n * n + n + m;
	nms = n * (n + 1) / 2 + m;
	pi[i] += 2 * std::real(Lj[nms] * Ynm[nm]);
	spherical[0] += 2 * std::real(Lj[nms] * Ynm[nm]) / r * n;
	spherical[1] += 2 * std::real(Lj[nms] * YnmTheta[nm]);
	spherical[2] += 2 * std::real(Lj[nms] * Ynm[nm] * I) * m;
      }
    }
    sph2cart(r, theta, phi, spherical, cartesian);
    Fi[i] += cartesian;
  }
}
