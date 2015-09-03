void evaluate(complex_t wavek, int numBodies, vec3 * Xj, complex_t * qj, complex_t * pi, cvec3 * Fi,
	      complex_t (* Multipole)[(P+1)*(P+1)], complex_t (* Local)[(P+1)*(P+1)], int numCells,
	      int numLevels, real_t * scale, real_t R0) {
  int list[189];
  real_t xquad[2*P], wquad[2*P];
  for (int i=0; i<numBodies; i++) {
    pi[i] = 0;
    Fi[i] = 0;
  }
  getAnm();
  for (int icell=0; icell<numCells; icell++) {
    for (int n=0; n<=P; n++) {
      for (int m=-n; m<=n; m++) {
	int nm = n * n + n + m;
	Multipole[icell][nm] = 0;
	Local[icell][nm] = 0;
      }
    }
  }

  logger::startTimer("P2M");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      if (cells[icell][6] == 0) {
	int ibegin = cells[icell][7];
	int isize = cells[icell][8];
	P2M(wavek, scale[level], &Xj[ibegin], &qj[ibegin], isize,
	    centers[icell], Multipole[icell]);
      }
    }
  }
  logger::stopTimer("P2M");

  logger::startTimer("M2M");
  for (int level=numLevels; level>2; level--) {
    real_t radius = 2 * R0 / (1 << level) * sqrt(3.0);
    int nquad = fmax(6, 2 * P);
    legendre(nquad, xquad, wquad);
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells[icell][6]; ilist++) {
	int jcell = cells[icell][5] + ilist;
	M2M(wavek, scale[level], centers[jcell], Multipole[jcell],
	    scale[level-1], centers[icell], Multipole[icell],
	    radius, xquad, wquad, nquad);
      }
    }
  }
  logger::stopTimer("M2M");

  logger::startTimer("M2L");
  real_t coef1 = P * 1.65 - 15.5;
  real_t coef2 = P * 0.25 + 3.0;
  for (int level=2; level<=numLevels; level++) {
    real_t diameter = 2 * R0 / (1 << level);
    real_t radius = diameter * sqrt(3.0) * 0.5;
    int nquad = fmax(6, P);
    legendre(nquad, xquad, wquad);
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      int nlist;
      getList(1, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	real_t dx = fabs(cells[jcell][1] - cells[icell][1]);
	real_t dy = fabs(cells[jcell][2] - cells[icell][2]);
	real_t dz = fabs(cells[jcell][3] - cells[icell][3]);
	if (dx > 0) dx -= .5;
	if (dy > 0) dy -= .5;
	if (dz > 0) dz -= .5;
	real_t rr = sqrt(dx * dx + dy * dy + dz * dz);
	int Popt = coef1 / (rr * rr) + coef2;
	M2L(wavek, scale[level], centers[jcell], Multipole[jcell],
	    scale[level], centers[icell], Local[icell],
	    Popt, radius, xquad, wquad, nquad);
      }
    }
  }
  logger::stopTimer("M2L");

  logger::startTimer("L2L");
  for (int level=3; level<=numLevels; level++) {
    real_t radius = 2 * R0 / (1 << level) * sqrt(3.0);
    int nquad = fmax(6, P);
    legendre(nquad, xquad, wquad);
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells[icell][6]; ilist++) {
	int jcell = cells[icell][5]+ilist;
	L2L(wavek, scale[level-1], centers[icell], Local[icell],
	    scale[level], centers[jcell], Local[jcell],
	    radius, xquad, wquad, nquad);
      }
    }
  }
  logger::stopTimer("L2L");

  logger::startTimer("L2P");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      if (cells[icell][6] == 0) {
	int ibegin = cells[icell][7];
        int isize = cells[icell][8];
        L2P(wavek, scale[level], centers[icell], Local[icell], &Xj[ibegin], isize,
	    &pi[ibegin], &Fi[ibegin]);
      }
    }
  }
  logger::stopTimer("L2P");

  logger::startTimer("P2P");
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    if (cells[icell][6] == 0) {
      P2P(cells[icell], pi, Fi, cells[icell], Xj, qj, wavek);
      int nlist;
      getList(0, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	P2P(cells[icell], pi, Fi, cells[jcell], Xj, qj, wavek);
      }
    }
  }
  logger::stopTimer("P2P");
}

void fmm(complex_t wavek, int numBodies, vec3 * Xj, complex_t * qj, complex_t * pi, cvec3 * Fi) {
  int * permutation = new int [numBodies];
  real_t * scale = new real_t [maxLevel];
  vec3 * Xjd = new vec3 [numBodies]; 
  complex_t * qjd = new complex_t [numBodies];
  complex_t * pid = new complex_t [numBodies];
  cvec3 * Fid = new cvec3 [numBodies];
  levelOffset = new int [maxLevel];
  vec3 X0;
  real_t R0;
  logger::startTimer("Tree");
  getBounds(Xj, numBodies, X0, R0);
  int numCells, numLevels;
  buildTree(Xj, numBodies, numCells, permutation, numLevels, X0, R0);
  for (int level=0; level<=numLevels; level++) {
    scale[level] = (2 * R0 / (1 << level)) * abs(wavek);
    if (scale[level]>=1.0) scale[level] = 1.0;
  }
  for (int i=0; i<numBodies; i++) {
    Xjd[i] = Xj[permutation[i]];
    qjd[i] = qj[permutation[i]];
  }
  logger::stopTimer("Tree");
  complex_t (* Multipole)[(P+1)*(P+1)] = new complex_t [numCells][(P+1)*(P+1)]();
  complex_t (* Local)[(P+1)*(P+1)] = new complex_t [numCells][(P+1)*(P+1)]();
  evaluate(wavek, numBodies, Xjd, qjd, pid, Fid, Multipole, Local, numCells, numLevels, scale, R0);
  for (int i=0; i<numBodies; i++) {
    pi[permutation[i]] = pid[i];
    Fi[permutation[i]] = Fid[i];
  }
  delete[] permutation;
  delete[] scale;
  delete[] Xjd;
  delete[] qjd;
  delete[] pid;
  delete[] Fid;
  delete[] Multipole;
  delete[] Local;
}
