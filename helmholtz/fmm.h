void evaluate(int numBodies, vec3 * Xj, complex_t * qj, complex_t * pi, cvec3 * Fi,
	      complex_t (* Multipole)[P*P], complex_t (* Local)[P*P], int numCells,
	      int numLevels, real_t * scale) {
  int list[189];
  for (int i=0; i<numBodies; i++) {
    pi[i] = 0;
    Fi[i] = 0;
  }
  getAnm();
  C_iter C = cells.begin();
  for (int icell=0; icell<numCells; icell++,C++) {
    C->M = 0;
    C->L = 0;
    for (int n=0; n<P; n++) {
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
      C_iter C = cells.begin() + icell;
      if (C->NCHILD == 0) {
	P2M(Multipole[icell], C);
      }
    }
  }
  logger::stopTimer("P2M");

  logger::startTimer("M2M");
  for (int level=numLevels; level>2; level--) {
    nquad = fmax(6, 2 * P);
    legendre();
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells2[icell][6]; ilist++) {
	int jcell = cells2[icell][5] + ilist;
	M2M(scale[level], centers[jcell], Multipole[jcell],
	    scale[level-1], centers[icell], Multipole[icell]);
      }
    }
  }
  logger::stopTimer("M2M");

  logger::startTimer("M2L");
  for (int level=2; level<=numLevels; level++) {
    nquad = fmax(6, P);
    legendre();
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      int nlist;
      getList(1, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	M2L(scale[level], centers[jcell], Multipole[jcell],
	    scale[level], centers[icell], Local[icell]);
      }
    }
  }
  logger::stopTimer("M2L");

  logger::startTimer("L2L");
  for (int level=3; level<=numLevels; level++) {
    nquad = fmax(6, P);
    legendre();
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells2[icell][6]; ilist++) {
	int jcell = cells2[icell][5]+ilist;
	L2L(scale[level-1], centers[icell], Local[icell],
	    scale[level], centers[jcell], Local[jcell]);
      }
    }
  }
  logger::stopTimer("L2L");

  logger::startTimer("L2P");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      if (cells2[icell][6] == 0) {
	int ibegin = cells2[icell][7];
        int isize = cells2[icell][8];
        L2P(scale[level], centers[icell], Local[icell], &Xj[ibegin], isize,
	    &pi[ibegin], &Fi[ibegin]);
      }
    }
  }
  logger::stopTimer("L2P");

  logger::startTimer("P2P");
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    if (cells2[icell][6] == 0) {
      P2P(cells2[icell], pi, Fi, cells2[icell], Xj, qj);
      int nlist;
      getList(0, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	P2P(cells2[icell], pi, Fi, cells2[jcell], Xj, qj);
      }
    }
  }
  logger::stopTimer("P2P");
}

void fmm(int numBodies, vec3 * Xj, complex_t * qj, complex_t * pi, cvec3 * Fi) {
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
    scale[level] = (2 * R0 / (1 << level));
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      cells[icell].R = scale[level];
    }
  }
  Bodies buffer(numBodies);
  for (int i=0; i<numBodies; i++) {
    buffer[i] = bodies[permutation[i]];
    Xjd[i] = Xj[permutation[i]];
    qjd[i] = qj[permutation[i]];
  }
  B_iter B = buffer.begin();
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    C->BODY = B + C->IBODY;
  }
  logger::stopTimer("Tree");
  complex_t (* Multipole)[P*P] = new complex_t [numCells][P*P]();
  complex_t (* Local)[P*P] = new complex_t [numCells][P*P]();
  evaluate(numBodies, Xjd, qjd, pid, Fid, Multipole, Local, numCells, numLevels, scale);
  for (int i=0; i<numBodies; i++) {
    bodies[permutation[i]].TRG = buffer[i].TRG;
    pi[permutation[i]] = pid[i];
    Fi[permutation[i]] = Fid[i];
  }
  delete[] listOffset;
  delete[] lists;
  delete[] levelOffset;
  delete[] centers;
  delete[] permutation;
  delete[] scale;
  delete[] Xjd;
  delete[] qjd;
  delete[] pid;
  delete[] Fid;
  delete[] Multipole;
  delete[] Local;
}
