void evaluate(int numBodies, vec3 * Xj, real_t * qj, real_t * pi, vec3 * Fi,
	      complex_t (* Multipole)[NTERM], complex_t (* Local)[NTERM], int numCells,
	      int numLevels) {
  int list[189];
  for (int i=0; i<numBodies; i++) {
    pi[i] = 0;
    Fi[i] = 0;
  }
  for (int icell=0; icell<numCells; icell++) {
    for (int n=0; n<=NTERM; n++) {
      Multipole[icell][n] = 0;
      Local[icell][n] = 0;
    }
  }

  logger::startTimer("P2M");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      if (cells[icell][6] == 0) {
	int ibegin = cells[icell][7];
	int isize = cells[icell][8];
	P2M(&Xj[ibegin], &qj[ibegin], isize, centers[icell], Multipole[icell]);
      }
    }
  }
  logger::stopTimer("P2M");

  logger::startTimer("M2M");
  for (int level=numLevels; level>2; level--) {
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells[icell][6]; ilist++) {
	int jcell = cells[icell][5] + ilist;
	M2M(centers[jcell], Multipole[jcell], centers[icell], Multipole[icell]);
      }
    }
  }
  logger::stopTimer("M2M");

  logger::startTimer("M2L");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      int nlist;
      getList(1, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	M2L(centers[jcell], Multipole[jcell], centers[icell], Local[icell]);
      }
    }
  }
  logger::stopTimer("M2L");

  logger::startTimer("L2L");
  for (int level=3; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level-1]; icell<levelOffset[level]; icell++) {
      for (int ilist=0; ilist<cells[icell][6]; ilist++) {
	int jcell = cells[icell][5]+ilist;
	L2L(centers[icell], Local[icell], centers[jcell], Local[jcell]);
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
        L2P(centers[icell], Local[icell], &Xj[ibegin], isize,
	    &pi[ibegin], &Fi[ibegin]);
      }
    }
  }
  logger::stopTimer("L2P");

  logger::startTimer("P2P");
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    if (cells[icell][6] == 0) {
      P2P(cells[icell], pi, Fi, cells[icell], Xj, qj);
      int nlist;
      getList(0, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	P2P(cells[icell], pi, Fi, cells[jcell], Xj, qj);
      }
    }
  }
  logger::stopTimer("P2P");
}

void fmm(int numBodies, vec3 * Xj, real_t * qj, real_t * pi, vec3 * Fi) {
  int * permutation = new int [numBodies];
  vec3 * Xjd = new vec3 [numBodies]; 
  real_t * qjd = new real_t [numBodies];
  real_t * pid = new real_t [numBodies];
  vec3 * Fid = new vec3 [numBodies];
  levelOffset = new int [maxLevel];
  vec3 X0;
  real_t R0;
  logger::startTimer("Tree");
  getBounds(Xj, numBodies, X0, R0);
  int numCells, numLevels;
  buildTree(Xj, numBodies, numCells, permutation, numLevels, X0, R0);
  for (int i=0; i<numBodies; i++) {
    Xjd[i] = Xj[permutation[i]];
    qjd[i] = qj[permutation[i]];
  }
  logger::stopTimer("Tree");
  complex_t (* Multipole)[NTERM] = new complex_t [numCells][NTERM]();
  complex_t (* Local)[NTERM] = new complex_t [numCells][NTERM]();
  evaluate(numBodies, Xjd, qjd, pid, Fid, Multipole, Local, numCells, numLevels);
  for (int i=0; i<numBodies; i++) {
    pi[permutation[i]] = pid[i];
    Fi[permutation[i]] = Fid[i];
  }
  delete[] permutation;
  delete[] Xjd;
  delete[] qjd;
  delete[] pid;
  delete[] Fid;
  delete[] Multipole;
  delete[] Local;
}
