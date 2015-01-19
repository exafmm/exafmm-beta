void evaluate(complex_t wavek, int numBodies, vec3 * Xj, complex_t * qj, complex_t * pi, cvec3 * Fi,
	      complex_t (* Multipole)[P+1][2*P+1], complex_t (* Local)[P+1][2*P+1], int numCells,
	      int numLevels, real_t * scale, real_t R0) {
  int list[189];
  real_t xquad[2*P], wquad[2*P];
  real_t Anm1[P+1][P+1], Anm2[P+1][P+1];
  for (int i=0; i<numBodies; i++) {
    pi[i] = 0;
    Fi[i] = 0;
  }
  getAnm(Anm1,Anm2);
  for (int icell=0; icell<numCells; icell++) {
    for (int n=0; n<=P; n++) {
      for (int m=0; m<=2*P; m++) {
	Multipole[icell][n][m] = 0;
	Local[icell][n][m] = 0;
      }
    }
  }
  for (int level=2; level<=numLevels; level++) {
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      if (cells[icell][8] == 0) break;
      if (cells[icell][6] == 0) {
	int ibegin = cells[icell][7];
	int isize = cells[icell][8];
	P2M(wavek,scale[level],&Xj[ibegin],&qj[ibegin],isize,centers[icell],Multipole[icell],Anm1,Anm2);
      }
    }
  }
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
  complex_t (* Multipole)[P+1][2*P+1] = new complex_t [numCells][P+1][2*P+1]();
  complex_t (* Local)[P+1][2*P+1] = new complex_t [numCells][P+1][2*P+1]();
  evaluate(wavek, numBodies, Xjd, qjd, pid, Fid, Multipole, Local, numCells, numLevels, scale, R0);
  delete[] permutation;
  delete[] scale;
  delete[] Xjd;
  delete[] qjd;
  delete[] pid;
  delete[] Fid;
  delete[] Multipole;
  delete[] Local;
}
