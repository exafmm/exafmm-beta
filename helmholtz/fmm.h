void fmm(complex_t wavek, int numBodies, vec3 * Xj, complex_t * qj,
	 complex_t * pi, cvec3 * Fi) {
  int * permutation = new int [numBodies];
  real_t * scale = new real_t [maxLevel];
  vec3 * Xjd = new vec3 [numBodies]; 
  complex_t * qjd = new complex_t [numBodies];
  complex_t * pid = new complex_t [numBodies];
  cvec3 * Fid = new cvec3 [numBodies];
  permutation = new int [numBodies];
  levelOffset = new int [maxLevel];
  vec3 X0;
  real_t R0;
  getBounds(Xj,numBodies,X0,R0);

  delete[] permutation;
  delete[] scale;
  delete[] Xjd;
  delete[] qjd;
  delete[] pid;
  delete[] Fid;
}
