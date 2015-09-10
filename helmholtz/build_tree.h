void getBounds(vec3 * Xj, int numBodies, vec3 & X0, real_t & R0) {
  vec3 Xmin = Xj[1];
  vec3 Xmax = Xj[1];
  for (int i=0; i<numBodies; i++) {
    Xmin = min(Xj[i],Xmin);
    Xmax = max(Xj[i],Xmax);
  }
  real_t diameter = 0;
  for (int d=0; d<3; d++) {
    real_t dX = Xmax[d] - Xmin[d];
    diameter = dX > diameter ? dX : diameter;
    X0[d] = (Xmax[d] + Xmin[d]) * 0.5;
  }
  R0 = diameter * 0.5;
}

void reorder(vec3 X0, real_t R0, int level, int * iX, vec3 * Xj,
	     int * permutation, int n, int * iwork, int * nbody) {
  int offset[9];
  vec3 X;
  real_t R = R0 / (1 << level);
  for (int d=0; d<3; d++) {
    X[d] = X0[d] - R0 + iX[d] * R * 2 + R;
  }
  for (int i=0; i<8; i++) nbody[i] = 0;
  for (int i=0; i<n; i++) {
    int j = permutation[i];
    int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);
    nbody[octant]++;
  }
  offset[0] = 0;
  for (int i=0; i<8; i++) {
    offset[i+1] = offset[i] + nbody[i];
    nbody[i] = 0;
  }
  for (int i=0; i<n; i++) {
    int j = permutation[i];
    int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);
    iwork[offset[octant]+nbody[octant]] = permutation[i];
    nbody[octant]++;
  }
  for (int i=0; i<n; i++) {
    permutation[i] = iwork[i];
  }
}

uint64_t getKey(ivec3 iX, int level) {
  uint64_t index = ((1 << 3 * level) - 1) / 7;
  for (int l=0; l<level; l++) {
    for (int d=0; d<3; d++) {
      index += (iX[d] & 1) << (3 * l + d);
      iX[d] >>= 1;
    }
  }
  return index;
}

int getLevel(uint64_t key) {
  int level = -1;
  while( int(key) >= 0 ) {
    level++;
    key -= 1 << 3*level;
  }
  return level;
}

ivec3 getIndex(uint64_t key) {
  int level = -1;
  while( int(key) >= 0 ) {
    level++;
    key -= 1 << 3*level;
  }
  key += 1 << 3*level;
  level = 0;
  ivec3 iX = 0;
  int d = 0;
  while( key > 0 ) {
    iX[d] += (key % 2) * (1 << level);
    key >>= 1;
    d = (d+1) % 3;
    if( d == 0 ) level++;
  }
  return iX;
}

void growTree(vec3 * Xj, int numBodies, int (* nodes)[10], int & numCells,
	      int * permutation, int & numLevels, vec3 X0, real_t R0) {
  const int maxLevel = 30;
  int nbody8[8];
  int * iwork = new int [numBodies];
  int * levelOffset = new int [maxLevel];
  nodes[0][0] = 0;
  nodes[0][1] = 0;
  nodes[0][2] = 0;
  nodes[0][3] = 0;
  nodes[0][4] = 0;
  nodes[0][5] = 0;
  nodes[0][6] = 0;
  nodes[0][7] = 0;
  nodes[0][8] = numBodies;
  levelOffset[0] = 0;
  levelOffset[1] = 1;
  for (int i=0; i<numBodies; i++) {
    permutation[i] = i;
  }
  int ncrit = 1000;
  if (P < 40) ncrit = 500;
  if (P < 30) ncrit = 200;
  if (P < 20) ncrit = 100;
  numCells = 1;
  numLevels = 0;
  for (int level=0; level<maxLevel; level++) {
    for (int iparent=levelOffset[level]; iparent<levelOffset[level+1]; iparent++) {
      int nbody = nodes[iparent][8];
      if (nbody > ncrit) {
	int ibody = nodes[iparent][7];
	reorder(X0, R0, level, &nodes[iparent][1], Xj, &permutation[ibody], nbody, iwork, nbody8);
	int nchild = 0;
	int offset = ibody;
	nodes[iparent][5] = numCells;
	for (int i=0; i<8; i++) {
	  nodes[numCells][0] = level + 1;
	  nodes[numCells][1] = nodes[iparent][1] * 2 + i % 2;
	  nodes[numCells][2] = nodes[iparent][2] * 2 + (i / 2) % 2;
	  nodes[numCells][3] = nodes[iparent][3] * 2 + i / 4;
	  nodes[numCells][4] = iparent;
	  nodes[numCells][5] = 0;
	  nodes[numCells][6] = 0;
	  nodes[numCells][7] = offset;
	  nodes[numCells][8] = nbody8[i];
	  nchild++;
	  offset += nbody8[i];
	  numCells++;
	  numLevels=level+1;
	}
	nodes[iparent][6] = nchild;
      }
    }
    levelOffset[level+2] = numCells;
    if (levelOffset[level+1] == levelOffset[level+2]) break;
  }
  delete[] levelOffset;
  delete[] iwork;
}

Cells buildTree(vec3 * Xj, int numBodies, int & numCells, int * permutation,
	       int & numLevels, vec3 X0, real_t R0) {
  int (* nodes)[10] = new int [numBodies][10]();
  growTree(Xj, numBodies, nodes, numCells, permutation, numLevels, X0, R0);
  Cells cells(numCells);
  C_iter C = cells.begin();
  ivec3 iX;
  for (int i=0; i<numCells; i++,C++) {
    int level = nodes[i][0];
    iX[0]      = nodes[i][1];
    iX[1]      = nodes[i][2];
    iX[2]      = nodes[i][3];
    C->ICELL   = getKey(iX, level);
    C->IPARENT = nodes[i][4];
    C->ICHILD  = nodes[i][5];
    C->NCHILD  = nodes[i][6];
    C->IBODY   = nodes[i][7];
    C->NBODY   = nodes[i][8];
    real_t R = R0 / (1 << level);
    C->R = 2 * R;
    for (int d=0; d<3; d++) {
      C->X[d] = X0[d] - R0 + iX[d] * R * 2 + R;
    }
  }
  delete[] nodes;
  return cells;
}

