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

void growTree(vec3 * Xj, int numBodies, int (* nodes)[10], int & numCells,
	      int * permutation, int & numLevels, vec3 X0, real_t R0) {
  int nbody8[8];
  int * iwork = new int [numBodies];
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
  for (int level=0; level<198; level++) {
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
  delete[] iwork;
}

void getList(int itype, int icell, int * list, int & nlist) {
  int ilast = listOffset[icell][itype];
  nlist = 0;
  while (ilast >= 0) {
    if (lists[ilast][1] > 0) {
      list[nlist] = lists[ilast][1];
      nlist++;
    }
    ilast = lists[ilast][0];
  }
}

void setList(int itype, int icell, int list) {
  lists[numele][0] = listOffset[icell][itype];
  lists[numele][1] = list;
  listOffset[icell][itype] = numele;
  numele++;
}

void setLists(int numCells) {
  int childs[216], neighbors[27];
  for (int i=0; i<numCells; i++) {
    for (int j=0; j<3; j++) {
      listOffset[i][j] = -1;
    }
  }
  for (int icell=1; icell<numCells; icell++) {
    int iparent = nodes[icell][4];
    neighbors[0] = iparent;
    int numNeighbors;
    getList(2,iparent,&neighbors[1],numNeighbors);
    numNeighbors++;
    int nchilds = 0;
    for (int i=0; i<numNeighbors; i++) {
      int jparent = neighbors[i];
      for (int j=0; j<nodes[jparent][6]; j++) {
	int jcell = nodes[jparent][5]+j;
	if (jcell != icell) {
	  childs[nchilds] = jcell;
	  nchilds++;
	}
      }
    }
    for (int i=0; i<nchilds; i++) {
      int jcell = childs[i];
      if (nodes[icell][1]-1 <= nodes[jcell][1] && nodes[jcell][1] <= nodes[icell][1]+1 &&
	  nodes[icell][2]-1 <= nodes[jcell][2] && nodes[jcell][2] <= nodes[icell][2]+1 &&
	  nodes[icell][3]-1 <= nodes[jcell][3] && nodes[jcell][3] <= nodes[icell][3]+1) {
	setList(2,icell,jcell);
      }	else {
	setList(1,icell,jcell);
      }
    }
  }
  for (int icell=0; icell<numCells; icell++) {
    if (nodes[icell][5] == 0) {
      int numNeighbors;
      getList(2,icell,neighbors,numNeighbors);
      for (int j=0; j<numNeighbors; j++) {
	int jcell = neighbors[j];
	if (nodes[jcell][5] == 0) {
	  setList(0,icell,jcell);
	}
      }
    }
  }
}

void buildTree(vec3 * Xj, int numBodies, int & numCells, int * permutation,
	       int & numLevels, vec3 X0, real_t R0) {
  nodes = new int [numBodies][10]();
  growTree(Xj, numBodies, nodes, numCells, permutation, numLevels, X0, R0);
  listOffset = new int [numCells][3]();
  lists = new int [189*numCells][2]();
  cells.resize(numCells);
  cells2 = new int [numCells][10]();
  centers = new vec3 [numCells];
  C_iter C = cells.begin();
  for (int i=0; i<numCells; i++,C++) {
    C->LEVEL   = nodes[i][0];
    C->IPARENT = nodes[i][4];
    C->ICHILD  = nodes[i][5];
    C->NCHILD  = nodes[i][6];
    C->IBODY   = nodes[i][7];
    C->NBODY   = nodes[i][8];
    for (int j=0; j<10; j++) {
      cells2[i][j] = nodes[i][j];
    }
    real_t R = R0 / (1 << nodes[i][0]);
    for (int d=0; d<3; d++) {
      C->X[d] = X0[d] - R0 + nodes[i][d+1] * R * 2 + R;
      centers[i][d] = X0[d] - R0 + cells2[i][d+1] * R * 2 + R;
    }
  }
  setLists(numCells);
  delete[] nodes;
}

