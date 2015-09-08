void evaluate(int numCells, int numLevels) {
  vec3 Xperiodic = 0;
  bool mutual = false;
  int list[189];
  getAnm();
  C_iter C0 = cells.begin();
  C_iter C = C0;
  for (int icell=0; icell<numCells; icell++,C++) {
    C->M = 0;
    C->L = 0;
  }

  logger::startTimer("P2M");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      C_iter C = C0 + icell;
      if (C->NCHILD == 0) {
	kernel::P2M(C);
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
      C_iter Ci = C0 + icell;
      kernel::M2M(Ci, C0);
    }
  }
  logger::stopTimer("M2M");

  logger::startTimer("M2L");
  for (int level=2; level<=numLevels; level++) {
    nquad = fmax(6, P);
    legendre();
#pragma omp parallel for private(list) schedule(dynamic)
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      C_iter Ci = C0 + icell;
      int nlist;
      getList(1, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	C_iter Cj = C0 + jcell;
	kernel::M2L(Ci, Cj, Xperiodic, mutual);
      }
    }
  }
  logger::stopTimer("M2L");

  logger::startTimer("L2L");
  for (int level=3; level<=numLevels; level++) {
    nquad = fmax(6, P);
    legendre();
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      C_iter Ci = C0 + icell;
      kernel::L2L(Ci, C0);
    }
  }
  logger::stopTimer("L2L");

  logger::startTimer("L2P");
  for (int level=2; level<=numLevels; level++) {
#pragma omp parallel for
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      C_iter Ci = C0 + icell;
      if (Ci->NCHILD == 0) {
	kernel::L2P(Ci);
      }
    }
  }
  logger::stopTimer("L2P");

  logger::startTimer("P2P");
  real_t eps2 = 0;
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    if (Ci->NCHILD == 0) {
      kernel::P2P(Ci, Ci, eps2, Xperiodic, mutual);
      int nlist;
      getList(0, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	C_iter Cj = C0 + jcell;
	kernel::P2P(Ci, Cj, eps2, Xperiodic, mutual);
      }
    }
  }
  logger::stopTimer("P2P");
}

void fmm(int numBodies, vec3 * Xj) {
  int * permutation = new int [numBodies];
  levelOffset = new int [maxLevel];
  vec3 X0;
  real_t R0;
  logger::startTimer("Tree");
  getBounds(Xj, numBodies, X0, R0);
  int numCells, numLevels;
  buildTree(Xj, numBodies, numCells, permutation, numLevels, X0, R0);
  for (int level=0; level<=numLevels; level++) {
    real_t scale = (2 * R0 / (1 << level));
    for (int icell=levelOffset[level]; icell<levelOffset[level+1]; icell++) {
      cells[icell].R = scale;
    }
  }
  Bodies buffer(numBodies);
  for (int i=0; i<numBodies; i++) {
    buffer[i] = bodies[permutation[i]];
  }
  B_iter B = buffer.begin();
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    C->BODY = B + C->IBODY;
  }
  logger::stopTimer("Tree");
  evaluate(numCells, numLevels);
  for (int i=0; i<numBodies; i++) {
    bodies[permutation[i]].TRG = buffer[i].TRG;
  }
  delete[] listOffset;
  delete[] lists;
  delete[] levelOffset;
  delete[] permutation;
}
