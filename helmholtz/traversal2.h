#include <omp.h>
#include "kernel.h"

int (* listOffset)[3];
int (* lists)[2];

void resetCellRadius(C_iter C, C_iter C0, real_t R0, int level) {
  C->R = R0 / (1 << level);
  for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {
    resetCellRadius(CC, C0, R0, level+1);
  }
}

void getList(int itype, int icell, int * list, int & numList) {
  int ilast = listOffset[icell][itype];
  numList = 0;
  while (ilast >= 0) {
    if (lists[ilast][1] > 0) {
      list[numList] = lists[ilast][1];
      numList++;
    }
    ilast = lists[ilast][0];
  }
}

void setList(int itype, int icell, int list, int & numLists) {
  lists[numLists][0] = listOffset[icell][itype];
  lists[numLists][1] = list;
  listOffset[icell][itype] = numLists;
  numLists++;
}

void setLists(Cells cells) {
  int numCells = cells.size();
  int childs[216], neighbors[27];
  C_iter C0 = cells.begin();
  for (int i=0; i<numCells; i++) {
    for (int j=0; j<3; j++) {
      listOffset[i][j] = -1;
    }
  }
  int numLists = 0;
  for (int icell=1; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    int iparent = Ci->IPARENT;
    neighbors[0] = iparent;
    int numNeighbors;
    getList(2, iparent, &neighbors[1], numNeighbors);
    numNeighbors++;
    ivec3 iX = getIndex(Ci->ICELL);
    int nchilds = 0;
    for (int i=0; i<numNeighbors; i++) {
      int jparent = neighbors[i];
      C_iter Cj = C0 + jparent;
      for (int j=0; j<Cj->NCHILD; j++) {
	int jcell = Cj->ICHILD+j;
	if (jcell != icell) {
	  childs[nchilds] = jcell;
	  nchilds++;
	}
      }
    }
    for (int i=0; i<nchilds; i++) {
      int jcell = childs[i];
      C_iter Cj = C0 + jcell;
      ivec3 jX = getIndex(Cj->ICELL);
      if (iX[0]-1 <= jX[0] && jX[0] <= iX[0]+1 &&
	  iX[1]-1 <= jX[1] && jX[1] <= iX[1]+1 &&
	  iX[2]-1 <= jX[2] && jX[2] <= iX[2]+1) {
	setList(2, icell, jcell, numLists);
      }	else {
	setList(1, icell, jcell, numLists);
      }
    }
  }
  for (int icell=0; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    if (Ci->ICHILD == 0) {
      int numNeighbors;
      getList(2, icell, neighbors, numNeighbors);
      for (int j=0; j<numNeighbors; j++) {
	int jcell = neighbors[j];
	C_iter Cj = C0 + jcell;
	if (Cj->ICHILD == 0) {
	  setList(0, icell, jcell, numLists);
	}
      }
    }
  }
}

void listBasedTraversal(Cells & cells) {
  int numCells = cells.size();
  C_iter C0 = cells.begin();
  real_t R0 = C0->R;
  vec3 Xperiodic = 0;
  bool mutual = false;
  int list[189];
  listOffset = new int [numCells][3]();
  lists = new int [189*numCells][2]();
  resetCellRadius(C0, C0, R0, 0);
  setLists(cells);

  logger::startTimer("M2L");
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    int nlist;
    getList(1, icell, list, nlist);
    for (int ilist=0; ilist<nlist; ilist++) {
      int jcell = list[ilist];
      C_iter Cj = C0 + jcell;
      kernel::M2L(Ci, Cj, Xperiodic, mutual);
    }
  }
  logger::stopTimer("M2L");

  logger::startTimer("P2P");
#pragma omp parallel for private(list) schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    if (Ci->NCHILD == 0) {
      kernel::P2P(Ci, Ci, Xperiodic, mutual);
      int nlist;
      getList(0, icell, list, nlist);
      for (int ilist=0; ilist<nlist; ilist++) {
	int jcell = list[ilist];
	C_iter Cj = C0 + jcell;
	kernel::P2P(Ci, Cj, Xperiodic, mutual);
      }
    }
  }
  logger::stopTimer("P2P");
  delete[] listOffset;
  delete[] lists;
}
