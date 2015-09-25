#include "kernel.h"

void upwardPass(Cells & cells) {
  int numCells = cells.size();
  C_iter C0 = cells.begin();
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    C->M = 0;
    C->L = 0;
  }

  logger::startTimer("P2M");
#pragma omp parallel for schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    C_iter C = C0 + icell;
    if (C->NCHILD == 0) {
      kernel::P2M(C);
    }
  }
  logger::stopTimer("P2M");

  logger::startTimer("M2M");
#pragma omp parallel for schedule(dynamic)
  for (int icell=numCells-1; icell>=0; icell--) {
    C_iter Ci = C0 + icell;
    kernel::M2M(Ci, C0);
  }
  logger::stopTimer("M2M");
}

void downwardPass(Cells & cells) {
  int numCells = cells.size();
  C_iter C0 = cells.begin();

  logger::startTimer("L2L");
#pragma omp parallel for schedule(dynamic)
  for (int icell=1; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    kernel::L2L(Ci, C0);
  }
  logger::stopTimer("L2L");

  logger::startTimer("L2P");
#pragma omp parallel for schedule(dynamic)
  for (int icell=0; icell<numCells; icell++) {
    C_iter Ci = C0 + icell;
    if (Ci->NCHILD == 0) {
      kernel::L2P(Ci);
    }
  }
  logger::stopTimer("L2P");
}
