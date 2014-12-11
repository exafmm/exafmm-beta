#include "ewald2.h"
#include "fmm.h"
#include "logger2.h"

#ifdef SAKURA
#define DIM 3
#define NN 26
#define FN 37
#define CN 152
void formInteractionStencil(int *common_stencil, int *far_stencil, int *near_stencil);
#endif

int main() {
  Fmm FMM;
  const int numBodies = 10000;
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int numNeighbors = 1;
  const int numImages = 3;
  const real cycle = 10 * M_PI;
  real potDif = 0, potNrm = 0, accDif = 0, accNrm = 0;
  printf("Height: %d\n", maxLevel);

#ifdef SAKURA
  int common_stencil[DIM*CN] = {0};
  int far_stencil[8*DIM*FN] = {0};
  int near_stencil[8*DIM*NN] = {0};
  formInteractionStencil(common_stencil, far_stencil, near_stencil);
#endif

  logger::verbose = true;
  logger::printTitle("FMM Profiling");

  logger::startTimer("Allocate");
  FMM.allocate(numBodies, maxLevel, numNeighbors, numImages);
  logger::stopTimer("Allocate");

  logger::startTimer("Init bodies");
  FMM.initBodies(cycle);
  logger::stopTimer("Init bodies");
  
  logger::startTimer("Sort bodies");
  FMM.sortBodies();
  logger::stopTimer("Sort bodies");


  logger::startTimer("Fill leafs");
  FMM.fillLeafs();
  logger::stopTimer("Fill leafs");
  
  logger::startTimer("P2M");
  FMM.P2M();
  logger::stopTimer("P2M");

  logger::startTimer("M2M");
  FMM.M2M();
  logger::stopTimer("M2M");

  logger::startTimer("M2L");
#ifndef SAKURA
  FMM.M2L();
#else
  FMM.M2L(common_stencil, far_stencil);
#endif
  logger::stopTimer("M2L");

  logger::startTimer("L2L");
  FMM.L2L();
  logger::stopTimer("L2L");

  logger::startTimer("L2P");
  FMM.L2P();
  logger::stopTimer("L2P");

  logger::startTimer("P2P");
#ifndef SAKURA
  FMM.P2P();
#else
  FMM.P2P(near_stencil);
#endif
  logger::stopTimer("P2P");

  logger::startTimer("Verify");

#if 0
  FMM.direct();
  FMM.verify(100, potDif, potNrm, accDif, accNrm);
#else
  const int ksize = 11;
  const real alpha = 10 / cycle;
  const real sigma = .25 / M_PI;
  const real cutoff = 10;
  Ewald ewald(numBodies, maxLevel, cycle, ksize, alpha, sigma, cutoff);
  ewald.wavePart(FMM.Ibodies2, FMM.Jbodies);
  ewald.realPart(FMM.Ibodies2, FMM.Jbodies, FMM.Leafs);
  FMM.dipoleCorrection();
  FMM.verify(numBodies, potDif, potNrm, accDif, accNrm);
#endif

  logger::stopTimer("Verify");

  logger::startTimer("Deallocate");
  FMM.deallocate();
  logger::stopTimer("Deallocate");

  logger::printTitle("FMM vs. direct");
  logger::printError("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
  logger::printError("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
}
