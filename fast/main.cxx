#include "dataset.h"
#include "tree.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, bodies2;
  Cells cells;
  Dataset DATA;
  DATA.kernelName = "Laplace";
  SerialFMM FMM;
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
  FMM.printNow = true;
#if BUILD
  for ( int it=32; it<33; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
#endif
  numBodies = int(pow(10,(it+24)/8.0));
  std::cout << "N             : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.random(bodies);
  FMM.startTimer("FMM          ");
  FMM.topdown(bodies,cells);
//  FMM.bottomup(bodies,cells);
#if BUILD
#else
  FMM.approximate(cells);
  FMM.stopTimer("FMM          ",true);
  FMM.eraseTimer("FMM          ");
  FMM.writeTime();
  FMM.resetTimer();

#ifdef DIRECT
  bodies2 = bodies;
  DATA.initTarget(bodies);
  FMM.startTimer("Direct sum   ");
  FMM.direct(bodies,cells);
  FMM.stopTimer("Direct sum   ",true);
  FMM.eraseTimer("Direct sum   ");

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  DATA.printError(diff1,norm1,diff2,norm2);
#endif
#endif
  }
}
