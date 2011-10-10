#include "dataset.h"
#include "tree.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, bodies2;
  Dataset D;
  D.kernelName = "Laplace";
  TreeConstructor T;
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
#if BUILD
  T.printNow = true;
  for ( int it=32; it<33; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
#endif
  numBodies = int(pow(10,(it+24)/8.0));
  std::cout << "N             : " << numBodies << std::endl;
  bodies.resize(numBodies);
  D.random(bodies);
  T.startTimer("FMM          ");
  T.topdown(bodies);
//  T.bottomup(bodies);
#if BUILD
#else
  T.approximate(bodies);
  T.stopTimer("FMM          ",true);
  T.eraseTimer("FMM          ");
  T.writeTime();
  T.resetTimer();

#ifdef DIRECT
  bodies2 = bodies;
  T.startTimer("Direct sum   ");
  T.exact(bodies2);
  T.stopTimer("Direct sum   ",true);
  T.eraseTimer("Direct sum   ");

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  D.printError(diff1,norm1,diff2,norm2);
#endif
#endif
  }
}
