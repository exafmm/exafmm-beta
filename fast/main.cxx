#include "dataset.h"
#include "parallelfmm.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, bodies2;
  Cells cells;
  Dataset DATA;
  ParallelFMM FMM;
#if HYBRID
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
  numBodies = int(pow(10,(it+24)/8.0));
#else
  {
  FMM.printNow = true;
#if BUILD
  numBodies = 10000000;
#else
  numBodies = 1000000;
#endif
#endif // MANY
  std::cout << "N                    : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.cube(bodies);
  FMM.startTimer("FMM");
#if BOTTOMUP
  FMM.bottomup(bodies,cells);
#else
  FMM.topdown(bodies,cells);
#endif
#if BUILD
#else
  FMM.startPAPI();
#if IneJ
  FMM.evaluate(cells,cells);
#else
  FMM.evaluate(cells);
#endif
  FMM.stopPAPI();
  FMM.stopTimer("FMM",true);
  FMM.eraseTimer("FMM");
  FMM.writeTime();
  FMM.resetTimer();

  bodies2 = bodies;
  if (bodies2.size() > 100) bodies2.resize(100);
  DATA.initTarget(bodies2);
  FMM.startTimer("Direct sum");
  FMM.direct(bodies2,bodies);
  FMM.stopTimer("Direct sum",true);
  FMM.eraseTimer("Direct sum");
  if (bodies.size() > 100) bodies.resize(100);
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  DATA.printError(diff1,norm1,diff2,norm2);
#endif // BUILD
  }
}
