#include "dataset.h"
#include "serialfmm.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset DATA;
  SerialFMM FMM;
  FMM.printNow = true;
#if HYBRID
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
  numBodies = int(pow(10,(it+24)/8.0));
#else
  {
  numBodies = 1000000;
#endif // MANY
  if(FMM.printNow) std::cout << "N                    : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.cube(bodies);
  FMM.startTimer("FMM");
  FMM.setBounds(bodies);
  FMM.buildTree(bodies,cells);
  FMM.upwardPass(cells);
  FMM.startPAPI();
#if IneJ
  FMM.evaluate(cells,cells);
#else
  FMM.evaluate(cells);
#endif
  FMM.stopPAPI();
  FMM.downwardPass(cells);
  FMM.stopTimer("FMM",FMM.printNow);
  FMM.eraseTimer("FMM");
  FMM.writeTime();
  FMM.resetTimer();

  jbodies = bodies;
  if (bodies.size() > 100) bodies.resize(100);
  Bodies bodies2 = bodies;
  DATA.initTarget(bodies2);
  FMM.startTimer("Direct sum");
  FMM.direct(bodies2,jbodies);
  FMM.normalize(bodies2);
  FMM.stopTimer("Direct sum",FMM.printNow);
  FMM.eraseTimer("Direct sum");
  real_t diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  if(FMM.printNow) DATA.printError(diff1,norm1,diff2,norm2);
  }
}
