#include "kernels.h"
#include "logger.h"

class Evaluator : public Kernel {
public:
  bool printNow;

public:
  Evaluator() : printNow(true) {}

  double getTime() const {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec+tv.tv_usec*1e-6);
  }

  void upwardPass() const {
    int rankOffset = 13 * numCells;
    for( int i=0; i<numCells; i++ ) {
      for_m Multipole[i+rankOffset][m] = 0;
      for_l Local[i][l] = 0;
    }
    P2M();
    M2M();
  }

  void downwardPass(int *common_stencil, int *far_stencil,
                    int *near_stencil) {
    logger::startTimer("Traverse");
    //M2L();
    M2L(common_stencil, far_stencil);
    logger::stopTimer("Traverse", 0);

    logger::startTimer("Downward pass");
    L2L();
    L2P();
    logger::stopTimer("Downward pass");

    logger::startTimer("Traverse");
    //P2P();
    P2P(near_stencil);
    logger::stopTimer("Traverse");
  }
};
