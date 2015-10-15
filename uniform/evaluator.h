#include "kernels.h"
//#include "kernels_smooth.h"
#include "logger.h"

namespace exafmm {
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

    void downwardPass() {
      logger::startTimer("Traverse");
      M2L();
      logger::stopTimer("Traverse", 0);

      logger::startTimer("Downward pass");
      L2L();
      L2P();
      logger::stopTimer("Downward pass");

      logger::startTimer("Traverse");
      P2P();
      logger::stopTimer("Traverse");
    }
  };
}
