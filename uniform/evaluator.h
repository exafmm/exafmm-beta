#include "kernel.h"
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
    double tic, toc;

    int rankOffset = 13 * numCells;
    for( int i=0; i<numCells; i++ ) {
      for_m Multipole[i+rankOffset][m] = 0;
      for_l Local[i][l] = 0;
    }

    tic = getTime();
    P2M();
    toc = getTime();
    if( printNow ) printf("P2M     : %lf : %f GFlops\n",toc-tic,148.*numBodies/(toc-tic)*1e-9);

    tic = getTime();
    M2M();
    toc = getTime();
    if( printNow ) printf("M2M     : %lf : %f GFlops\n",toc-tic,955.*numCells/(toc-tic)*1e-9);
  }

  void downwardPass() {
    double tic, toc;

    logger::startTimer("Traverse");
    tic = getTime();
    M2L();
    toc = getTime();
    if( printNow ) printf("M2L     : %lf : %f GFlops\n",toc-tic,2417.*numCells*189/(toc-tic)*1e-9);
    logger::stopTimer("Traverse");

    logger::startTimer("Downward pass");
    tic = getTime();
    L2L();
    toc = getTime();
    if( printNow ) printf("L2L     : %lf : %f GFlops\n",toc-tic,1902.*numCells/(toc-tic)*1e-9);

    tic = getTime();
    L2P();
    toc = getTime();
    if( printNow ) printf("L2P     : %lf : %f GFlops\n",toc-tic,558.*numBodies/(toc-tic)*1e-9);
    logger::stopTimer("Downward pass");

#if DO_P2P 
    logger::startTimer("Traverse");
    tic = getTime();
    P2P();
    toc = getTime();
    if( printNow ) printf("P2P     : %lf : %f GFlops\n",toc-tic,22.*numBodies*numBodies*27/(toc-tic)/numLeafs*1e-9);
    logger::stopTimer("Traverse");
#endif
  }
};
