#ifndef partition_h
#define partition_h
#include "logger.h"
#include "sort.h"

//! Handles all the partitioning of domains
class Partition {
private:
  const int mpirank;
  const int mpisize;
  const int numBins;
  int numLevels;
  int level;
  int numPartitions;
  int ipart;
  int rank;
  int binRefine;
  int direction;
  int rankSplit;
  int oldRankCount;
  int bodyBegin;
  int bodyEnd;
  float localWeightSum;
  float globalWeightSum;
  float globalSplit;
  float globalOffset;
  real_t xmax;
  real_t xmin;
  real_t dx;
  int * rankDispl;
  int * rankCount;
  int * rankColor;
  int * rankKey;
  int * rankMap;
  int * sendDispl;
  int * sendCount;
  int * scanHist;
  int * countHist;
  float * weightHist;
  float * globalHist;
  Bounds bounds;
  Bounds * rankBounds;
  Bodies buffer;
  B_iter B00, B0;

private:

public:
  //! Constructor
  Partition(int _mpirank, int _mpisize) : mpirank(_mpirank), mpisize(_mpisize), numBins(16) {
    rankDispl  = new int [mpisize];
    rankCount  = new int [mpisize];
    rankColor  = new int [mpisize];
    rankKey    = new int [mpisize];
    rankMap    = new int [mpisize];
    sendDispl  = new int [mpisize];
    sendCount  = new int [mpisize];
    scanHist   = new int [numBins];
    countHist  = new int [numBins];
    weightHist = new float [numBins];
    globalHist = new float [numBins];
    rankBounds = new Bounds [mpisize];
    numLevels = 0;
    int size = mpisize - 1;
    while (size > 0) {
      size >>= 1;
      numLevels++;
    }
  }

//! Destructor
  ~Partition() {
    delete[] rankDispl;
    delete[] rankCount;
    delete[] rankColor;
    delete[] rankKey;
    delete[] rankMap;
    delete[] sendDispl;
    delete[] sendCount;
    delete[] scanHist;
    delete[] countHist;
    delete[] weightHist;
    delete[] globalHist;
    delete[] rankBounds;
  }

  //! Partitioning by orthogonal recursive bisection
  Bounds bisection(Bodies & bodies, Bounds globalBounds) {
    B00 = bodies.begin();
    logger::startTimer("Partition");
    for (int irank=0; irank<mpisize; irank++) {
      rankDispl[irank] = 0;
      rankCount[irank] = mpisize;
      rankColor[irank] = 0;
      rankKey[irank] = 0;
      rankMap[irank] = 0;
      sendDispl[irank] = 0;
      sendCount[irank] = bodies.size();
      rankBounds[irank] = globalBounds;
    }
    buffer = bodies;
    loopLevels();
    B_iter B = bodies.begin();
    for (int irank=0; irank<mpisize; irank++) {
      int bodyBegin = sendDispl[irank];
      int bodyEnd = bodyBegin + sendCount[irank];
      for (int b=bodyBegin; b<bodyEnd; b++, B++) {
	B->IRANK = irank;
      }
    }
    logger::stopTimer("Partition");
    return rankBounds[mpirank];
  }

  void loopLevels() {
    for (level=0; level<numLevels; level++) {
      numPartitions = rankColor[mpisize-1] + 1;
      for (ipart=0; ipart<numPartitions; ipart++) {
	rank = rankMap[ipart];
	bounds = rankBounds[rank];
	direction = 0;
	real_t length = 0;
	for (int d=0; d<3; d++) {
	  if (length < (bounds.Xmax[d] - bounds.Xmin[d])) {
	    direction = d;
	    length = (bounds.Xmax[d] - bounds.Xmin[d]);
	  }
	}
	rankSplit = rankCount[rank] / 2;
	oldRankCount = rankCount[rank];
	bodyBegin = 0;
	bodyEnd = sendCount[rank];
	B0 = B00 + sendDispl[rank];
	localWeightSum = 0;
	for (int b=bodyBegin; b<bodyEnd; b++) {
	  localWeightSum += B0[b].WEIGHT;
	}
	MPI_Allreduce(&localWeightSum, &globalWeightSum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	globalSplit = globalWeightSum * rankSplit / oldRankCount;
	globalOffset = 0;
	xmax = bounds.Xmax[direction];
	xmin = bounds.Xmin[direction];
	dx = (xmax - xmin) / numBins;
	if (rankSplit > 0) {
	  for (binRefine=0; binRefine<3; binRefine++) {
	    for (int ibin=0; ibin<numBins; ibin++) {
	      countHist[ibin] = 0;
	      weightHist[ibin] = 0;
	    }
	    for (int b=bodyBegin; b<bodyEnd; b++) {
	      real_t x = B0[b].X[direction];
	      int ibin = (x - xmin + EPS) / (dx + EPS);
	      countHist[ibin]++;
	      weightHist[ibin] += B0[b].WEIGHT;
	    }
	    scanHist[0] = countHist[0];
	    for (int ibin=1; ibin<numBins; ibin++) {
	      scanHist[ibin] = scanHist[ibin-1] + countHist[ibin];
	    }
	    for (int b=bodyEnd-1; b>=bodyBegin; b--) {
	      real_t x = B0[b].X[direction];
	      int ibin = (x - xmin + EPS) / (dx + EPS);
	      scanHist[ibin]--;
	      int bnew = scanHist[ibin] + bodyBegin;
	      buffer[bnew] = B0[b];
	    }
	    for (int b=bodyBegin; b<bodyEnd; b++) {
	      B0[b] = buffer[b];
	    }
	    MPI_Allreduce(weightHist, globalHist, numBins, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	    int splitBin = 0;
	    while (globalOffset < globalSplit) {
	      globalOffset += globalHist[splitBin];
	      splitBin++;
	    }
	    splitBin--;
	    globalOffset -= globalHist[splitBin];
	    xmax = xmin + (splitBin + 1) * dx;
	    xmin = xmin + splitBin * dx;
	    dx = (xmax - xmin) / numBins;
	    scanHist[0] = 0;
	    for (int ibin=1; ibin<numBins; ibin++) {
	      scanHist[ibin] = scanHist[ibin-1] + countHist[ibin-1];
	    }
	    bodyBegin += scanHist[splitBin];
	    bodyEnd = bodyBegin + countHist[splitBin];
	  }
	}
	int rankBegin = rankDispl[rank];
	int rankEnd = rankBegin + rankCount[rank];
	for (rank=rankBegin; rank<rankEnd; rank++) {
	  rankSplit = rankCount[rank] / 2;
	  if (rank - rankDispl[rank] < rankSplit) {
	    rankCount[rank] = rankSplit;
	    rankColor[rank] = rankColor[rank] * 2;
	    rankBounds[rank].Xmax[direction] = xmin;
	    sendCount[rank] = bodyBegin;
	  } else {
	    rankDispl[rank] += rankSplit;
	    rankCount[rank] -= rankSplit;
	    rankColor[rank] = rankColor[rank] * 2 + 1;
	    rankBounds[rank].Xmin[direction] = xmin;
	    sendDispl[rank] += bodyBegin;
	    sendCount[rank] -= bodyBegin;
	  }
	  if (level == numLevels-1) rankColor[rank] = rankDispl[rank];
	  rankKey[rank] = rank - rankDispl[rank];
	}
      }
      ipart = 0;
      for (int irank=0; irank<mpisize; irank++) {
	if (rankKey[irank] == 0) {
	  rankMap[ipart] = rankDispl[irank];
	  ipart++;
	}
      }
    }
  }
};
#endif
