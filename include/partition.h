#ifndef partition_h
#define partition_h
#include "logger.h"
#include "sort.h"

//! Handles all the partitioning of domains
class Partition {
private:
  const int mpirank;                                            //!< Rank of MPI communicator
  const int mpisize;                                            //!< Size of MPI communicator
  const int numBins;                                            //!< Number of sampling bins
  int numLevels;                                                //!< Levels of MPI rank binary tree
  int * rankDispl;                                              //!< Displacement of MPI rank group
  int * rankCount;                                              //!< Size of MPI rank group
  int * rankColor;                                              //!< Color of MPI rank group
  int * rankKey;                                                //!< Key of MPI rank group
  int * rankMap;                                                //!< Map to MPI rank group
  int * sendDispl;                                              //!< Displacement of bodies to send per rank
  int * sendCount;                                              //!< Count of bodies to send per rank
  int * scanHist;                                               //!< Scan of histogram
  int * localHist;                                              //!< Local histogram
  int * globalHist;                                             //!< Global histogram
  Bounds * rankBounds;                                          //!< Bounds of each rank
  Bodies buffer;                                                //!< MPI communication buffer for bodies

private:

public:
  //! Constructor
  Partition(int _mpirank, int _mpisize) : mpirank(_mpirank), mpisize(_mpisize), numBins(16) {
    rankDispl  = new int [mpisize];                             // Allocate displacement of MPI rank group
    rankCount  = new int [mpisize];                             // Allocate size of MPI rank group
    rankColor  = new int [mpisize];                             // Allocate color of MPI rank group
    rankKey    = new int [mpisize];                             // Allocate key of MPI rank group
    rankMap    = new int [mpisize];                             // Allocate map to MPI rank group
    sendDispl  = new int [mpisize];                             // Allocate displacement of bodies to send per rank
    sendCount  = new int [mpisize];                             // Allocate count of bodies to send per rank
    scanHist   = new int [numBins];                             // Allocate scan of histogram
    localHist  = new int [numBins];                             // Allocate local histogram
    globalHist = new int [numBins];                             // Allocate global histogram
    rankBounds = new Bounds [mpisize];                          // Allocate bounds of each rank
    numLevels = 0;                                              // Initialize levels of MPI rank binary tree
    int size = mpisize - 1;                                     // Initialize size of MPI rank binary tree
    while (size > 0) {                                          // While there are bits to process
      size >>= 1;                                               //  Bitshift MPI size
      numLevels++;                                              //  Increment levels of MPI rank binary tree
    }                                                           // End while for bits to process
  }

//! Destructor
  ~Partition() {
    delete[] rankDispl;                                         // Deallocate displacement of MPI rank group
    delete[] rankCount;                                         // Deallocate size of MPI rank group
    delete[] rankColor;                                         // Deallocate color of MPI rank group
    delete[] rankKey;                                           // Deallocate key of MPI rank group
    delete[] rankMap;                                           // Deallocate ,ap to MPI rank group
    delete[] sendDispl;                                         // Deallocate displacement of bodies to send per rank
    delete[] sendCount;                                         // Deallocate count of bodies to send per rank
    delete[] scanHist;                                          // Deallocate scan of histogram
    delete[] localHist;                                         // Deallocate local histogram
    delete[] globalHist;                                        // Deallocate global histogram
    delete[] rankBounds;                                        // Deallocate bounds of each rank
  }

  //! Partitioning by recursive bisection
  Bounds bisection(Bodies & bodies, Bounds globalBounds) {
    logger::startTimer("Partition");                            // Start timer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over MPI ranks
      rankDispl[irank] = 0;                                     //  Initialize displacement of MPI rank group
      rankCount[irank] = mpisize;                               //  Initialize size of MPI rank group
      rankColor[irank] = 0;                                     //  Initialize color of MPI rank group
      rankKey[irank] = 0;                                       //  Initialize key of MPI rank group
      rankMap[irank] = 0;                                       //  Initialize map to MPI rank group
      sendDispl[irank] = 0;                                     //  Initialize displacement of bodies to send per rank
      sendCount[irank] = bodies.size();                         //  Initialize count of bodies to send per rank
      rankBounds[irank] = globalBounds;                         //  Initialize bounds of each rank
    }                                                           // End loop over MPI ranks
    buffer = bodies;                                            // Resize sort buffer
    for (int level=0; level<numLevels; level++) {
      int numPartitions = rankColor[mpisize-1] + 1;
      for (int ipart=0; ipart<numPartitions; ipart++) {
	int irank = rankMap[ipart];
	Bounds bounds = rankBounds[irank];
	int direction = 0;
	real_t length = 0;
	for (int d=0; d<3; d++) {
	  if (length < (bounds.Xmax[d] - bounds.Xmin[d])) {
	    direction = d;
	    length = (bounds.Xmax[d] - bounds.Xmin[d]);
	  }
	}
	int rankSplit = rankCount[irank] / 2;
	int oldRankCount = rankCount[irank];
	int numLocalBodies = sendCount[irank];
	int bodyBegin = 0;
	int bodyEnd = numLocalBodies;
	int numGlobalBodies;
	MPI_Allreduce(&numLocalBodies, &numGlobalBodies, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	int globalSplit = numGlobalBodies * rankSplit / oldRankCount;
	int globalOffset = 0;
	real_t xmax = bounds.Xmax[direction];
	real_t xmin = bounds.Xmin[direction];
	real_t dx = (xmax - xmin) / numBins;
	B_iter B = bodies.begin() + sendDispl[irank];
	if (globalSplit > 0) {
	  for (int binRefine=0; binRefine<3; binRefine++) {
	    for (int ibin=0; ibin<numBins; ibin++) {
	      localHist[ibin] = 0;
	    }
	    for (int b=bodyBegin; b<bodyEnd; b++) {
	      real_t x = B[b].X[direction];
	      int ibin = (x - xmin + EPS) / (dx + EPS);
	      localHist[ibin]++;
	    }
	    scanHist[0] = localHist[0];
	    for (int ibin=1; ibin<numBins; ibin++) {
	      scanHist[ibin] = scanHist[ibin-1] + localHist[ibin];
	    }
	    for (int b=bodyEnd-1; b>=bodyBegin; b--) {
	      real_t x = B[b].X[direction];
	      int ibin = (x - xmin + EPS) / (dx + EPS);
	      scanHist[ibin]--;
	      int bnew = scanHist[ibin] + bodyBegin;
	      buffer[bnew] = B[b];
	    }
	    for (int b=bodyBegin; b<bodyEnd; b++) {
	      B[b] = buffer[b];
	    }
	    MPI_Allreduce(localHist, globalHist, numBins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
	      scanHist[ibin] = scanHist[ibin-1] + localHist[ibin-1];
	    }
	    bodyBegin += scanHist[splitBin];
	    bodyEnd = bodyBegin + localHist[splitBin];
	  }
	}
	int rankBegin = rankDispl[irank];
	int rankEnd = rankBegin + rankCount[irank];
	for (irank=rankBegin; irank<rankEnd; irank++) {
	  rankSplit = rankCount[irank] / 2;
	  if (irank - rankDispl[irank] < rankSplit) {
	    rankCount[irank] = rankSplit;
	    rankColor[irank] = rankColor[irank] * 2;
	    rankBounds[irank].Xmax[direction] = xmin;
	    sendCount[irank] = bodyBegin;
	  } else {
	    rankDispl[irank] += rankSplit;
	    rankCount[irank] -= rankSplit;
	    rankColor[irank] = rankColor[irank] * 2 + 1;
	    rankBounds[irank].Xmin[direction] = xmin;
	    sendDispl[irank] += bodyBegin;
	    sendCount[irank] -= bodyBegin;
	  }
	  if (level == numLevels-1) rankColor[irank] = rankDispl[irank];
	  rankKey[irank] = irank - rankDispl[irank];
	}
      }
      int ipart = 0;
      for (int irank=0; irank<mpisize; irank++) {
	if (rankKey[irank] == 0) {
	  rankMap[ipart] = rankDispl[irank];
	  ipart++;
	}
      }
    }
    B_iter B = bodies.begin();
    for (int irank=0; irank<mpisize; irank++) {
      int bodyBegin = sendDispl[irank];
      int bodyEnd = bodyBegin + sendCount[irank];
      for (int b=bodyBegin; b<bodyEnd; b++, B++) {
	B->IRANK = irank;
      }
    }
    logger::stopTimer("Partition");                             // Stop timer
    return rankBounds[mpirank];                                 // Return local bounds
  }

  //! Partition bodies with geometric octsection
  Bounds octsection(Bodies & bodies, Bounds global) {
    logger::startTimer("Partition");                            // Start timer
    int size = mpisize;                                         // Initialize MPI size counter
    int Npartition[3] = {1, 1, 1};                              // Number of partitions in each direction
    int d = 0;                                                  // Initialize dimension counter
    while (size != 1) {                                         // Divide domain while counter is not one
      Npartition[d] <<= 1;                                      //  Divide this dimension
      d = (d+1) % 3;                                            //  Increment dimension
      size >>= 1;                                               //  Right shift the bits of counter
    }                                                           // End while loop for domain subdivision
    real_t Xpartition[3];                                       // Size of partitions in each direction
    for (d=0; d<3; d++) {                                       // Loop over dimensions
      Xpartition[d] = (global.Xmax[d] - global.Xmin[d]) / Npartition[d];//  Size of partition in each direction
    }                                                           // End loop over dimensions
    int iX[3];                                                  // Index vector
    iX[0] = mpirank % Npartition[0];                            // x index of partition
    iX[1] = mpirank / Npartition[0] % Npartition[1];            // y index
    iX[2] = mpirank / Npartition[0] / Npartition[1];            // z index
    Bounds local;                                               // Local bounds
    for (d=0; d<3; d++) {                                       // Loop over dimensions
      local.Xmin[d] = global.Xmin[d] + iX[d] * Xpartition[d];   // Xmin of local domain at current rank
      local.Xmax[d] = global.Xmin[d] + (iX[d] + 1) * Xpartition[d];// Xmax of local domain at current rank
    }                                                           // End loop over dimensions
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      for (d=0; d<3; d++) {                                     //  Loop over dimensions
        iX[d] = int((B->X[d] - global.Xmin[d]) / Xpartition[d]);//   Index vector of partition
      }                                                         //  End loop over dimensions
      B->IRANK = iX[0] + Npartition[0] * (iX[1] + iX[2] * Npartition[1]);//  Set send rank
      assert(0 <= B->IRANK && B->IRANK < mpisize);
    }                                                           // End loop over bodies
    logger::stopTimer("Partition");                             // Stop timer
    logger::startTimer("Sort");                                 // Start timer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.irank(bodies);                                // Sort bodies according to IRANK
    logger::stopTimer("Sort");                                  // Stop timer
    return local;
  }

  //! Send bodies back to where they came from
  void unpartition(Bodies & bodies) {
    logger::startTimer("Sort");                                 // Start timer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.irank(bodies);                                // Sort bodies according to IRANK
    logger::stopTimer("Sort");                                  // Stop timer
  }
};
#endif
