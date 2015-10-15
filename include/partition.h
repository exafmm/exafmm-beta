#ifndef partition_h
#define partition_h
#include "logger.h"
#include "sort.h"

namespace exafmm {
  //! Handles all the partitioning of domains
  class Partition {
  private:
    const int mpirank;                                          //!< Rank of MPI communicator
    const int mpisize;                                          //!< Size of MPI communicator
    const int numBins;                                          //!< Number of sampling bins
    int numLevels;                                              //!< Levels of MPI rank binary tree
    int * rankDispl;                                            //!< Displacement of MPI rank group
    int * rankCount;                                            //!< Size of MPI rank group
    int * rankColor;                                            //!< Color of MPI rank group
    int * rankKey;                                              //!< Key of MPI rank group
    int * rankMap;                                              //!< Map partition to MPI rank group
    int * sendDispl;                                            //!< Displacement of bodies to send per rank
    int * sendCount;                                            //!< Count of bodies to send per rank
    int * scanHist;                                             //!< Scan of histogram
    int * countHist;                                            //!< Body count histogram
    float * weightHist;                                         //!< Body weight histogram
    float * globalHist;                                         //!< Global body weight histogram
    Bounds * rankBounds;                                        //!< Bounds of each rank
    Bodies buffer;                                              //!< MPI communication buffer for bodies

  private:

  public:
    //! Constructor
    Partition(int _mpirank, int _mpisize) : mpirank(_mpirank), mpisize(_mpisize), numBins(16) {
      rankDispl  = new int [mpisize];                           // Allocate displacement of MPI rank group
      rankCount  = new int [mpisize];                           // Allocate size of MPI rank group
      rankColor  = new int [mpisize];                           // Allocate color of MPI rank group
      rankKey    = new int [mpisize];                           // Allocate key of MPI rank group
      rankMap    = new int [mpisize];                           // Allocate map for partition to MPI rank group
      sendDispl  = new int [mpisize];                           // Allocate displacement of bodies to send per rank
      sendCount  = new int [mpisize];                           // Allocate count of bodies to send per rank
      scanHist   = new int [numBins];                           // Allocate scan of histogram
      countHist  = new int [numBins];                           // Allocate body count histogram
      weightHist = new float [numBins];                         // Allocate body weight histogram
      globalHist = new float [numBins];                         // Allocate global body weight histogram
      rankBounds = new Bounds [mpisize];                        // Allocate bounds of each rank
      numLevels = 0;                                            // Initialize levels of MPI rank binary tree
      int size = mpisize - 1;                                   // Initialize size of MPI rank binary tree
      while (size > 0) {                                        // While there are bits to process
	size >>= 1;                                             //  Bitshift MPI size
	numLevels++;                                            //  Increment levels of MPI rank binary tree
      }                                                         // End while for bits to process
    }

    //! Destructor
    ~Partition() {
      delete[] rankDispl;                                       // Deallocate displacement of MPI rank group
      delete[] rankCount;                                       // Deallocate size of MPI rank group
      delete[] rankColor;                                       // Deallocate color of MPI rank group
      delete[] rankKey;                                         // Deallocate key of MPI rank group
      delete[] rankMap;                                         // Deallocate map for partition to MPI rank group
      delete[] sendDispl;                                       // Deallocate displacement of bodies to send per rank
      delete[] sendCount;                                       // Deallocate count of bodies to send per rank
      delete[] scanHist;                                        // Deallocate scan of histogram
      delete[] countHist;                                       // Deallocate body count histogram
      delete[] weightHist;                                      // Deallocate body weight histogram
      delete[] globalHist;                                      // Deallocate global body weight histogram
      delete[] rankBounds;                                      // Deallocate bounds of each rank
    }

    //! Partitioning by orthogonal recursive bisection
    Bounds bisection(Bodies & bodies, Bounds globalBounds) {
      logger::startTimer("Partition");                          // Start timer
      for (int irank=0; irank<mpisize; irank++) {               // Loop over MPI ranks
	rankDispl[irank] = 0;                                   //  Initialize displacement of MPI rank group
	rankCount[irank] = mpisize;                             //  Initialize size of MPI rank group
	rankColor[irank] = 0;                                   //  Initialize color of MPI rank group
	rankKey[irank] = 0;                                     //  Initialize key of MPI rank group
	rankMap[irank] = 0;                                     //  Initialize map to MPI rank group
	sendDispl[irank] = 0;                                   //  Initialize displacement of bodies to send per rank
	sendCount[irank] = bodies.size();                       //  Initialize count of bodies to send per rank
	rankBounds[irank] = globalBounds;                       //  Initialize bounds of each rank
      }                                                         // End loop over MPI ranks
      buffer = bodies;                                          // Resize sort buffer
      for (int level=0; level<numLevels; level++) {             // Loop over levels of MPI rank binary tree
	int numPartitions = rankColor[mpisize-1] + 1;           //  Number of partitions in current level
	for (int ipart=0; ipart<numPartitions; ipart++) {       //  Loop over partitions
	  int irank = rankMap[ipart];                           //   Map partition to MPI rank
	  Bounds bounds = rankBounds[irank];                    //   Bounds of current rank
	  int direction = 0;                                    //   Initialize direction of partitioning
	  real_t length = 0;                                    //   Initialize length of partition
	  for (int d=0; d<3; d++) {                             //   Loop over dimensions
	    if (length < (bounds.Xmax[d] - bounds.Xmin[d])) {   //    If present dimension is longer
	      direction = d;                                    //     Update direction to current dimension
	      length = (bounds.Xmax[d] - bounds.Xmin[d]);       //     Update length to current dimension
	    }                                                   //    End if for longer dimension
	  }                                                     //   End loop over dimensions
	  int rankSplit = rankCount[irank] / 2;                 //   MPI rank splitter
	  int oldRankCount = rankCount[irank];                  //   Old MPI rank count
	  int bodyBegin = 0;                                    //   Initialize body begin index
	  int bodyEnd = sendCount[irank];                       //   Initialize body end index
	  B_iter B = bodies.begin() + sendDispl[irank];         //   Body begin iterator of current partition
	  float localWeightSum = 0;                             //   Initialize local sum of weights in current partition
	  for (int b=bodyBegin; b<bodyEnd; b++) {               //    Loop over bodies in current partition
	    localWeightSum += B[b].WEIGHT;                      //     Add weights of body to local sum
	  }                                                     //    End loop over bodies in current partition
	  float globalWeightSum;                                //   Declare global sum of weights in current partition
	  MPI_Allreduce(&localWeightSum, &globalWeightSum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Reduce sum of weights
	  float globalSplit = globalWeightSum * rankSplit / oldRankCount;// Global weight splitter index
	  float globalOffset = 0;                               //   Initialize global weight offset
	  real_t xmax = bounds.Xmax[direction];                 //   Upper bound of partition
	  real_t xmin = bounds.Xmin[direction];                 //   Lower bound of partition
	  real_t dx = (xmax - xmin) / numBins;                  //   Length of bins
	  if (rankSplit > 0) {                                  //   If the partition requires splitting
	    for (int binRefine=0; binRefine<3; binRefine++) {   //    Loop for bin refinement
	      for (int ibin=0; ibin<numBins; ibin++) {          //     Loop over bins
		countHist[ibin] = 0;                            //      Initialize body count histogram
		weightHist[ibin] = 0;                           //      Initialize body weight histogram
	      }                                                 //     End loop over bins
	      for (int b=bodyBegin; b<bodyEnd; b++) {           //     Loop over bodies
		real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
		int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
		countHist[ibin]++;                              //      Increment body count histogram
		weightHist[ibin] += B[b].WEIGHT;                //      Increment body weight histogram
	      }                                                 //     End loop over bodies
	      scanHist[0] = countHist[0];                       //     Initialize scan array
	      for (int ibin=1; ibin<numBins; ibin++) {          //     Loop over bins
		scanHist[ibin] = scanHist[ibin-1] + countHist[ibin]; //   Inclusive scan for body count histogram
	      }                                                 //     End loop over bins
	      for (int b=bodyEnd-1; b>=bodyBegin; b--) {        //     Loop over bodies backwards
		real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
		int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
		scanHist[ibin]--;                               //      Decrement scanned histogram
		int bnew = scanHist[ibin] + bodyBegin;          //      New index of sorted body
		buffer[bnew] = B[b];                            //      Copy body to sort buffer
	      }                                                 //     End loop over bodies backwards
	      for (int b=bodyBegin; b<bodyEnd; b++) {           //     Loop over bodies
		B[b] = buffer[b];                               //      Copy back bodies from buffer
	      }                                                 //     End loop over bodies
	      MPI_Allreduce(weightHist, globalHist, numBins, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Reduce weight histogram
	      int splitBin = 0;                                 //     Initialize bin splitter
	      while (globalOffset < globalSplit) {              //     While scan of global histogram is less than splitter
		globalOffset += globalHist[splitBin];           //      Scan global histogram
		splitBin++;                                     //      Increment bin count
	      }                                                 //     End while for global histogram scan
	      splitBin--;                                       //     Move back one bin
	      globalOffset -= globalHist[splitBin];             //     Decrement offset accordingly
	      xmax = xmin + (splitBin + 1) * dx;                //     Zoom in to current bin by redefining upper
	      xmin = xmin + splitBin * dx;                      //     and lower bounds of partition
	      dx = (xmax - xmin) / numBins;                     //     Update length of partition accordingly
	      scanHist[0] = 0;                                  //     Initialize scan array
	      for (int ibin=1; ibin<numBins; ibin++) {          //     Loop over bins
		scanHist[ibin] = scanHist[ibin-1] + countHist[ibin-1]; // Exclusive scan of body count histogram
	      }                                                 //     End loop over bins
	      bodyBegin += scanHist[splitBin];                  //     Update body begin index
	      bodyEnd = bodyBegin + countHist[splitBin];        //     Update body end index
	    }                                                   //    End loop for bin refinement
	  }                                                     //   End if for splitting partition
	  int rankBegin = rankDispl[irank];                     //   Save current range of MPI ranks 
	  int rankEnd = rankBegin + rankCount[irank];           //   so that they don't get overwritten
	  for (irank=rankBegin; irank<rankEnd; irank++) {       //   Loop over current range of MPI ranks
	    rankSplit = rankCount[irank] / 2;                   //    MPI rank splitter
	    if (irank - rankDispl[irank] < rankSplit) {         //    If on left side of splitter
	      rankCount[irank] = rankSplit;                     //     Count is the splitter index
	      rankColor[irank] = rankColor[irank] * 2;          //     Color is doubled
	      rankBounds[irank].Xmax[direction] = xmin;         //     Update right bound with splitter
	      sendCount[irank] = bodyBegin;                     //     Send body count is the begin index
	    } else {                                            //    If on right side of splitter
	      rankDispl[irank] += rankSplit;                    //     Update displacement with splitter index
	      rankCount[irank] -= rankSplit;                    //     Count is remainder of splitter index
	      rankColor[irank] = rankColor[irank] * 2 + 1;      //     Color is doubled plus one
	      rankBounds[irank].Xmin[direction] = xmin;         //     Update left bound
	      sendDispl[irank] += bodyBegin;                    //     Increment send body displacement
	      sendCount[irank] -= bodyBegin;                    //     Decrement send body count
	    }                                                   //    End if for side of splitter
	    if (level == numLevels-1) rankColor[irank] = rankDispl[irank]; // Special case for final rank color
	    rankKey[irank] = irank - rankDispl[irank];          //    Rank key is determined from rank displacement
	  }                                                     //   End loop over current range of MPI ranks
	}                                                       //  End loop over partitions
	int ipart = 0;                                          //  Initialize partition index
	for (int irank=0; irank<mpisize; irank++) {             //  Loop over MPI ranks
	  if (rankKey[irank] == 0) {                            //   If rank key is zero it's a new group
	    rankMap[ipart] = rankDispl[irank];                  //    Update rank map with rank displacement
	    ipart++;                                            //    Increment partition index
	  }                                                     //   End if for new group
	}                                                       //  End loop over MPI ranks
      }                                                         // End loop over levels of MPI rank binary tree
      B_iter B = bodies.begin();                                // Body begin iterator
      for (int irank=0; irank<mpisize; irank++) {               // Loop over MPI ranks
	int bodyBegin = sendDispl[irank];                       //  Body begin index for current rank
	int bodyEnd = bodyBegin + sendCount[irank];             //  Body end index for current rank
	for (int b=bodyBegin; b<bodyEnd; b++, B++) {            //  Loop over bodies in current rank
	  B->IRANK = irank;                                     //   Copy MPI rank to body
	}                                                       //  End loop over bodies in current rank
      }                                                         // End loop over MPI ranks
      logger::stopTimer("Partition");                           // Stop timer
      return rankBounds[mpirank];                               // Return local bounds
    }

    //! Partition bodies with geometric octsection
    Bounds octsection(Bodies & bodies, Bounds global) {
      logger::startTimer("Partition");                          // Start timer
      int size = mpisize;                                       // Initialize MPI size counter
      int Npartition[3] = {1, 1, 1};                            // Number of partitions in each direction
      int d = 0;                                                // Initialize dimension counter
      while (size != 1) {                                       // Divide domain while counter is not one
	Npartition[d] <<= 1;                                    //  Divide this dimension
	d = (d+1) % 3;                                          //  Increment dimension
	size >>= 1;                                             //  Right shift the bits of counter
      }                                                         // End while loop for domain subdivision
      real_t Xpartition[3];                                     // Size of partitions in each direction
      for (d=0; d<3; d++) {                                     // Loop over dimensions
	Xpartition[d] = (global.Xmax[d] - global.Xmin[d]) / Npartition[d];//  Size of partition in each direction
      }                                                         // End loop over dimensions
      int iX[3];                                                // Index vector
      iX[0] = mpirank % Npartition[0];                          // x index of partition
      iX[1] = mpirank / Npartition[0] % Npartition[1];          // y index
      iX[2] = mpirank / Npartition[0] / Npartition[1];          // z index
      Bounds local;                                             // Local bounds
      for (d=0; d<3; d++) {                                     // Loop over dimensions
	local.Xmin[d] = global.Xmin[d] + iX[d] * Xpartition[d]; // Xmin of local domain at current rank
	local.Xmax[d] = global.Xmin[d] + (iX[d] + 1) * Xpartition[d];// Xmax of local domain at current rank
      }                                                         // End loop over dimensions
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	for (d=0; d<3; d++) {                                   //  Loop over dimensions
	  iX[d] = int((B->X[d] - global.Xmin[d]) / Xpartition[d]);//   Index vector of partition
	}                                                       //  End loop over dimensions
	B->IRANK = iX[0] + Npartition[0] * (iX[1] + iX[2] * Npartition[1]);//  Set send rank
	assert(0 <= B->IRANK && B->IRANK < mpisize);
      }                                                         // End loop over bodies
      logger::stopTimer("Partition");                           // Stop timer
      logger::startTimer("Sort");                               // Start timer
      Sort sort;                                                // Instantiate sort class
      bodies = sort.irank(bodies);                              // Sort bodies according to IRANK
      logger::stopTimer("Sort");                                // Stop timer
      return local;
    }

    //! Send bodies back to where they came from
    void unpartition(Bodies & bodies) {
      logger::startTimer("Sort");                               // Start timer
      Sort sort;                                                // Instantiate sort class
      bodies = sort.irank(bodies);                              // Sort bodies according to IRANK
      logger::stopTimer("Sort");                                // Stop timer
    }
  };
}
#endif
