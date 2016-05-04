#ifndef partition_h
#define partition_h
#include "logger.h"
#include "sort.h"
#include "hilbert_key.h"
#include <numeric>

namespace exafmm {
//! Handles all the partitioning of domains
class Partition {
  typedef std::vector<int64_t> VHilbert;                        //!< Type of Hilbert key vectors
  typedef std::pair<int64_t, int64_t> KeyPair;                  //!< Type of Hilbert Key Pair
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
  VHilbert range;                                             //!< Key ranges for ranks
  VHilbert::iterator rangeBegin;                              //!< begin iterator of range
  VHilbert::iterator rangeEnd;                                //!< end iterator of range
  KeyPair globalKeyBounds;                                    //!< Global key bounds

private:
  //! Generate Space-filling order (Hilbert) from given the particles structure  and returns min/max orders
  KeyPair assignSFCtoBodies(Bodies& bodies, Bounds const& bounds, uint32_t& order) {
    const int accuracy = 10000;
    real_t const& _min = min(bounds.Xmin);                    // get min of all dimensions
    real_t const& _max = max(bounds.Xmax);                    // get max of all dimensions
    int64_t min_h = 0ull;                                     // initialize min Hilbert order
    int64_t max_h = 0ull;                                     // initialize max Hilbert order
    //order = cast_coord(ceil(log(accuracy) / log(2)));         // get needed bits to represent Hilbert order in each dimension
    order = 21;
    real_t diameter = _max - _min;                            // set domain's diameter
    assert(order <= 21);                                      // maximum key is 63 bits
    B_iter begin = bodies.begin();                            // bodies begin iterator
    for (int i = 0; i < bodies.size(); ++i) {                 // loop over bodies
      B_iter B = begin + i;                                   // Capture body iterator
      int position[3] = {                                     // initialize shifted position
        cast_coord((B->X[0] - _min) / diameter * accuracy),
        cast_coord((B->X[1] - _min) / diameter * accuracy),
        cast_coord((B->X[2] - _min) / diameter * accuracy)
      };
#if 1
      axesToTranspose(position, order);                       // get transposed Hilbert order
      B->ICELL = flattenTransposedKey(position, order);       // generate 1 flat Hilbert order (useful for efficient sorting)
#else
      B->ICELL = getHilbert(position, order);
#endif
      if (min_h != 0ull) {
        if (B->ICELL > max_h)                                 // update min/max hilbert orders
          max_h = B->ICELL;
        else if (B->ICELL < min_h)
          min_h = B->ICELL;
      }
      else {
        min_h = max_h = B->ICELL;
        assert(min_h != 0ull);
      }
    }                                                         // end loop
    return std::make_pair(min_h, max_h);                      // return min/max tuple
  }

  //! Gets the rank based on the key given the range
  inline int getPartitionNumber(int64_t const& key) {
    VHilbert::iterator low = std::lower_bound(rangeBegin, rangeEnd, key);//!< get iterator of first element that is >= input key
    if (low > rangeBegin) {
      return (low - rangeBegin) - 1;
    }
    return 0;
  }

  //! updates range parameters
  void updateRangeParameters(VHilbert const& rg) {
    range = rg;
    rangeBegin = range.begin();
    rangeEnd = range.end();
    globalKeyBounds.first = range[0];
    globalKeyBounds.second = range[mpisize];
  }

  //! applies partition sort algorithm to sort through particles based on Hilbert index
  void hilbertPartitionSort(VHilbert& keys, int64_t lbound, int64_t rbound, VHilbert&rq, int64_t pivot = 0, size_t depth = 0) {
    VHilbert right;
    VHilbert left;
    for (int i = 0; i < keys.size(); ++i) {
      if (keys[i] >= pivot)
        right.push_back(keys[i]);
      else left.push_back(keys[i]);
    }
    size_t sendSize[2];
    sendSize[0] = right.size();
    sendSize[1] = left.size();
    size_t size_d[2];
    MPI_Allreduce(sendSize, size_d, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); // Reduce domain SUM
    float diff = 0.0f;
    bool recurse = true;
    size_t max = std::max(size_d[0], size_d[1]);
    size_t min = std::min(size_d[0], size_d[1]);
    diff = (max - min) / float(max);
    size_t templ = lbound;
    size_t tempr = rbound;
    while (diff > 0.05f) {
      if (size_d[0] > size_d[1]) {
        templ = pivot;
        pivot += (tempr - pivot) >> 1;
      }
      else if (size_d[1] > size_d[0]) {
        tempr = pivot;
        pivot -= (pivot - templ) >> 1;
      }
      right.clear();
      left.clear();
      for (int i = 0; i < keys.size(); ++i) {
        if (keys[i] >= pivot)
          right.push_back(keys[i]);
        else left.push_back(keys[i]);
      }
      sendSize[0] = right.size();
      sendSize[1] = left.size();
      MPI_Allreduce(sendSize, size_d, 2, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); // Reduce domain SUM
      max = std::max(size_d[0], size_d[1]);
      min = std::min(size_d[0], size_d[1]);
      diff = (max - min) / float(max);
    }
    if (recurse) {
      rq.push_back(pivot);
      if (2 << depth < mpisize) {
        hilbertPartitionSort(left , lbound, pivot, rq, lbound + ((pivot - lbound) >> 1), depth + 1);
        hilbertPartitionSort(right, pivot, rbound, rq, pivot  + ((rbound - pivot) >> 1), depth + 1);
      }
    }
  }

public:
  //! Constructor
  Partition(int _mpirank, int _mpisize): mpirank(_mpirank), mpisize(_mpisize), numBins(16) {
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
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over MPI ranks
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
    for (int level = 0; level < numLevels; level++) {         // Loop over levels of MPI rank binary tree
      int numPartitions = rankColor[mpisize - 1] + 1;         //  Number of partitions in current level
      for (int ipart = 0; ipart < numPartitions; ipart++) {   //  Loop over partitions
        int irank = rankMap[ipart];                           //   Map partition to MPI rank
        Bounds bounds = rankBounds[irank];                    //   Bounds of current rank
        int direction = 0;                                    //   Initialize direction of partitioning
        real_t length = 0;                                    //   Initialize length of partition
        for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
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
        for (int b = bodyBegin; b < bodyEnd; b++) {           //    Loop over bodies in current partition
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
          for (int binRefine = 0; binRefine < 3; binRefine++) { //    Loop for bin refinement
            for (int ibin = 0; ibin < numBins; ibin++) {      //     Loop over bins
              countHist[ibin] = 0;                            //      Initialize body count histogram
              weightHist[ibin] = 0;                           //      Initialize body weight histogram
            }                                                 //     End loop over bins
            for (int b = bodyBegin; b < bodyEnd; b++) {       //     Loop over bodies
              real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
              int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
              if (ibin >= numBins) ibin = numBins - 1;        //      Account for floating point errors
              countHist[ibin]++;                              //      Increment body count histogram
              weightHist[ibin] += B[b].WEIGHT;                //      Increment body weight histogram
            }                                                 //     End loop over bodies
            scanHist[0] = countHist[0];                       //     Initialize scan array
            for (int ibin = 1; ibin < numBins; ibin++) {      //     Loop over bins
              scanHist[ibin] = scanHist[ibin - 1] + countHist[ibin]; //   Inclusive scan for body count histogram
            }                                                 //     End loop over bins
            for (int b = bodyEnd - 1; b >= bodyBegin; b--) {  //     Loop over bodies backwards
              real_t x = B[b].X[direction];                   //      Coordinate of body in current direction
              int ibin = (x - xmin + EPS) / (dx + EPS);       //      Assign bin index to body
              if (ibin >= numBins) ibin = numBins - 1;        //      Account for floating point errors
              scanHist[ibin]--;                               //      Decrement scanned histogram
              int bnew = scanHist[ibin] + bodyBegin;          //      New index of sorted body
              buffer[bnew] = B[b];                            //      Copy body to sort buffer
            }                                                 //     End loop over bodies backwards
            for (int b = bodyBegin; b < bodyEnd; b++) {       //     Loop over bodies
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
            for (int ibin = 1; ibin < numBins; ibin++) {      //     Loop over bins
              scanHist[ibin] = scanHist[ibin - 1] + countHist[ibin - 1]; // Exclusive scan of body count histogram
            }                                                 //     End loop over bins
            bodyBegin += scanHist[splitBin];                  //     Update body begin index
            bodyEnd = bodyBegin + countHist[splitBin];        //     Update body end index
          }                                                   //    End loop for bin refinement
        }                                                     //   End if for splitting partition
        int rankBegin = rankDispl[irank];                     //   Save current range of MPI ranks
        int rankEnd = rankBegin + rankCount[irank];           //   so that they don't get overwritten
        for (irank = rankBegin; irank < rankEnd; irank++) {   //   Loop over current range of MPI ranks
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
          if (level == numLevels - 1) rankColor[irank] = rankDispl[irank]; // Special case for final rank color
          rankKey[irank] = irank - rankDispl[irank];          //    Rank key is determined from rank displacement
        }                                                     //   End loop over current range of MPI ranks
      }                                                       //  End loop over partitions
      int ipart = 0;                                          //  Initialize partition index
      for (int irank = 0; irank < mpisize; irank++) {         //  Loop over MPI ranks
        if (rankKey[irank] == 0) {                            //   If rank key is zero it's a new group
          rankMap[ipart] = rankDispl[irank];                  //    Update rank map with rank displacement
          ipart++;                                            //    Increment partition index
        }                                                     //   End if for new group
      }                                                       //  End loop over MPI ranks
    }                                                         // End loop over levels of MPI rank binary tree
    B_iter B = bodies.begin();                                // Body begin iterator
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over MPI ranks
      int bodyBegin = sendDispl[irank];                       //  Body begin index for current rank
      int bodyEnd = bodyBegin + sendCount[irank];             //  Body end index for current rank
      for (int b = bodyBegin; b < bodyEnd; b++, B++) {        //  Loop over bodies in current rank
        B->IRANK = irank;                                     //   Copy MPI rank to body
      }                                                       //  End loop over bodies in current rank
    }                                                         // End loop over MPI ranks
    logger::stopTimer("Partition");                           // Stop timer
    return rankBounds[mpirank];                               // Return local bounds
  }

#if EXAFMM_COUNT_LIST
  template <typename T>
  void rebalance(Cells const& cells, Bodies& bodies, real_t const& numP2P, T& remoteInteractionList) {
    logger::startTimer("Partition");                           // Start timer
    for (int i = 0; i < bodies.size(); ++i) {
      bodies[i].IRANK = mpirank;
    }
    std::vector<double> remoteLoad(mpisize);
    const double imbalanceRate = 0.05;
    double sumWork = static_cast<double>(numP2P);
    MPI_Allgather(&sumWork, 1, MPI_DOUBLE, (double*)&remoteLoad[0], 1, MPI_DOUBLE, MPI_COMM_WORLD);// Gather all domain bounds
    double avgLoad = std::accumulate(remoteLoad.begin(), remoteLoad.end(), 0.0) / mpisize;
    double upperImbBound = avgLoad + avgLoad * imbalanceRate;
    double lowerImbBound = avgLoad - avgLoad * imbalanceRate;
    MPI_Comm balance_comm;
    int balance_rank, balance_size;
    if (remoteLoad[mpirank] > lowerImbBound)  {
      MPI_Comm_split(MPI_COMM_WORLD, 0, mpirank, &balance_comm);
      MPI_Comm_free(&balance_comm);
      logger::stopTimer("Partition");                           // Stop timer
      logger::startTimer("Sort");                               // Start timer
      logger::stopTimer("Sort");                                // Stop timer
      return;
    } else {
      MPI_Comm_split(MPI_COMM_WORLD, 1, mpirank, &balance_comm);
      MPI_Comm_rank(balance_comm, &balance_rank);
      MPI_Comm_size(balance_comm, &balance_size);
    }
    std::vector<bool> overloadedRanks(mpisize, false);
    int overloadedCount = 0;
    for (int irank = 1; irank < mpisize; ++irank) {
      if (remoteLoad[irank] > upperImbBound) {
        overloadedRanks[irank] = true;
        overloadedCount++;
      }
    }
    std::vector<bool> underloadedRanks(mpisize, false);
    for (int irank = 1; irank < mpisize; ++irank) {
      if (remoteLoad[irank] < lowerImbBound)
        underloadedRanks[irank] = true;
    }
    // the overloaded should send the underloaded ranks
    double loadDificit = remoteLoad[mpirank] - avgLoad;
    int diff = mpisize - mpirank;
    float** globalInteractionList = new float*[balance_size];
    float* pool = new float[balance_size * mpisize]; // allocate pool
    for (int irank = 0; irank < balance_size; ++irank, pool += mpisize) {
      globalInteractionList[irank] = pool;
    }
    float* interactionList = new float[mpisize];
    for (int i = 0; i < mpisize; ++i) interactionList[i] = 0;
    int cellCount = remoteInteractionList.size();
    for (int cc = 0; cc < cellCount; ++cc) {
      for (int irank = 0; irank < mpisize; ++irank) {
        if (overloadedRanks[irank]) {
          real_t localinteraction = remoteInteractionList[cc][mpirank];
          real_t remoteinteraction = remoteInteractionList[cc][irank];
          if (remoteinteraction > localinteraction)
            interactionList[irank] += remoteinteraction;
        }

      }
    }
    std::vector<int> ranksToBalance;
    MPI_Allgather(interactionList, mpisize, MPI_FLOAT, globalInteractionList[0], mpisize, MPI_FLOAT, balance_comm);// Gather interactions from balance communicators
    for (int irank = 0; irank < mpisize; ++irank) {
      if (overloadedRanks[irank]) {
        float max = 0;
        int maxRank = -1;
        for (int jrank = 0; jrank < balance_size; ++jrank) {
          if (globalInteractionList[jrank][irank] > max) {
            max = globalInteractionList[jrank][irank];
            maxRank = jrank;
          }
        }
        if (maxRank == balance_rank && maxRank >= 0)
          ranksToBalance.push_back(irank);
      }
    }
    B_iter B0 = bodies.begin();
    //std::cout<<"avgLoad: " << avgLoad << std::endl;
    int rebalanceSize = ranksToBalance.size();
    bool repartition = (rebalanceSize > 0);
    for (int i = 0; i < ranksToBalance.size(); ++i) {
      repartition = true;
      int irecv = ranksToBalance[i];
      int interactionCount = remoteInteractionList.size();
      double rightLoad = remoteLoad[irecv];
      //std::cout<<" rank " << mpirank << " is balancing " << irecv <<" with load " << rightLoad << std::endl;
      for (int cc = 0; cc < interactionCount; ++cc) {
        std::vector<int> const& cellInteraction = remoteInteractionList[cc];
        int maxDiff = 0;
        int index = 0;
        if (cellInteraction[irecv] > cellInteraction[mpirank]) {
          double cellLoad = std::accumulate(cellInteraction.begin(), cellInteraction.end(), 0.0);
          loadDificit += cellLoad - cellInteraction[mpirank];
          rightLoad -= cellLoad - cellInteraction[irecv];
          if (loadDificit > 0 || rightLoad <= avgLoad)
            break;
          int ibody = cells[cc].IBODY;
          int nbody = cells[cc].NBODY;
          for (int ii = ibody; ii < nbody; ++ii) {
            (B0 + ii)->IRANK = irecv;
          }
          remoteInteractionList.erase(remoteInteractionList.begin() + cc);
          interactionCount--;
        }
      }
      //std::cout<<" rank " << mpirank << " reduced load of " << irecv <<" to  " << rightLoad << std::endl;
    }
    logger::stopTimer("Partition");                           // Start timer
    logger::startTimer("Sort");                               // Start timer
    if (repartition) {
      Sort sort;                                              // Instantiate sort class
      bodies = sort.irank(bodies);                            // Sort bodies according to IRANK
    }
    logger::stopTimer("Sort");                                // Stop timer
    MPI_Comm_free(&balance_comm);
    delete[] globalInteractionList[0];
    delete[] globalInteractionList;
    delete[] interactionList;
  }
#endif
  //! Partition bodies with geometric octsection
  Bounds octsection(Bodies & bodies, Bounds global) {
    logger::startTimer("Partition");                          // Start timer
    int size = mpisize;                                       // Initialize MPI size counter
    int Npartition[3] = {1, 1, 1};                            // Number of partitions in each direction
    int d = 0;                                                // Initialize dimension counter
    while (size != 1) {                                       // Divide domain while counter is not one
      Npartition[d] <<= 1;                                    //  Divide this dimension
      d = (d + 1) % 3;                                        //  Increment dimension
      size >>= 1;                                             //  Right shift the bits of counter
    }                                                         // End while loop for domain subdivision
    real_t Xpartition[3];                                     // Size of partitions in each direction
    for (d = 0; d < 3; d++) {                                 // Loop over dimensions
      Xpartition[d] = (global.Xmax[d] - global.Xmin[d]) / Npartition[d];//  Size of partition in each direction
    }                                                         // End loop over dimensions
    int iX[3];                                                // Index vector
    iX[0] = mpirank % Npartition[0];                          // x index of partition
    iX[1] = mpirank / Npartition[0] % Npartition[1];          // y index
    iX[2] = mpirank / Npartition[0] / Npartition[1];          // z index
    Bounds local;                                             // Local bounds
    for (d = 0; d < 3; d++) {                                 // Loop over dimensions
      local.Xmin[d] = global.Xmin[d] + iX[d] * Xpartition[d]; // Xmin of local domain at current rank
      local.Xmax[d] = global.Xmin[d] + (iX[d] + 1) * Xpartition[d];// Xmax of local domain at current rank
    }                                                         // End loop over dimensions
    for (B_iter B = bodies.begin(); B != bodies.end(); B++) { // Loop over bodies
      int ic = 0;                                             //  Residual index
      if (B->ICELL < 0) ic = B->ICELL;                        //  Use first body in group
      for (d = 0; d < 3; d++) {                               //  Loop over dimensions
        iX[d] = int(((B + ic)->X[d] - global.Xmin[d]) / Xpartition[d]); //   Index vector of partition
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

  void partitionHilbert(Bodies & bodies, Bounds globalBounds) {
    logger::startTimer("Partition");                          // Start timer
    if (mpisize > 1) {
      uint32_t depth;
      KeyPair localHilbertBounds = assignSFCtoBodies(bodies, globalBounds, depth);
      int64_t min, max;
      MPI_Allreduce(&localHilbertBounds.first,  &min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);// Reduce domain Xmin
      MPI_Allreduce(&localHilbertBounds.second, &max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmax
      globalKeyBounds = std::make_pair(min, max);
      int64_t lbound = globalKeyBounds.first;
      int64_t rbound = globalKeyBounds.second;
      VHilbert rq;
      rq.reserve(mpisize + 1);
      rq.push_back(lbound); rq.push_back(rbound);
      int64_t startPivot = lbound + ((rbound - lbound) >> 1);
      VHilbert keys(bodies.size());
      int index = 0;
      for (B_iter B = bodies.begin(); B != bodies.end(); ++B)
        keys[index++] = B->ICELL;
      hilbertPartitionSort(keys, lbound, rbound, rq, startPivot);
      std::sort(rq.begin(), rq.end());
      updateRangeParameters(rq);
      for (B_iter B = bodies.begin(); B != bodies.end(); ++B) {
        B->IRANK = getPartitionNumber(B->ICELL);                  // Get partition
        assert(B->IRANK >= 0 && B->IRANK < mpisize);              // Make sure rank is within communication size
      }
    } else {
      for (B_iter B = bodies.begin(); B != bodies.end(); ++B)
        B->IRANK = mpirank;
    }
    logger::stopTimer("Partition");                             // Stop timer
    logger::startTimer("Sort");
    Sort sort;
    bodies = sort.irank(bodies);
    logger::stopTimer("Sort");
  }

  //! Wrapper for partitioning
  Bounds partition(Bodies& bodies, Bounds globalBounds, const char* partitioingType) {
    Bounds localBounds;
    switch (partitioingType[0]) {
    case 'h':
      partitionHilbert(bodies, globalBounds);
      break;
    case 'b':
      localBounds = bisection(bodies, globalBounds);
      break;
    case 'o':
      localBounds = octsection(bodies, globalBounds);
      break;
    default:
      localBounds = bisection(bodies, globalBounds);
    }
    return localBounds;
  }

  //! Send bodies back to where they came from
  void unpartition(Bodies & bodies) {
    logger::startTimer("Sort");                                 // Start timer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.irank(bodies);                                // Sort bodies according to IRANK
    logger::stopTimer("Sort");                                  // Stop timer
  }
};
}
#endif
