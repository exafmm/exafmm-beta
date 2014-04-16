#ifndef partition_h
#define partition_h
#include <algorithm>
#include "logger.h"
#include "sort.h"

//! Handles all the partitioning of domains
class Partition {
private:
  typedef std::vector<int> Ints;                                //!< Vector of int
  const int mpirank;                                            //!< Rank of MPI communicator
  const int mpisize;                                            //!< Size of MPI communicator
  int mpirankSplit;                                             //!< Rank of split MPI communicator
  int mpisizeSplit;                                             //!< Size of split MPI communicator
  int numBins;                                                  //!< Number of cells in one dimension (leaf level)
  int numLevels;                                                //!< Level of the MPI process binary tree
  int rankCount[64][2];                                         //!< Number of processes in the two split groups
  int rankDispl[64][2];                                         //!< Offset of body in the two split groups
  int rankColor[64][3];                                         //!< Color of hypercube communicators
  int   rankKey[64][3];                                         //!< RankKey of hypercube communicators
  MPI_Comm MPI_COMM_SPLIT[64][3];                               //!< Hypercube communicators
  Bounds * bounds;                                              //!< Array of bounds of domain at each level
  Bodies buffer;                                                //!< MPI communication buffer for bodies

private:
  //! Split range and return partial range
  void splitRange(int & begin, int & end, int iSplit, int numSplit) {
    assert(end > begin);                                        // Check that size > 0
    int size = end - begin;                                     // Size of range
    int increment = size / numSplit;                            // Increment of splitting
    int remainder = size % numSplit;                            // Remainder of splitting
    begin += iSplit * increment + std::min(iSplit,remainder);   // Increment the begin counter
    end = begin + increment;                                    // Increment the end counter
    if (remainder > iSplit) end++;                              // Adjust the end counter for remainder
  }

  //! Max level for bottom up tree build
  int getMaxLevel(Bodies &bodies) {
    const int NCRIT = 16;                                       // Fake NCRIT
    const long N = bodies.size() * mpisize;                     // Number of bodies
    int level;                                                  // Max level
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;   // Decide max level from N/Ncrit
    int MPIlevel = int(log(mpisize-1) / M_LN2 / 3) + 1;         // Level of local root cell
    if( mpisize == 1 ) MPIlevel = 0;                            // For serial execution local root cell is root cell
    if( MPIlevel > level ) {                                    // If process hierarchy is deeper than tree
      level = MPIlevel;                                         //  Increase level to process hierarchy depth
    }                                                           // End if for process hierarchy deeper than tree
    return level;                                               // Return max level
  }

//! Turn positions into indices of bins
  void binBodies(Bodies & bodies, Bounds globalBounds, int d) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->IRANK = int((B->X[d] - globalBounds.Xmin[d])           // Bin body positions into integers
        / (globalBounds.Xmax[d] - globalBounds.Xmin[d]) * numBins);
    }                                                           // End loop over bodies
  }

//! Split bodies according to iSplit
  int splitBodies(Bodies & bodies, int iSplit) {
    int nth = 0;                                                // Initialize splitting index
    while (int(bodies[nth].IRANK) <= iSplit && nth < int(bodies.size())) nth++;// Determine end index of left side
    return nth;                                                 // Return end index of left side
  }

//! Split domain according to iSplit
  void splitDomain(int iSplit, int l, int d) {
    real_t X = (iSplit + 1) * (bounds[l].Xmax[d] - bounds[l].Xmin[d]) / numBins + bounds[l].Xmin[d];// Coordinate corresponding to iSplit
    bounds[l+1].Xmax = bounds[l].Xmax;                          // Set Xmax for next subdivision
    bounds[l+1].Xmin = bounds[l].Xmin;                          // Set Xmin for next subdivision
    if ((rankColor[l+1][0] & 1) == 0) {                         // If on left side
      bounds[l+1].Xmax[d] = X;                                  //  Set max to X
    } else {                                                    // If on right side
      bounds[l+1].Xmin[d] = X;                                  //  Set min to X
    }                                                           // End if for sides
  }

//! Get global bucket bodies for parallel nth_element
  int getBucket(Bodies & bodies, int numBodies, int lOffset, Ints & send, Ints & recv, MPI_Comm MPI_COMM) {
    int maxBucket = send.size();                                // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample = std::min(maxBucket/mpisizeSplit, numBodies);// Number of local samples
    int * recvCount = new int [mpisizeSplit];                   // MPI recv count
    int * recvDispl = new int [mpisizeSplit];                   // MPI recv displacement
    for (int i=0; i<numSample; i++) {                           // Loop over local samples
      int stride = numBodies / numSample;                       //  Sampling stride
      send[i] = bodies[lOffset+i*stride].IRANK;                 //  Put sampled bodies in send buffer
    }                                                           // End loop over samples
    MPI_Gather(&numSample, 1, MPI_INT, recvCount, 1, MPI_INT, 0, MPI_COMM); // Gather size of sample bodies to rank 0
    if (mpirankSplit == 0) {                                    // Only rank 0 operates on gathered info
      numBucket = 0;                                            //  Initialize number of buckets
      for (int irank=0; irank<mpisizeSplit; irank++) {          //  Loop over processes
        recvDispl[irank] = numBucket;                           //   Calculate recv displacement
        numBucket += recvCount[irank];                          //   Accumulate recv count to get number of buckets
      }                                                         //  End loop over processes
      recv.resize(numBucket);                                   //  Resize recv so that end() is valid
    }                                                           // End if for rank 0
    MPI_Gatherv(&send[0], numSample, MPI_INT,                   // Gather sample bodies to rank 0
                &recv[0], recvCount, recvDispl, MPI_INT, 0, MPI_COMM);
    if (mpirankSplit == 0) {                                    // Only rank 0 operates on gathered info
      std::sort(recv.begin(), recv.end());                      //  Sort the bucket bodies
      numBucket = std::unique(recv.begin(), recv.end())-recv.begin();// Remove duplicate bucket bodies
      recv.resize(numBucket);                                   //  Resize recv again
    }                                                           // End if for rank 0
    MPI_Bcast(&numBucket, 1, MPI_INT, 0, MPI_COMM);             // Broadcast number of buckets
    MPI_Bcast(&recv[0], numBucket, MPI_INT, 0, MPI_COMM);       // Broadcast bucket bodies
    delete[] recvCount;                                         // Delete recv count
    delete[] recvDispl;                                         // Delete recv displacement
    return numBucket;                                           // Return number of buckets
  }

//! Parallel global nth_element on distributed memory
  int nth_element(Bodies & bodies, int n, MPI_Comm MPI_COMM=0) {
    if( MPI_COMM == 0 ) {                                       // If MPI_COMM is not specified
      MPI_Comm_split(MPI_COMM_WORLD, 0, mpirank, &MPI_COMM);    //  Create an artificial MPI_COMM
    }
    MPI_Comm_size(MPI_COMM, &mpisizeSplit);                     // Get number of MPI processes for split comm
    MPI_Comm_rank(MPI_COMM, &mpirankSplit);                     // Get index of current MPI process for split comm
    int numBodies = bodies.size();                              // Total size of bodies to perform nth_element
    int maxBucket = 1000;                                       // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int lOffset = 0;                                            // Local offset of region being considered
    int gOffset = 0;                                            // Global offset of region being considered
    int * recvCount = new int [mpisizeSplit];                   // MPI recv count
    int * isend = new int [maxBucket];                          // MPI send buffer for index
    int * irecv = new int [maxBucket];                          // MPI recv buffer for index
    int * iredu = new int [maxBucket];                          // Local scan of index buffer
    Ints send(maxBucket);                                       // MPI send buffer for bodies
    Ints recv(maxBucket);                                       // MPI recv buffer for bodies
    numBucket = getBucket(bodies, numBodies, lOffset, send, recv, MPI_COMM);// Get global bucket bodies

    while (numBucket > 1) {                                     // While there are multipole candidates
      int ic=0, nth=0;                                          //  Initialize counters
      for (int i=0; i<maxBucket; i++) isend[i] = 0;             //  Initialize bucket counter
      for (int i=0; i<numBodies; i++) {                         //  Loop over range of bodies
        while (int(bodies[lOffset + i].IRANK) > recv[ic] && ic < numBucket-1) ic++;// Set counter to current bucket
        isend[ic]++;                                            //   Increment bucket counter
      }                                                         //  End loop over bodies
      MPI_Reduce(isend, irecv, numBucket, MPI_INT, MPI_SUM, 0, MPI_COMM); // Reduce bucket counter
      if (mpirankSplit == 0) {                                  //  Only rank 0 operates on reduced bodies
        iredu[0] = 0;                                           //   Initialize global scan index
        for (int i=0; i<numBucket-1; i++) {                     //   Loop over buckets
          iredu[i+1] = iredu[i] + irecv[i];                     //    Increment global scan index
        }                                                       //   End loop over buckets
        nth = 0;                                                //   Initialize index for bucket containing nth element
        while (n-gOffset > iredu[nth] && nth < numBucket) nth++;//   Find index for bucket containing nth element
        --nth;                                                  //   Account for overshoot (do while?)
        if (nth == -1) nth = 0;                                 //   If nth is -1 don't split
        gOffset += iredu[nth];                                  //   Increment global offset of region being considered
      }                                                         //  End if for rank 0
      MPI_Bcast(&nth, 1, MPI_INT, 0, MPI_COMM);                 //  Broadcast index for bucket containing nth element
      MPI_Bcast(&gOffset, 1, MPI_INT, 0, MPI_COMM);             //  Broadcast global offset
      iredu[0] = 0;                                             //  Initialize local scan index
      for (int i=0; i<numBucket-1; i++) {                       //  Loop over buckets
        iredu[i+1] = iredu[i] + isend[i];                       //   Increment local scan index
      }                                                         //  End loop over buckets
      if (nth == numBucket-1) {                                 //  If nth is last bucket
        numBodies = numBodies-iredu[nth];                       //   Region of interest is that bucket
      } else {                                                  //  If nth is not the last bucket
        numBodies = iredu[nth+1]-iredu[nth];                    //   Region of interest is that bucket
      }                                                         //  End if for last bucket
      lOffset += iredu[nth];                                    //  Increment local offset to that bucket
      numBucket = getBucket(bodies, numBodies, lOffset, send, recv, MPI_COMM);// Get global bucket bodies
    }                                                           // End while loop
    delete[] recvCount;                                         // Delete recv count
    delete[] isend;                                             // Delete send buffer for index
    delete[] irecv;                                             // Delete recv buffer for index
    delete[] iredu;                                             // Delete local scan of index buffer
    return recv[0]-1;                                           // Return nth element
  }

//! Split the MPI communicator into N-D hypercube
  void bisectionGetComm(int l) {
    for (int i=1; i>=0; i--) {                                  // Loop over 1 and 0
      int pSplit = (rankCount[l][0] + i) / 2;                   //  Splitting criteria
      if (mpirank - rankDispl[l][0] < pSplit) {                 //  If on left side
        rankCount[l+1][i] = pSplit;                             //   Group size is the splitting criteria
        rankDispl[l+1][i] = rankDispl[l][0];                    //   Offset is the same as previous
         rankColor[l+1][i] =  rankColor[l][0] * 2;              //   Color is twice the previous
      } else {                                                  //  If on right side
        rankCount[l+1][i] = rankCount[l][0] - pSplit;           //   Group size is what remains from the left side
        rankDispl[l+1][i] = rankDispl[l][0] + pSplit;           //   Offset is incremented by splitting criteria
         rankColor[l+1][i] =  rankColor[l][0] * 2 + 1;          //   Color is twice the previous plus one
      }                                                         //  End if for sides
      rankKey[l+1][i] = mpirank - rankDispl[l+1][i];            //  Key is determined from offset
    }                                                           // End loop over 1 and 0
    rankKey[l+1][2] = rankColor[l+1][1] & 1;                    // The third type of key is determined from color[1]
    rankColor[l+1][2] = rankKey[l+1][1] + rankColor[l][0] * (1 << (numLevels - l - 1));// The corresponding color is determined from key[1]
    for (int i=0; i<3; i++) {                                   // Loop over the three types of colors and keys
      MPI_Comm_split(MPI_COMM_WORLD, rankColor[l+1][i],         // Split the MPI communicator
		     rankKey[l+1][i], &MPI_COMM_SPLIT[l+1][i]);
    }                                                           // End loop over three types
#ifdef DEBUG
    BaseMPI baseMPI;
    baseMPI.print("level : ",0);                                // Print identifier
    baseMPI.print(l+1,0);                                       // Print current level
    baseMPI.print("\n",0);                                      // New line
    baseMPI.print("key   : \n",0);                              // Print identifier
    baseMPI.print(rankKey[l+1],0,3);                            // Print rankKey
    baseMPI.print("color : \n",0);                              // Print identifier
    baseMPI.print(rankColor[l+1],0,3);                          // Print rankColor
#endif
  }

//! One-to-one MPI_Alltoallv
  void bisectionAlltoall(Bodies & bodies, int nthLocal, int numLocal, int & newSize, int l) {
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int sendCount[2] = {nthLocal, numLocal - nthLocal};         // Set send count to right and left size
    int recvCount[2] = {0, 0};                                  // Initialize recv count
    MPI_Alltoall(sendCount, 1, MPI_INT,                         // Communicate the send count to get recv count
		 recvCount, 1, MPI_INT, MPI_COMM_SPLIT[l+1][2]);
    int sendDispl[2] = {0, sendCount[0]};                       // Send displacement
    int recvDispl[2] = {0, recvCount[0]};                       // Recv displacement
    newSize = recvCount[0] + recvCount[1];                      // Get new size from recv count
    if (rankColor[l+1][0] != rankColor[l+1][1]) newSize = numLocal; // Unless it's a leftover process of an odd group
    buffer.resize(newSize);                                     // Resize recv buffer

    for (int i=0; i<2; i++) {                                   // Loop over 0 and 1
      sendCount[i] *= bytes;                                    // Multiply send count by byte size of data
      sendDispl[i] *= bytes;                                    // Multiply send displacement by byte size of data
      recvCount[i] *= bytes;                                    // Multiply recv count by byte size of data
      recvDispl[i] *= bytes;                                    // Multiply recv displacement by byte size of data
    }                                                           // End loop over 0 and 1
    MPI_Alltoallv(&bodies[0], sendCount, sendDispl, MPI_BYTE,   // Communicate bodies
                  &buffer[0], recvCount, recvDispl, MPI_BYTE,   // Using same buffer as sort buffer
                  MPI_COMM_SPLIT[l+1][2]);                      // MPI_COMM_SPLIT[2] is for the one-to-one pair
    if (rankColor[l+1][0] == rankColor[l+1][1]) bodies = buffer;// Don't update if leftover process
    buffer.resize(bodies.size());                               // Resize sort buffer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.irank(bodies);                                // Sort bodies in ascending order
  }

//! Scattering from leftover processes
  void bisectionScatter(Bodies & bodies, int nthLocal, int & newSize, int l) {
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numScatter = rankCount[l+1][1] - 1;                     // Number of processes to scatter to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int * sendCount = new int [rankCount[l+1][1]];              // Send count
    int * sendDispl = new int [rankCount[l+1][1]];              // Send displacement
    int recvCount;                                              // Recv count
    if (rankKey[l+1][1] == numScatter) {                        // If this is the leftover proc to scatter from
      sendDispl[0] = 0;                                         //  Initialize send displacement
      for (int i=0; i<numScatter; i++) {                        //  Loop over processes to send to
        int begin = 0, end = nthLocal;                          //  Set begin and end of range to send
        splitRange(begin, end, i, numScatter);                  //  Split range into numScatter and get i-th range
        sendCount[i] = end - begin;                             //  Set send count based on range
        sendDispl[i+1] = sendDispl[i] + sendCount[i];           //  Set send displacement based on send count
      }                                                         // End loop over processes to send to
      sendCount[numScatter] = 0;                                // Send count to self should be 0
      oldSize = 0;                                              // Reset oldSize to account for sent bodies
      newSize -= sendDispl[numScatter];                         // Set newSize to account for sent bodies
      buffer.erase(buffer.begin(), buffer.begin()+sendDispl[numScatter]);// Erase from recv buffer the part that was sent
    }                                                           // End if for leftover proc
    MPI_Scatter(sendCount,  1, MPI_INT,                         // Scatter the send count to get recv count
		&recvCount, 1, MPI_INT,numScatter, MPI_COMM_SPLIT[l+1][1]);
    if (rankKey[l+1][1] != numScatter) {                        // If it's one the receiving ranks
      newSize += recvCount;                                     //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize the recv buffer
    }                                                           // End if for receiving ranks
    for (int i=0; i<rankCount[l+1][1]; i++) {                   // Loop over group of processes
      sendCount[i] *= bytes;                                    //  Multiply send count by byte size of data
      sendDispl[i] *= bytes;                                    //  Multiply send displacement by byte size of data
    }                                                           // End loop over group of processes
    recvCount *= bytes;                                         // Multiply recv count by byte size of data
    MPI_Scatterv(&bodies[0], sendCount, sendDispl, MPI_BYTE,    // Communicate bodies via MPI_Scatterv
                 &buffer[oldSize], recvCount, MPI_BYTE,         // Offset recv buffer by oldSize
                 numScatter, MPI_COMM_SPLIT[l+1][1]);           // MPI_COMM_SPLIT[1] is used for scatter
    bodies = buffer;                                            // Copy recv buffer to bodies
    if (rankKey[l+1][1] != numScatter) {                        // If it's one the receiving ranks
      Sort sort;                                                //  Instantiate sort class
      bodies = sort.irank(bodies);                              //  Sort bodies in ascending order
    }                                                           // End if for receiving ranks
    delete[] sendCount;                                         // Delete send count
    delete[] sendDispl;                                         // Delete send displacement
  }

//! Gathering to leftover processes
  void bisectionGather(Bodies & bodies, int nthLocal, int numLocal, int & newSize, int l) {
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numGather = rankCount[l+1][0] - 1;                      // Number of processes to gather to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int sendCount;                                              // Send count
    int * recvCount = new int [rankCount[l+1][0]];              // Recv count
    int * recvDispl = new int [rankCount[l+1][0]];              // Recv displacement
    if (rankKey[l+1][0] > 0) {                                  // If this is not the leftover proc
      int begin = 0, end = numLocal - nthLocal;                 //  Set begin and end of range to send
      splitRange(begin, end, rankKey[l+1][0]-1, rankCount[l+1][0]); // Split range into rankCount[0]
      sendCount = end - begin;                                  //  Set send count based on range
      newSize -= sendCount;                                     //  Set newSize to account for sent data
      buffer.erase(buffer.begin(), buffer.begin()+sendCount);   //  Erase from recv buffer the part that was sent
    }                                                           // End if for leftover proc
    MPI_Barrier(MPI_COMM_SPLIT[l+1][0]);                        // Sync processes
    MPI_Gather(&sendCount, 1, MPI_INT,                          // Gather the send count to get recv count
	       recvCount,  1, MPI_INT, 0, MPI_COMM_SPLIT[l+1][0]);
    if (rankKey[l+1][0] == 0) {                                 // If this is the leftover proc
      recvDispl[0] = 0;                                         //  Initialize the recv displacement
      for (int i=0; i<numGather; i++) {                         //  Loop over processes to gather from
        recvDispl[i+1] = recvDispl[i] + recvCount[i];           //   Set recv displacement based on recv count
      }                                                         //  End loop over processes
      newSize += recvDispl[numGather] + recvCount[numGather];   //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize recv buffer
    }                                                           // End if for leftover proc
    sendCount *= bytes;                                         // Multiply send count by byte size of data
    for (int i=0; i<rankCount[l+1][0]; i++) {                   // Loop over group of processes
      recvCount[i] *= bytes;                                    //  Multiply recv count by byte size of data
      recvDispl[i] *= bytes;                                    //  Multiply recv displacement by byte size of data
    }                                                           // End loop over group of processes
    MPI_Gatherv(&bodies[0], sendCount, MPI_BYTE,                // Communicate bodies via MPI_Gatherv
                &buffer[oldSize], recvCount, recvDispl, MPI_BYTE, // Offset recv buffer by oldSize
                0, MPI_COMM_SPLIT[l+1][0]);                     // MPI_COMM_SPLIT[0] is used for gather
    bodies = buffer;                                            // Copy recv buffer to bodies
    delete[] recvCount;                                         // Delete recv count
    delete[] recvDispl;                                         // Delete send count
    if (rankKey[l+1][0] == 0) {                                 // If this is the leftover proc
      Sort sort;                                                //  Instantiate sort class
      bodies = sort.irank(bodies);                              //  Sort bodies in ascending order
    }                                                           // End if for leftover proc
  }

public:
  //! Constructor
  Partition(int _mpirank, int _mpisize) : mpirank(_mpirank), mpisize(_mpisize) {
    numLevels = int(log(mpisize) / M_LN2 - 1e-5) + 1;           // Level of the process binary tree
    if (mpisize == 1) numLevels = 0;                            // Level is 0 for a serial execution
    bounds = new Bounds [numLevels+2];                          // Allocate bounds array
    logger::verbose = mpirank == 0;                             // Print the split comm time
    logger::startTimer("Split comm");                           // Start timer
    rankCount[0][0] = rankCount[0][1] = mpisize;                // Initialize number of processes in groups
    rankDispl[0][0] = rankDispl[0][1] = 0;                      // Initialize offset of body in groups
    rankColor[0][0] = rankColor[0][1] = rankColor[0][2] = 0;    // Initialize color of communicators
      rankKey[0][0] =   rankKey[0][1] =   rankKey[0][2] = 0;    // Initialize key of communicators
    for (int l=0; l<numLevels; l++) {                           // Loop over levels of N-D hypercube communication
      bisectionGetComm(l);                                      //  Split the MPI communicator for that level
    }                                                           // End loop over levels of N-D hypercube communication
    logger::stopTimer("Split comm");                            // Stop timer
  }

//! Destructor
  ~Partition() {
    delete[] bounds;
  }

  //! Partitioning by recursive bisection
  Bounds bisection(Bodies & bodies, Bounds global) {
    logger::startTimer("Partition");                            // Start timer
    int newSize;                                                // New size of recv buffer
    int numLocal = bodies.size();                               // Local bodies size
    int numGlobal;                                              // Global bodies size
    MPI_Allreduce(&numLocal, &numGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);// Reduce Global bodies size
    int nthGlobal = (numGlobal * (rankCount[0][0] / 2)) / rankCount[0][0];// Split at nth global element
    numBins = 1 << getMaxLevel(bodies);                         // Get number of cells in 1D
    bounds[0].Xmin = global.Xmin;                               // Copy global bounds Xmin
    bounds[0].Xmax = global.Xmax;                               // Copy global bounds Xmax
    binBodies(bodies, bounds[0], 2);                            // Bin bodies into leaf level cells
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.irank(bodies);                                // Sort bodies in ascending order
    int iSplit = nth_element(bodies, nthGlobal);                // Get cell index of nth global element
    int nthLocal = splitBodies(bodies, iSplit);                 // Split bodies based on iSplit
    for (int l=0; l<numLevels; l++) {                           // Loop over levels of N-D hypercube communication
      splitDomain(iSplit, l, 2-l%3);                            //  Split the domain according to iSplit
      bisectionAlltoall(bodies, nthLocal, numLocal, newSize, l);//  Communicate bodies by one-to-one MPI_Alltoallv
      if ((rankCount[l][0] & 1) == 1 && rankCount[l][0] != 1 && rankCount[l+1][0] <= rankCount[l+1][1])// If scatter is necessary
        bisectionScatter(bodies, nthLocal, newSize, l);         //  Communicate bodies by scattering from leftover proc
      if ((rankCount[l][0] & 1) == 1 && rankCount[l][0] != 1 && rankCount[l+1][0] >= rankCount[l+1][1])// If gather is necessary
        bisectionGather(bodies, nthLocal, numLocal, newSize, l);//  Communicate bodies by gathering to leftover proc
      numLocal = newSize;                                       //  Update local bodies size
      MPI_Allreduce(&numLocal, &numGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_SPLIT[l+1][0]);// Reduce global bodies size
      nthGlobal = (numGlobal * (rankCount[l+1][0] / 2)) / rankCount[l+1][0];//  Split at nth global element
      binBodies(bodies, bounds[0], 2-(l+1)%3);                  //  Bin bodies into leaf level cells
      bodies = sort.irank(bodies);                              //  Sort bodies in ascending order
      iSplit = nth_element(bodies, nthGlobal, MPI_COMM_SPLIT[l+1][0]); // Get cell index of nth global element
      nthLocal = splitBodies(bodies, iSplit);                   //  Split bodies based on iSplit
    }                                                           // End loop over levels of N-D hypercube communication
    logger::stopTimer("Partition");                             //  Stop timer
    return bounds[numLevels-1];                                 // Return local bounds
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
