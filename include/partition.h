#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "construct.h"

namespace {
class Partition : public MyMPI, public TreeConstructor {        // Handles all the partitioning of domains
private:
  int numCells1D;                                               // Number of cells in one dimension (leaf level)

protected:
  int LEVEL;                                                    // Level of the MPI process binary tree
  std::vector<vect> XMIN;                                       // Minimum position vector of bodies
  std::vector<vect> XMAX;                                       // Maximum position vector of bodies
  int oldnprocs;                                                // Old number of processes in the group
  int nprocs[2];                                                // Number of processes in the two split groups
  int offset[2];                                                // Offset of body in the two split groups
  int  color[3];                                                // Color of Gather, Scatter, and Alltoall communicators
  int    key[3];                                                // Key of Gather, Scatter, and Alltoall communicators
  MPI_Comm MPI_COMM[3];                                         // Communicators for Gather, Scatter, and Alltoall

private:
  void splitDomain(bigint iSplit, int l, int d) {               // Split domain accoring to iSplit
    real X = (iSplit + 1) * (XMAX[l][d] - XMIN[l][d]) / numCells1D + XMIN[l][d];// Coordinate corresponding to iSplit
    XMAX[l+1] = XMAX[l];                                        // Set XMAX for next subdivision
    XMIN[l+1] = XMIN[l];                                        // Set XMIN for next subdivision
    if( color[0] % 2 == 0 ) {                                   // If on left side
      XMAX[l+1][d] = X;                                         //  Set max to X
    } else {                                                    // If on right side
      XMIN[l+1][d] = X;                                         //  Set min to X
    }                                                           // Endif for sides
  }

  template<typename T>
  int getBucket(T &data, int numData, int lOffset, Bigints &send, Bigints &recv, MPI_Comm MPI_COMM0) {
    int maxBucket = send.size();                                // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample = std::min(maxBucket/SIZES,numData);          // Number of local samples
    MPI_Datatype MPI_TYPE = getType(data[0].I);                 // Get MPI data type
    int *rcnt = new int [SIZES];                                // MPI recv count
    int *rdsp = new int [SIZES];                                // MPI recv displacement
    for( int i=0; i!=numSample; ++i ) {                         // Loop over local samples
      int stride = numData/numSample;                           //  Sampling stride
      send[i] = data[lOffset + i * stride].I;                   //  Put sampled data in send buffer
    }                                                           // End loop over samples
    MPI_Gather(&numSample,1,MPI_INT,                            // Gather size of sample data to rank 0
               rcnt,      1,MPI_INT,
               0,MPI_COMM0);
    if( RANKS == 0 ) {                                          // Only rank 0 operates on gathered info
      numBucket = 0;                                            //  Initialize number of buckets
      for( int irank=0; irank!=SIZES; ++irank ) {               //  Loop over processes
        rdsp[irank] = numBucket;                                //   Calculate recv displacement
        numBucket += rcnt[irank];                               //   Accumulate recv count to get number of buckets
      }                                                         //  End loop over processes
      recv.resize(numBucket);                                   //  Resize recv so that end() is valid
    }                                                           // Endif for rank 0
    MPI_Gatherv(&send[0],numSample,MPI_TYPE,                    // Gather sample data to rank 0
                &recv[0],rcnt,rdsp,MPI_TYPE,
                0,MPI_COMM0);
    if( RANKS == 0 ) {                                          // Only rank 0 operates on gathered info
      std::sort(recv.begin(),recv.end());                       //  Sort the bucket data
      numBucket = std::unique(recv.begin(),recv.end()) - recv.begin();// Remove duplicate bucket data
      recv.resize(numBucket);                                   //  Resize recv again
    }                                                           // Endif for rank 0
    MPI_Bcast(&numBucket,1,MPI_INT,0,MPI_COMM0);                // Broadcast number of buckets
    MPI_Bcast(&recv[0],numBucket,MPI_TYPE,0,MPI_COMM0);         // Broadcast bucket data
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    return numBucket;                                           // Return number of buckets
  }

protected:
  void bisectionGetComm(int l) {                                // Split the MPI communicator into N-D hypercube
    int  oldcolor = color[0];                                   // Save color of previous level
    oldnprocs = nprocs[0];                                      // Save group size of previous level
    for( int i=1; i>=0; --i ) {                                 // Loop over 1 and 0
      int pSplit = (oldnprocs + i) / 2;                         //  Splitting criteria
      if( RANK - offset[0] < pSplit ) {                         //  If on left side
        nprocs[i] = pSplit;                                     //   Group size is the splitting criteria
        offset[i] = offset[0];                                  //   Offset is the same as previous
         color[i] =  color[0] * 2;                              //   Color is twice the previous
      } else {                                                  //  If on right side
        nprocs[i] = oldnprocs - pSplit;                         //   Group size is what remains from the left side
        offset[i] = offset[0] + pSplit;                         //   Offset is incremented by splitting criteria
         color[i] =  color[0] * 2 + 1;                          //   Color is twice the previous plus one
      }                                                         //  Endif for sides
      key[i] = RANK - offset[i];                                //  Key is determined from offset
    }                                                           // End loop over 1 and 0
    key[2] = color[1] % 2;                                      // The third type of key is determined from color[1]
    color[2] = key[1] + oldcolor * (1 << (LEVEL - l - 1));      // The corresponding color is determined from key[1]
    for( int i=0; i!=3; ++i ) {                                 // Loop over the three types of colors and keys
      MPI_Comm_split(MPI_COMM_WORLD,color[i],key[i],&MPI_COMM[i]);// Split the MPI communicator
    }                                                           // End loop over three types
#ifdef DEBUG
    print("level : ",0);                                        // Print identifier
    print(l,0);                                                 // Print current level
    print("\n",0);                                              // New line
    print("key   : \n",0);                                      // Print identifier
    print(key,0,3);                                             // Print key
    print("color : \n",0);                                      // Print identifier
    print(color,0,3);                                           // Print color
#endif
  }

  void bisectionAlltoall(Bodies &bodies, int nthLocal, int numLocal, int &newSize) {// One-to-one MPI_Alltoallv
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int scnt[2] = {nthLocal, numLocal - nthLocal};              // Set send count to right and left size
    int rcnt[2] = {0, 0};                                       // Initialize recv count
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[2]);    // Communicate the send count to get recv count
    int sdsp[2] = {0, scnt[0]};                                 // Send displacement
    int rdsp[2] = {0, rcnt[0]};                                 // Recv displacement
    newSize = rcnt[0] + rcnt[1];                                // Get new size from recv count
    if( color[0] != color[1] ) newSize = numLocal;              // Unless it's a leftover process of an odd group
    buffer.resize(newSize);                                     // Resize recv buffer

    for(int i=0; i!=2; ++i ) {                                  // Loop over 0 and 1
      scnt[i] *= bytes;                                         // Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         // Multiply send displacement by byte size of data
      rcnt[i] *= bytes;                                         // Multiply recv count by byte size of data
      rdsp[i] *= bytes;                                         // Multiply recv displacement by byte size of data
    }                                                           // End loop over 0 and 1
    MPI_Alltoallv(&bodies[0].I,scnt,sdsp,MPI_BYTE,              // Communicate bodies
                  &buffer[0].I,rcnt,rdsp,MPI_BYTE,              // Using same buffer as sort buffer
                  MPI_COMM[2]);                                 // MPI_COMM[2] is for the one-to-one pair
    if( color[0] == color[1] ) bodies = buffer;                 // Don't update if leftover process
    buffer.resize(bodies.size());                               // Resize sort buffer
    sort(bodies,buffer);                                        // Sort bodies in ascending order
  }

  void bisectionScatter(Bodies &bodies, int nthLocal, int &newSize) {// Scattering from leftover proc
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numScatter = nprocs[1] - 1;                             // Number of processes to scatter to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int *scnt = new int [nprocs[1]];                            // Send count
    int *sdsp = new int [nprocs[1]];                            // Send displacement
    int rcnt;                                                   // Recv count
    if( key[1] == numScatter ) {                                // If this is the leftover proc to scatter from
      sdsp[0] = 0;                                              //  Initialize send displacement
      for(int i=0; i!=numScatter; ++i ) {                       //  Loop over processes to send to
        int begin = 0, end = nthLocal;                          //  Set begin and end of range to send
        splitRange(begin,end,i,numScatter);                     //  Split range into numScatter and get i-th range
        scnt[i] = end - begin;                                  //  Set send count based on range
        sdsp[i+1] = sdsp[i] + scnt[i];                          //  Set send displacement based on send count
      }                                                         // End loop over processes to send to
      scnt[numScatter] = 0;                                     // Send count to self should be 0
      oldSize = 0;                                              // Reset oldSize to account for sent data
      newSize -= sdsp[numScatter];                              // Set newSize to account for sent data
      buffer.erase(buffer.begin(),buffer.begin()+sdsp[numScatter]);// Erase from recv buffer the part that was sent
    }                                                           // Endif for leftover proc
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[1]);// Scatter the send count to get recv count
    if( key[1] != numScatter ) {                                // If it's one the receiving procs
      newSize += rcnt;                                          //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize the recv buffer
    }                                                           // Endif for receiving procs
    for(int i=0; i!= nprocs[1]; ++i ) {                         // Loop over group of processes
      scnt[i] *= bytes;                                         //  Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         //  Multiply send displacement by byte size of data
    }                                                           // End loop over group of processes
    rcnt *= bytes;                                              // Multiply recv count by byte size of data
    MPI_Scatterv(&bodies[0].I,      scnt,sdsp,MPI_BYTE,         // Communicate bodies via MPI_Scatterv
                 &buffer[oldSize].I,rcnt,     MPI_BYTE,         // Offset recv buffer by oldSize
                 numScatter,MPI_COMM[1]);                       // MPI_COMM[1] is used for scatter
    bodies = buffer;                                            // Copy recv buffer to bodies
    if( key[1] != numScatter ) sort(bodies,buffer);             // Sort bodies in ascending order
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
  }

  void bisectionGather(Bodies &bodies, int nthLocal, int numLocal, int &newSize) {// Gathering to leftover proc
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numGather = nprocs[0] - 1;                              // Number of processes to gather to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int scnt;                                                   // Send count
    int *rcnt = new int [nprocs[0]];                            // Recv count
    int *rdsp = new int [nprocs[0]];                            // Recv displacement
    if( key[0] != 0 ) {                                         // If this is not the leftover proc
      int begin = 0, end = numLocal - nthLocal;                 //  Set begin and end of range to send
      splitRange(begin,end,key[0]-1,nprocs[0]);                 //  Split range into nprocs[0]
      scnt = end - begin;                                       //  Set send count based on range
      newSize -= scnt;                                          //  Set newSize to account for sent data
      buffer.erase(buffer.begin(),buffer.begin()+scnt);         //  Erase from recv buffer the part that was sent
    }                                                           // Endif for leftover proc
    MPI_Barrier(MPI_COMM[0]);                                   // Sync processes
    MPI_Gather(&scnt,1,MPI_INT,rcnt,1,MPI_INT,0,MPI_COMM[0]);   // Gather the send count to get recv count
    if( key[0] == 0 ) {                                         // If this is the leftover proc
      rdsp[0] = 0;                                              //  Initialize the recv displacement
      for(int i=0; i!=numGather; ++i ) {                        //  Loop over processes to gather from
        rdsp[i+1] = rdsp[i] + rcnt[i];                          //   Set recv displacement based on recv count
      }                                                         //  End loop over processes
      newSize += rdsp[numGather] + rcnt[numGather];             //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize recv buffer
    }                                                           // Endif for leftover proc
    scnt *= bytes;                                              // Multiply send count by byte size of data
    for(int i=0; i!= nprocs[0]; ++i ) {                         // Loop over group of processes
      rcnt[i] *= bytes;                                         //  Multiply recv count by byte size of data
      rdsp[i] *= bytes;                                         //  Multiply recv displacement by byte size of data
    }                                                           // End loop over group of processes
    MPI_Gatherv(&bodies[0].I,      scnt,     MPI_BYTE,          // Communicate bodies via MPI_Gatherv
                &buffer[oldSize].I,rcnt,rdsp,MPI_BYTE,          // Offset recv buffer by oldSize
                0,MPI_COMM[0]);                                 // MPI_COMM[0] is used for gather
    bodies = buffer;                                            // Copy recv buffer to bodies
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete send count
    if( key[0] == 0 ) sort(bodies,buffer);                      // Sort bodies in ascending order
  }

public:
  Partition() : TreeConstructor() {                             // Constructor
    LEVEL = int(log(SIZE) / M_LN2 - 1e-5) + 1;                  // Level of the process binary tree
    XMIN.resize(LEVEL+1);                                       // Minimum position vector at each level
    XMAX.resize(LEVEL+1);                                       // Maximum position vector at each level
  }
  ~Partition() {}                                               // Destructor

  void setGlobDomain(Bodies &bodies) {                          // Set bounds of domain to be partitioned
    numCells1D = 1 << getMaxLevel(bodies);                      // Set initial number of bodies
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN[0] = XMAX[0] = B->X;                                   // Initialize xmin,xmax
    MPI_Datatype MPI_TYPE = getType(XMIN[0][0]);                // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < XMIN[0][d]) XMIN[0][d] = B->X[d];     //   Determine xmin
        else if(B->X[d] > XMAX[0][d]) XMAX[0][d] = B->X[d];     //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over bodies
    vect X;                                                     // Recv buffer
    MPI_Allreduce(XMAX[0],X,3,MPI_TYPE,MPI_MAX,MPI_COMM_WORLD); // Reduce global maximum
    XMAX[0] = X;                                                // Get data from buffer
    MPI_Allreduce(XMIN[0],X,3,MPI_TYPE,MPI_MIN,MPI_COMM_WORLD); // Reduce global minimum
    XMIN[0] = X;                                                // Get data from buffer
    R0 = 0;                                                     // Initialize root radius
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = (XMAX[0][d] + XMIN[0][d]) / 2;                    //  Calculate center of domain
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(XMAX[0][d] - X0[d], R0);                    //  Calculate max distance from center
      R0 = std::max(X0[d] - XMIN[0][d], R0);                    //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1. + log(R0) / M_LN2));                     // Add some leeway to root radius
    XMAX[0] = X0 + R0;                                          // Reposition global maximum
    XMIN[0] = X0 - R0;                                          // Reposition global minimum
  }

  void binBodies(Bodies &bodies, int d) {                       // Turn positions into indices of bins
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->I = bigint((B->X[d] - XMIN[0][d])                      // Bin body positions into integers
        / (XMAX[0][d] - XMIN[0][d]) * numCells1D);
    }                                                           // End loop over bodies
  }

  int splitBodies(Bodies &bodies, bigint iSplit) {              // Split bodies according to iSplit
    int nth = 0;                                                // Initialize splitting index
    while( bodies[nth].I <= iSplit && nth < int(bodies.size()) ) nth++;// Determine end index of left side
    return nth;                                                 // Return end index of left side
  }

  void shiftBodies(Bodies &bodies) {                            // Send bodies to next rank (wrap around)
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    const int isend = (RANK + 1       ) % SIZE;                 // Send to next rank (wrap around)
    const int irecv = (RANK - 1 + SIZE) % SIZE;                 // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, recv request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0].I,oldSize*bytes,MPI_BYTE,irecv,        // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0].I,newSize*bytes,MPI_BYTE,isend,        // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

  template<typename T, typename T2>
  T2 nth_element(T &data, T2 n, MPI_Comm MPI_COMM0=0) {         // Find nth element of global data
    if( MPI_COMM0 == 0 ) {                                      // If MPI_COMM is not specified
      MPI_Comm_split(MPI_COMM_WORLD,0,RANK,&MPI_COMM0);         //  Create an artificial MPI_COMM
    }
    MPI_Comm_size(MPI_COMM0,&SIZES);                            // Get number of MPI processes for split comm
    MPI_Comm_rank(MPI_COMM0,&RANKS);                            // Get index of current MPI process for split comm
    int numData = data.size();                                  // Total size of data to perform nth_element
    int maxBucket = 1000;                                       // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int lOffset = 0;                                            // Local offset of region being considered
    MPI_Datatype MPI_TYPE = getType(n);                         // Get MPI data type
    int *rcnt = new int [SIZES];                                // MPI recv count
    Bigints send(maxBucket);                                    // MPI send buffer for data
    Bigints recv(maxBucket);                                    // MPI recv buffer for data
    T2 gOffset = 0;                                             // Global offset of region being considered
    T2 *isend = new T2 [maxBucket];                             // MPI send buffer for index
    T2 *irecv = new T2 [maxBucket];                             // MPI recv buffer for index
    T2 *iredu = new T2 [maxBucket];                             // Local scan of index buffer
    numBucket = getBucket(data,numData,lOffset,send,recv,MPI_COMM0);// Get global bucket data

    while( numBucket > 1 ) {                                    // While there are multipole candidates
      int ic=0, nth=0;                                          //  Initialize counters
      for( int i=0; i!=maxBucket; ++i ) isend[i] = 0;           //  Initialize bucket counter
      for( int i=0; i!=numData; ++i ) {                         //  Loop over range of data
        while( data[lOffset + i].I > recv[ic] && ic < numBucket-1 ) ++ic;// Set counter to current bucket
        isend[ic]++;                                            //   Increment bucket counter
      }                                                         //  End loop over data
      MPI_Reduce(isend,irecv,numBucket,MPI_TYPE,                //  Reduce bucket counter
                 MPI_SUM,0,MPI_COMM0);
      if( RANKS == 0 ) {                                        //  Only rank 0 operates on reduced data
        iredu[0] = 0;                                           //   Initialize global scan index
        for( int i=0; i!=numBucket-1; ++i ) {                   //   Loop over buckets
          iredu[i+1] = iredu[i] + irecv[i];                     //    Increment global scan index
        }                                                       //   End loop over buckets
        nth = 0;                                                //   Initialize index for bucket containing nth element
        while( n - gOffset > iredu[nth] && nth < numBucket ) ++nth;// Find index for bucket containing nth element
        --nth;                                                  //   Account for overshoot (do while?)
        if( nth == -1 ) nth = 0;                                //   If nth is -1 don't split
        gOffset += iredu[nth];                                  //   Increment global offset of region being considered
      }                                                         //  Endif for rank 0
      MPI_Bcast(&nth,1,MPI_INT,0,MPI_COMM0);                    //  Broadcast index for bucket containing nth element
      MPI_Bcast(&gOffset,1,MPI_TYPE,0,MPI_COMM0);               //  Broadcast global offset
      iredu[0] = 0;                                             //  Initialize local scan index
      for( int i=0; i!=numBucket-1; ++i ) {                     //  Loop over buckets
        iredu[i+1] = iredu[i] + isend[i];                       //   Increment local scan index
      }                                                         //  End loop over buckets
      if( nth == numBucket-1 ) {                                //  If nth is last bucket
        numData = numData-iredu[nth];                           //   Region of interest is that bucket
      } else {                                                  //  If nth is not the last bucket
        numData = iredu[nth+1]-iredu[nth];                      //   Region of interest is that bucket
      }                                                         //  Endif for last bucket
      lOffset += iredu[nth];                                    //  Increment local offset to that bucket
      numBucket = getBucket(data,numData,lOffset,send,recv,MPI_COMM0);// Get global bucket data
    }                                                           // End while loop
    delete[] rcnt;                                              // Delete recv count
    delete[] isend;                                             // Delete send buffer for index
    delete[] irecv;                                             // Delete recv buffer for index
    delete[] iredu;                                             // Delete local scan of index buffer
    return recv[0]-1;                                           // Return nth element
  }

  void bisection(Bodies &bodies) {                              // Partitioning by recursive bisection
    startTimer("Partition    ");                                // Start timer
    MPI_Datatype MPI_TYPE = getType(bodies[0].I);               // Get MPI data type
    nprocs[0] = nprocs[1] = SIZE;                               // Initialize number of processes in groups
    offset[0] = offset[1] = 0;                                  // Initialize offset of body in groups
     color[0] =  color[1] =  color[2] = 0;                      // Initialize color of communicators
       key[0] =    key[1] =    key[2] = 0;                      // Initialize key of communicators
    int newSize;                                                // New size of recv buffer
    bigint numLocal = bodies.size();                            // Local data size
    bigint numGlobal;                                           // Global data size
    MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM_WORLD);// Reduce Global data size
    bigint nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];// Split at nth global element
    binBodies(bodies,2);                                        // Bin bodies into leaf level cells
    buffer.resize(numLocal);                                    // Resize sort buffer
    sort(bodies,buffer);                                        // Sort bodies in ascending order
    bigint iSplit = nth_element(bodies,nthGlobal);              // Get cell index of nth global element
    int nthLocal = splitBodies(bodies,iSplit);                  // Split bodies based on iSplit
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      bisectionGetComm(l);                                      //  Split the MPI communicator for that level
      splitDomain(iSplit,l,2-l%3);                              //  Split the domain according to iSplit
      bisectionAlltoall(bodies,nthLocal,numLocal,newSize);      //  Communicate bodies by one-to-one MPI_Alltoallv
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] <= nprocs[1] )// If scatter is necessary
        bisectionScatter(bodies,nthLocal,newSize);              //  Communicate bodies by scattering from leftover proc
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] >= nprocs[1] )// If gather is necessary
        bisectionGather(bodies,nthLocal,numLocal,newSize);      //  Communicate bodies by gathering to leftover proc
#ifdef DEBUG
      for(int j=0; j!=SIZE; ++j) {                              //  Loop over ranks
        MPI_Barrier(MPI_COMM_WORLD);                            //   Sync processes
        usleep(WAIT);                                           //   Wait for "WAIT" milliseconds
        if( RANK == j ) {                                       //   If it's my turn to print
          std::cout << "RANK : " << j << std::endl;             //    Print rank
          for(int i=0; i!=9; ++i) std::cout << bodies[numLocal/10*i].I << " ";// Print sampled body indices
          std::cout << std::endl;                               //    New line
        }                                                       //   Endif for my turn
      }                                                         //  End loop over ranks
#endif
      numLocal = newSize;                                       //  Update local data size
      MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM[0]);// Reduce global data size
      nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];    //  Split at nth global element
      binBodies(bodies,2-(l+1)%3);                              //  Bin bodies into leaf level cells
      buffer.resize(numLocal);                                  //  Resize sort buffer
      sort(bodies,buffer);                                      //  Sort bodies in ascending order
      iSplit = nth_element(bodies,nthGlobal,MPI_COMM[0]);       //  Get cell index of nth global element
      nthLocal = splitBodies(bodies,iSplit);                    //  Split bodies based on iSplit
    }                                                           // End loop over levels of N-D hypercube communication
    stopTimer("Partition    ",printNow);                        // Stop timer & print
  }

  void octsection(Bodies &bodies) {                             // Partition by recursive octsection
    startTimer("Partition    ");                                // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int level = int(log(SIZE + 1) / M_LN2 / 3);                 // Max level/3 of N-D hypercube communication
    BottomUp::setIndex(bodies,level);                           // Set index of bodies for that level
    buffer.resize(bodies.size());                               // Resize sort buffer
    sort(bodies,buffer);                                        // Sort bodies in ascending order
    int *scnt = new int [SIZE];                                 // Send count
    int *sdsp = new int [SIZE];                                 // Send displacement
    int *rcnt = new int [SIZE];                                 // Recv count
    int *rdsp = new int [SIZE];                                 // Recv displacement
    for( int i=0; i!=SIZE; ++i ) {                              // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int i = B->I - ((1 << 3*level) - 1) / 7;                  //  Get levelwise index
      scnt[i]++;                                                //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=SIZE-1; ++irank ) {                // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[SIZE-1]+rcnt[SIZE-1]);                   // Resize recv buffer
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    for( int l=0; l!=level; ++l ) {                             // Loop over level/3 of N-D hypercube
      for( int d=0; d!=3; ++d ) {                               //  Loop over 3 dimensions
        int i = 3 * l + d;                                      //   i = actual level of N-D hypercube
        XMIN[i+1] = XMIN[i];                                    //   Set XMAX for next subdivision
        XMAX[i+1] = XMAX[i];                                    //   Set XMIN for next subdivision
        if( (RANK >> (3 * (level - l - 1) + 2 - d)) % 2 ) {     //   If on left side
          XMIN[i+1][2-d] = (XMAX[i][2-d]+XMIN[i][2-d]) / 2;     //    Set XMIN to midpoint
        } else {                                                //   If on right side
          XMAX[i+1][2-d] = (XMAX[i][2-d]+XMIN[i][2-d]) / 2;     //    Set XMAX to midpoint
        }                                                       //   Endif for side
      }                                                         //  End loop over dimensions
    }                                                           // End loop over levels
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Partition    ",printNow);                        // Stop timer & print
  }
};
}

#endif
