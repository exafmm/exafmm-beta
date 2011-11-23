/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "serialfmm.h"

//! Handles all the partitioning of domains
template<Equation equation>
class Partition : public MyMPI, public SerialFMM<equation> {
private:
  int numCells1D;                                               //!< Number of cells in one dimension (leaf level)

protected:
  int LEVEL;                                                    //!< Level of the MPI process binary tree
  std::vector<vect> XMIN;                                       //!< Minimum position vector of bodies
  std::vector<vect> XMAX;                                       //!< Maximum position vector of bodies
  int nprocs[64][2];                                            //!< Number of processes in the two split groups
  int offset[64][2];                                            //!< Offset of body in the two split groups
  int  color[64][3];                                            //!< Color of hypercube communicators
  int    key[64][3];                                            //!< Key of hypercube communicators
  MPI_Comm MPI_COMM[64][3];                                     //!< Hypercube communicators

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using TreeStructure<equation>::buffer;                        //!< Buffer for MPI communication & sorting
  using TreeStructure<equation>::getLevel;                      //!< Get level from cell index
  using BottomUp<equation>::getMaxLevel;                        //!< Max level for bottom up tree build

private:
//! Split domain according to iSplit
  void splitDomain(bigint iSplit, int l, int d) {
    real X = (iSplit + 1) * (XMAX[l][d] - XMIN[l][d]) / numCells1D + XMIN[l][d];// Coordinate corresponding to iSplit
    XMAX[l+1] = XMAX[l];                                        // Set XMAX for next subdivision
    XMIN[l+1] = XMIN[l];                                        // Set XMIN for next subdivision
    if( color[l+1][0] % 2 == 0 ) {                              // If on left side
      XMAX[l+1][d] = X;                                         //  Set max to X
    } else {                                                    // If on right side
      XMIN[l+1][d] = X;                                         //  Set min to X
    }                                                           // Endif for sides
  }

//! Get global bucket data for parallel nth_element
  template<typename T>
  int getBucket(T &data, int numData, int lOffset, Bigints &send, Bigints &recv, MPI_Comm MPI_COMM0) {
    int maxBucket = send.size();                                // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample = std::min(maxBucket/MPISIZES,numData);       // Number of local samples
    MPI_Datatype MPI_TYPE = getType(data[0].ICELL);             // Get MPI data type
    int *rcnt = new int [MPISIZES];                             // MPI recv count
    int *rdsp = new int [MPISIZES];                             // MPI recv displacement
    for( int i=0; i!=numSample; ++i ) {                         // Loop over local samples
      int stride = numData/numSample;                           //  Sampling stride
      send[i] = data[lOffset + i * stride].ICELL;               //  Put sampled data in send buffer
    }                                                           // End loop over samples
    MPI_Gather(&numSample,1,MPI_INT,                            // Gather size of sample data to rank 0
               rcnt,      1,MPI_INT,
               0,MPI_COMM0);
    if( MPIRANKS == 0 ) {                                       // Only rank 0 operates on gathered info
      numBucket = 0;                                            //  Initialize number of buckets
      for( int irank=0; irank!=MPISIZES; ++irank ) {            //  Loop over processes
        rdsp[irank] = numBucket;                                //   Calculate recv displacement
        numBucket += rcnt[irank];                               //   Accumulate recv count to get number of buckets
      }                                                         //  End loop over processes
      recv.resize(numBucket);                                   //  Resize recv so that end() is valid
    }                                                           // Endif for rank 0
    MPI_Gatherv(&send[0],numSample,MPI_TYPE,                    // Gather sample data to rank 0
                &recv[0],rcnt,rdsp,MPI_TYPE,
                0,MPI_COMM0);
    if( MPIRANKS == 0 ) {                                       // Only rank 0 operates on gathered info
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
//! Split the MPI communicator into N-D hypercube
  void bisectionGetComm(int l) {
    for( int i=1; i>=0; --i ) {                                 // Loop over 1 and 0
      int pSplit = (nprocs[l][0] + i) / 2;                      //  Splitting criteria
      if( MPIRANK - offset[l][0] < pSplit ) {                   //  If on left side
        nprocs[l+1][i] = pSplit;                                //   Group size is the splitting criteria
        offset[l+1][i] = offset[l][0];                          //   Offset is the same as previous
         color[l+1][i] =  color[l][0] * 2;                      //   Color is twice the previous
      } else {                                                  //  If on right side
        nprocs[l+1][i] = nprocs[l][0] - pSplit;                 //   Group size is what remains from the left side
        offset[l+1][i] = offset[l][0] + pSplit;                 //   Offset is incremented by splitting criteria
         color[l+1][i] =  color[l][0] * 2 + 1;                  //   Color is twice the previous plus one
      }                                                         //  Endif for sides
      key[l+1][i] = MPIRANK - offset[l+1][i];                   //  Key is determined from offset
    }                                                           // End loop over 1 and 0
    key[l+1][2] = color[l+1][1] % 2;                            // The third type of key is determined from color[1]
    color[l+1][2] = key[l+1][1] + color[l][0] * (1 << (LEVEL - l - 1));// The corresponding color is determined from key[1]
    for( int i=0; i!=3; ++i ) {                                 // Loop over the three types of colors and keys
      MPI_Comm_split(MPI_COMM_WORLD,color[l+1][i],key[l+1][i],&MPI_COMM[l+1][i]);// Split the MPI communicator
    }                                                           // End loop over three types
#ifdef DEBUG
    print("level : ",0);                                        // Print identifier
    print(l+1,0);                                               // Print current level
    print("\n",0);                                              // New line
    print("key   : \n",0);                                      // Print identifier
    print(key[l+1],0,3);                                        // Print key
    print("color : \n",0);                                      // Print identifier
    print(color[l+1],0,3);                                      // Print color
#endif
  }

//! One-to-one MPI_Alltoallv
  void bisectionAlltoall(Bodies &bodies, int nthLocal, int numLocal, int &newSize, int l) {
    startTimer("Bi Alltoall  ");                                // Start timer
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int scnt[2] = {nthLocal, numLocal - nthLocal};              // Set send count to right and left size
    int rcnt[2] = {0, 0};                                       // Initialize recv count
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[l+1][2]);// Communicate the send count to get recv count
    int sdsp[2] = {0, scnt[0]};                                 // Send displacement
    int rdsp[2] = {0, rcnt[0]};                                 // Recv displacement
    newSize = rcnt[0] + rcnt[1];                                // Get new size from recv count
    if( color[l+1][0] != color[l+1][1] ) newSize = numLocal;    // Unless it's a leftover process of an odd group
    buffer.resize(newSize);                                     // Resize recv buffer

    for(int i=0; i!=2; ++i ) {                                  // Loop over 0 and 1
      scnt[i] *= bytes;                                         // Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         // Multiply send displacement by byte size of data
      rcnt[i] *= bytes;                                         // Multiply recv count by byte size of data
      rdsp[i] *= bytes;                                         // Multiply recv displacement by byte size of data
    }                                                           // End loop over 0 and 1
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,                // Communicate bodies
                  &buffer[0],rcnt,rdsp,MPI_BYTE,                // Using same buffer as sort buffer
                  MPI_COMM[l+1][2]);                            // MPI_COMM[2] is for the one-to-one pair
    if( color[l+1][0] == color[l+1][1] ) bodies = buffer;       // Don't update if leftover process
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Bi Alltoall  ",printNow);                        // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
  }

//! Scattering from leftover processes
  void bisectionScatter(Bodies &bodies, int nthLocal, int &newSize, int l) {
    startTimer("Bi Scatter   ");                                // Start timer
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numScatter = nprocs[l+1][1] - 1;                        // Number of processes to scatter to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int *scnt = new int [nprocs[l+1][1]];                       // Send count
    int *sdsp = new int [nprocs[l+1][1]];                       // Send displacement
    int rcnt;                                                   // Recv count
    if( key[l+1][1] == numScatter ) {                           // If this is the leftover proc to scatter from
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
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[l+1][1]);// Scatter the send count to get recv count
    if( key[l+1][1] != numScatter ) {                           // If it's one the receiving procs
      newSize += rcnt;                                          //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize the recv buffer
    }                                                           // Endif for receiving procs
    for(int i=0; i!= nprocs[l+1][1]; ++i ) {                    // Loop over group of processes
      scnt[i] *= bytes;                                         //  Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         //  Multiply send displacement by byte size of data
    }                                                           // End loop over group of processes
    rcnt *= bytes;                                              // Multiply recv count by byte size of data
    MPI_Scatterv(&bodies[0],      scnt,sdsp,MPI_BYTE,           // Communicate bodies via MPI_Scatterv
                 &buffer[oldSize],rcnt,     MPI_BYTE,           // Offset recv buffer by oldSize
                 numScatter,MPI_COMM[l+1][1]);                  // MPI_COMM[1] is used for scatter
    bodies = buffer;                                            // Copy recv buffer to bodies
    if( key[l+1][1] != numScatter ) sortBodies(bodies,buffer);  // Sort bodies in ascending order
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    stopTimer("Bi Scatter   ",printNow);                        // Stop timer 
  }

//! Gathering to leftover processes
  void bisectionGather(Bodies &bodies, int nthLocal, int numLocal, int &newSize, int l) {
    startTimer("Bi Gather    ");                                // Start timer
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int numGather = nprocs[l+1][0] - 1;                         // Number of processes to gather to
    int oldSize = newSize;                                      // Size of recv buffer before communication
    int scnt;                                                   // Send count
    int *rcnt = new int [nprocs[l+1][0]];                       // Recv count
    int *rdsp = new int [nprocs[l+1][0]];                       // Recv displacement
    if( key[l+1][0] != 0 ) {                                    // If this is not the leftover proc
      int begin = 0, end = numLocal - nthLocal;                 //  Set begin and end of range to send
      splitRange(begin,end,key[l+1][0]-1,nprocs[l+1][0]);       //  Split range into nprocs[0]
      scnt = end - begin;                                       //  Set send count based on range
      newSize -= scnt;                                          //  Set newSize to account for sent data
      buffer.erase(buffer.begin(),buffer.begin()+scnt);         //  Erase from recv buffer the part that was sent
    }                                                           // Endif for leftover proc
    MPI_Barrier(MPI_COMM[l+1][0]);                              // Sync processes
    MPI_Gather(&scnt,1,MPI_INT,rcnt,1,MPI_INT,0,MPI_COMM[l+1][0]);// Gather the send count to get recv count
    if( key[l+1][0] == 0 ) {                                    // If this is the leftover proc
      rdsp[0] = 0;                                              //  Initialize the recv displacement
      for(int i=0; i!=numGather; ++i ) {                        //  Loop over processes to gather from
        rdsp[i+1] = rdsp[i] + rcnt[i];                          //   Set recv displacement based on recv count
      }                                                         //  End loop over processes
      newSize += rdsp[numGather] + rcnt[numGather];             //  Increment newSize by the recv count
      buffer.resize(newSize);                                   //  Resize recv buffer
    }                                                           // Endif for leftover proc
    scnt *= bytes;                                              // Multiply send count by byte size of data
    for(int i=0; i!= nprocs[l+1][0]; ++i ) {                    // Loop over group of processes
      rcnt[i] *= bytes;                                         //  Multiply recv count by byte size of data
      rdsp[i] *= bytes;                                         //  Multiply recv displacement by byte size of data
    }                                                           // End loop over group of processes
    MPI_Gatherv(&bodies[0],      scnt,     MPI_BYTE,            // Communicate bodies via MPI_Gatherv
                &buffer[oldSize],rcnt,rdsp,MPI_BYTE,            // Offset recv buffer by oldSize
                0,MPI_COMM[l+1][0]);                            // MPI_COMM[0] is used for gather
    bodies = buffer;                                            // Copy recv buffer to bodies
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete send count
    if( key[l+1][0] == 0 ) sortBodies(bodies,buffer);           // Sort bodies in ascending order
    stopTimer("Bi Gather    ",printNow);                        // Stop timer 
  }

public:
//! Constructor
  Partition() : SerialFMM<equation>() {
    LEVEL = int(log(MPISIZE) / M_LN2 - 1e-5) + 1;               // Level of the process binary tree
    if(MPISIZE == 1) LEVEL = 0;                                 // Level is 0 for a serial execution
    XMIN.resize(LEVEL+1);                                       // Minimum position vector at each level
    XMAX.resize(LEVEL+1);                                       // Maximum position vector at each level
    startTimer("Split comm   ");                                // Start timer
    nprocs[0][0] = nprocs[0][1] = MPISIZE;                      // Initialize number of processes in groups
    offset[0][0] = offset[0][1] = 0;                            // Initialize offset of body in groups
     color[0][0] =  color[0][1] =  color[0][2] = 0;             // Initialize color of communicators
       key[0][0] =    key[0][1] =    key[0][2] = 0;             // Initialize key of communicators
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      bisectionGetComm(l);                                      //  Split the MPI communicator for that level
    }                                                           // End loop over levels of N-D hypercube communication
    stopTimer("Split comm   ",printNow);                        // Stop timer 
  }
//! Destructor
  ~Partition() {}

//! Set bounds of domain to be partitioned
  void setGlobDomain(Bodies &bodies, vect x0=0, real r0=M_PI) {
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
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = (XMAX[0][d] + XMIN[0][d]) / 2;                    //  Calculate center of domain
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(XMAX[0][d] - X0[d], R0);                    //  Calculate max distance from center
      R0 = std::max(X0[d] - XMIN[0][d], R0);                    //  Calculate max distance from center
    }                                                           // End loop over each dimension
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( X0[0]-R0 < x0[0]-r0 || x0[0]+r0 < X0[0]+R0            //  Check for outliers in x direction
       || X0[1]-R0 < x0[1]-r0 || x0[1]+r0 < X0[1]+R0            //  Check for outliers in y direction
       || X0[2]-R0 < x0[2]-r0 || x0[2]+r0 < X0[2]+R0 ) {        //  Check for outliers in z direction
        std::cout << "Error: Particles located outside periodic domain @ rank "
                  << MPIRANK << std::endl;                      //   Print error message
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is M_PI
    } else {                                                    // If not periodic boundary condition
      R0 += 1e-5;                                               //  Add some leeway to root radius
    }                                                           // Endif for periodic boundary condition
    XMAX[0] = X0 + R0;                                          // Reposition global maximum
    XMIN[0] = X0 - R0;                                          // Reposition global minimum
  }

//! Turn positions into indices of bins
  void binBodies(Bodies &bodies, int d) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = bigint((B->X[d] - XMIN[0][d])                  // Bin body positions into integers
        / (XMAX[0][d] - XMIN[0][d]) * numCells1D);
    }                                                           // End loop over bodies
  }

//! Split bodies according to iSplit
  int splitBodies(Bodies &bodies, bigint iSplit) {
    int nth = 0;                                                // Initialize splitting index
    while( bodies[nth].ICELL <= iSplit && nth < int(bodies.size()) ) nth++;// Determine end index of left side
    return nth;                                                 // Return end index of left side
  }

//! Send bodies to next rank (round robin)
  void shiftBodies(Bodies &bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    const int isend = (MPIRANK + 1          ) % MPISIZE;        // Send to next rank (wrap around)
    const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, recv request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0],oldSize*bytes,MPI_BYTE,irecv,          // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0],newSize*bytes,MPI_BYTE,isend,          // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

//! Parallel global nth_element on distributed memory
  template<typename T, typename T2>
  T2 nth_element(T &data, T2 n, MPI_Comm MPI_COMM0=0) {
    if( MPI_COMM0 == 0 ) {                                      // If MPI_COMM is not specified
      MPI_Comm_split(MPI_COMM_WORLD,0,MPIRANK,&MPI_COMM0);      //  Create an artificial MPI_COMM
    }
    MPI_Comm_size(MPI_COMM0,&MPISIZES);                         // Get number of MPI processes for split comm
    MPI_Comm_rank(MPI_COMM0,&MPIRANKS);                         // Get index of current MPI process for split comm
    int numData = data.size();                                  // Total size of data to perform nth_element
    int maxBucket = 1000;                                       // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int lOffset = 0;                                            // Local offset of region being considered
    MPI_Datatype MPI_TYPE = getType(n);                         // Get MPI data type
    int *rcnt = new int [MPISIZES];                             // MPI recv count
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
        while( data[lOffset + i].ICELL > recv[ic] && ic < numBucket-1 ) ++ic;// Set counter to current bucket
        isend[ic]++;                                            //   Increment bucket counter
      }                                                         //  End loop over data
      MPI_Reduce(isend,irecv,numBucket,MPI_TYPE,                //  Reduce bucket counter
                 MPI_SUM,0,MPI_COMM0);
      if( MPIRANKS == 0 ) {                                     //  Only rank 0 operates on reduced data
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

//! Partitioning by recursive bisection
  void bisection(Bodies &bodies) {
    startTimer("Bin bodies   ");                                // Start timer
    MPI_Datatype MPI_TYPE = getType(bodies[0].ICELL);           // Get MPI data type
    int newSize;                                                // New size of recv buffer
    bigint numLocal = bodies.size();                            // Local data size
    bigint numGlobal;                                           // Global data size
    MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM_WORLD);// Reduce Global data size
    bigint nthGlobal = (numGlobal * (nprocs[0][0] / 2)) / nprocs[0][0];// Split at nth global element
    binBodies(bodies,2);                                        // Bin bodies into leaf level cells
    buffer.resize(numLocal);                                    // Resize sort buffer
    stopTimer("Bin bodies   ",printNow);                        // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Split bodies ");                                // Start timer
    bigint iSplit = nth_element(bodies,nthGlobal);              // Get cell index of nth global element
    int nthLocal = splitBodies(bodies,iSplit);                  // Split bodies based on iSplit
    stopTimer("Split bodies ",printNow);                        // Stop timer 
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      splitDomain(iSplit,l,2-l%3);                              //  Split the domain according to iSplit
      bisectionAlltoall(bodies,nthLocal,numLocal,newSize,l);    //  Communicate bodies by one-to-one MPI_Alltoallv
      if( nprocs[l][0] % 2 == 1 && nprocs[l][0] != 1 && nprocs[l+1][0] <= nprocs[l+1][1] )// If scatter is necessary
        bisectionScatter(bodies,nthLocal,newSize,l);            //  Communicate bodies by scattering from leftover proc
      if( nprocs[l][0] % 2 == 1 && nprocs[l][0] != 1 && nprocs[l+1][0] >= nprocs[l+1][1] )// If gather is necessary
        bisectionGather(bodies,nthLocal,numLocal,newSize,l);    //  Communicate bodies by gathering to leftover proc
#ifdef DEBUG
      for(int j=0; j!=MPISIZE; ++j) {                           //  Loop over ranks
        MPI_Barrier(MPI_COMM_WORLD);                            //   Sync processes
        usleep(WAIT);                                           //   Wait for "WAIT" milliseconds
        if( MPIRANK == j ) {                                    //   If it's my turn to print
          std::cout << "MPIRANK : " << j << std::endl;          //    Print rank
          for(int i=0; i!=9; ++i) std::cout << bodies[numLocal/10*i].I << " ";// Print sampled body indices
          std::cout << std::endl;                               //    New line
        }                                                       //   Endif for my turn
      }                                                         //  End loop over ranks
#endif
      startTimer("Bin bodies   ");                              //  Start timer
      numLocal = newSize;                                       //  Update local data size
      MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM[l+1][0]);// Reduce global data size
      nthGlobal = (numGlobal * (nprocs[l+1][0] / 2)) / nprocs[l+1][0];//  Split at nth global element
      binBodies(bodies,2-(l+1)%3);                              //  Bin bodies into leaf level cells
      buffer.resize(numLocal);                                  //  Resize sort buffer
      stopTimer("Bin bodies   ",printNow);                      //  Stop timer 
      sortBodies(bodies,buffer);                                //  Sort bodies in ascending order
      startTimer("Split bodies ");                              //  Start timer
      iSplit = nth_element(bodies,nthGlobal,MPI_COMM[l+1][0]);  //  Get cell index of nth global element
      nthLocal = splitBodies(bodies,iSplit);                    //  Split bodies based on iSplit
      stopTimer("Split bodies ",printNow);                      //  Stop timer 
    }                                                           // End loop over levels of N-D hypercube communication
  }

//! Partition by recursive octsection
  void octsection(Bodies &bodies) {
    startTimer("Partition    ");                                // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // For serial execution local root cell is root cell
    BottomUp<equation>::setIndex(bodies,level);                 // Set index of bodies for that level
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Partition    ");                                 // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Partition    ");                                // Start timer
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int index = B->ICELL - ((1 << 3*level) - 1) / 7;          //  Get levelwise index
      int irank = index / (int(pow(8,level)) / MPISIZE);        //  Get rank which the cell belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over level of N-D hypercube
      int d = 2 - l % 3;                                        //  Dimension of subdivision
      XMIN[l+1] = XMIN[l];                                      //  Set XMAX for next subdivision
      XMAX[l+1] = XMAX[l];                                      //  Set XMIN for next subdivision
      if( (MPIRANK >> (LEVEL - l - 1)) % 2 ) {                  //  If on left side
        XMIN[l+1][d] = (XMAX[l][d]+XMIN[l][d]) / 2;             //   Set XMIN to midpoint
      } else {                                                  //  If on right side
        XMAX[l+1][d] = (XMAX[l][d]+XMIN[l][d]) / 2;             //   Set XMAX to midpoint
      }                                                         //  Endif for side
    }                                                           // End loop over levels
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Partition    ",printNow);                        // Stop timer 
  }

//! Send bodies back to where they came from
  void unpartition(Bodies &bodies) {
    startTimer("Unpartition  ");                                // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = B->IPROC;                                      //  Copy process rank to cell index for sorting
    }                                                           // End loop over bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Unpartition  ");                                 // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Unpartition  ");                                // Start timer
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int irank = B->IPROC;                                     //  Get rank which the body belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      assert(B->IPROC == MPIRANK);                              //  Make sure bodies are in the right rank
    }                                                           // End loop over bodies
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Unpartition  ",printNow);                        // Stop timer 
  }
};

#endif
