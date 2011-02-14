#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "construct.h"

class Partition : public MyMPI, public TreeConstructor {
private:
  int numBodies;                                                // Initial local number of bodies
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
public:
  Partition() : TreeConstructor() {                             // Constructor
    LEVEL = int(log(SIZE) / M_LN2 - 1e-5) + 1;                  // Level of the process binary tree
    XMIN.resize(LEVEL+1);                                       // Minimum position vector at each level
    XMAX.resize(LEVEL+1);                                       // Maximum position vector at each level
  }
  ~Partition() {}                                               // Destructor

  void setGlobDomain(Bodies &bodies) {                          // Set bounds of domain to be partitioned
    numBodies = bodies.size();                                  // Set initial number of bodies
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN[0] = XMAX[0] = B->pos;                                 // Initialize xmin,xmax
    int const MPI_TYPE = getType(XMIN[0][0]);                   // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < XMIN[0][d]) XMIN[0][d] = B->pos[d]; //   Determine xmin
        else if(B->pos[d] > XMAX[0][d]) XMAX[0][d] = B->pos[d]; //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over all bodies
    vect X;                                                     // Define recv buffer
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
  }

  void splitDomain(bigint iSplit, int l, int d) {
    real X = iSplit * (XMAX[l][d] - XMIN[l][d]) / numBodies + XMIN[l][d];
    XMAX[l+1] = XMAX[l];
    XMIN[l+1] = XMIN[l];
    if( color[0] % 2 == 0 ) {
      XMAX[l+1][d] = X;
    } else {
      XMIN[l+1][d] = X;
    }
  }

  void binBodies(Bodies &bodies, int d) {                       // Turn positions into indices of bins
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      B->I = bigint((B->pos[d] - XMIN[0][d])                    // Bin body positions into integers
        / (XMAX[0][d] - XMIN[0][d]) * numBodies);
    }                                                           // End loop over all bodies
  }

  void shiftBodies(Bodies &bodies) {                            // Send bodies to next rank (wrap around)
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    int const bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int const isend = (RANK + 1       ) % SIZE;                 // Send to next rank (wrap around)
    int const irecv = (RANK - 1 + SIZE) % SIZE;                 // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0].I,oldSize*bytes,MPI_BYTE,irecv,        // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0].I,newSize*bytes,MPI_BYTE,isend,        // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

  int splitBodies(Bodies &bodies, bigint iSplit) {
    int nth = 0;
    while( bodies[nth].I < iSplit && nth < int(bodies.size() - 1) ) nth++;
    return nth;
  }

  template<typename T>
  int getBucket(T &data, int numData, int lOffset, Bigints &send, Bigints &recv, MPI_Comm MPI_COMM0) {
    int maxBucket = send.size();                                // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample = std::min(maxBucket/SIZES,numData);          // Number of local samples
    int const MPI_TYPE = getType(data[0].I);                    // Get MPI data type
    int *rcnt = new int [SIZES];                                // MPI receive count
    int *rdsp = new int [SIZES];                                // MPI receive displacement
    for( int i=0; i!=numSample; ++i ) {                         // Loop over local samples
      int stride = numData/numSample;                           //  Sampling stride
      send[i] = data[lOffset + i * stride].I;                   //  Put sampled data in send buffer
    }                                                           // End loop over samples
    MPI_Gather(&numSample,1,MPI_INT,                            // Gather size of sample data to rank 0
               rcnt,      1,MPI_INT,
               0,MPI_COMM0);
    if( RANKS == 0 ) {                                          // Only rank 0 operates on gathered info
      numBucket = 0;                                            //  Initialize number of buckets
      for( int irank=0; irank!=SIZES; ++irank ) {               //  Loop over all processes
        rdsp[irank] = numBucket;                                //   Calculate receive displacement
        numBucket += rcnt[irank];                               //   Accumulate receive count to get number of buckets
      }                                                         //  End loop over all processes
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
    delete[] rcnt;                                              // Delete receive count
    delete[] rdsp;                                              // Delete receive displacement
    return numBucket;                                           // Return number of buckets
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
    int const MPI_TYPE = getType(n);                            // Get MPI data type
    int *rcnt = new int [SIZES];                                // MPI receive count
    Bigints send(maxBucket);                                    // MPI send buffer for data
    Bigints recv(maxBucket);                                    // MPI receive buffer for data
    T2 gOffset = 0;                                             // Global offset of region being considered
    T2 *isend = new T2 [maxBucket];                             // MPI send buffer for index
    T2 *irecv = new T2 [maxBucket];                             // MPI receive buffer for index
    T2 *iredu = new T2 [maxBucket];                             // Local scan of index buffer
    numBucket = getBucket(data,numData,lOffset,send,recv,MPI_COMM0);// Get global bucket data
    MPI_Barrier(MPI_COMM_WORLD);

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
    delete[] rcnt;                                              // Delete receive count
    delete[] isend;                                             // Delete send buffer for index
    delete[] irecv;                                             // Delete receive buffer for index
    delete[] iredu;                                             // Delete local scan of index buffer
    return recv[0];                                             // Return nth element
  }

  void bisectionGetComm(int l) {
    int  oldcolor = color[0];
    oldnprocs = nprocs[0];
    for( int i=1; i>=0; --i ) {
      int pSplit = (oldnprocs + i) / 2;
      if( RANK - offset[0] < pSplit ) {
        nprocs[i] = pSplit;
        offset[i] = offset[0];
         color[i] =  color[0] * 2;
      } else {
        nprocs[i] = oldnprocs - pSplit;
        offset[i] = offset[0] + pSplit;
         color[i] =  color[0] * 2 + 1;
      }
      key[i] = RANK - offset[i];
    }
    key[2] = color[1] % 2;
    color[2] = key[1] + oldcolor * (1 << (LEVEL - l - 1));
    for( int i=0; i!=3; ++i ) {
      MPI_Comm_split(MPI_COMM_WORLD,color[i],key[i],&MPI_COMM[i]);
    }
#ifdef DEBUG
    print("level : ",0);
    print(l,0);
    print("\n",0);
    print("key   : \n",0);
    print(key,0,3);
    print("color : \n",0);
    print(color,0,3);
#endif
  }

  void bisectionAlltoall(Bodies &bodies, int nthLocal, int numLocal, int &newSize) {
    int const bytes = sizeof(bodies[0]);
    int scnt[2] = {nthLocal, numLocal - nthLocal};
    int rcnt[2] = {0, 0};
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[2]);
    int sdsp[2] = {0, scnt[0]};
    int rdsp[2] = {0, rcnt[0]};
    newSize = rcnt[0] + rcnt[1];
    if( color[0] != color[1] ) newSize = numLocal;
    buffer.resize(newSize);

    for(int i=0; i!=2; ++i ) {
      scnt[i] *= bytes;
      sdsp[i] *= bytes;
      rcnt[i] *= bytes;
      rdsp[i] *= bytes;
    }
    MPI_Alltoallv(&bodies[0].I,scnt,sdsp,MPI_BYTE,
                  &buffer[0].I,rcnt,rdsp,MPI_BYTE,
                  MPI_COMM[2]);
    if( color[0] == color[1] && nthLocal != 0 ) bodies = buffer;
    buffer.resize(bodies.size());
    sort(bodies,buffer);
  }

  void bisectionScatter(Bodies &bodies, int nthLocal, int &newSize) {
    int const bytes = sizeof(bodies[0]);
    int numScatter = nprocs[1] - 1;
    int oldSize = newSize;
    int *scnt = new int [nprocs[1]];
    int *sdsp = new int [nprocs[1]];
    int rcnt;
    if( key[1] == numScatter ) {
      sdsp[0] = 0;
      for(int i=0; i!=numScatter; ++i ) {
        scnt[i] = nthLocal / numScatter;
        sdsp[i+1] = sdsp[i] + scnt[i];
      }
      scnt[numScatter] = 0;
      oldSize = 0;
      newSize -= sdsp[numScatter];
      buffer.erase( buffer.begin(), buffer.begin()+sdsp[numScatter]);
    }
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[1]);
    if( key[1] != numScatter ) {
      newSize += rcnt;
      buffer.resize(newSize);
    }

    for(int i=0; i!= nprocs[1]; ++i ) {
      scnt[i] *= bytes;
      sdsp[i] *= bytes;
    }
    rcnt *= bytes;

    MPI_Scatterv(&bodies[0].I,      scnt,sdsp,MPI_BYTE,
                 &buffer[oldSize].I,rcnt,     MPI_BYTE,
                 numScatter,MPI_COMM[1]);
    bodies = buffer;
    if( key[1] != numScatter ) sort(bodies,buffer);
    delete[] scnt;
    delete[] sdsp;
  }

  void bisectionGather(Bodies &bodies, int nthLocal, int &newSize) {
    int const bytes = sizeof(bodies[0]);
    int numGather = nprocs[0] - 1;
    int oldSize = newSize;
    int scnt;
    int *rcnt = new int [nprocs[0]];
    int *rdsp = new int [nprocs[0]];
    if( key[0] != 0 ) {
      scnt = nthLocal / numGather;
      newSize -= scnt;
       buffer.erase( buffer.begin(), buffer.begin()+scnt);
    }
    MPI_Gather(&scnt,1,MPI_INT,rcnt,1,MPI_INT,0,MPI_COMM[0]);
    if( key[0] == 0 ) {
      rdsp[0] = 0;
      for(int i=0; i!=numGather; ++i ) {
        rdsp[i+1] = rdsp[i] + rcnt[i];
      }
      newSize += rdsp[numGather] + rcnt[numGather];
       buffer.resize(newSize);
    }

    scnt *= bytes;
    for(int i=0; i!= nprocs[0]; ++i ) {
      rcnt[i] *= bytes;
      rdsp[i] *= bytes;
    }

    MPI_Gatherv(&bodies[0].I,      scnt,     MPI_BYTE,
                &buffer[oldSize].I,rcnt,rdsp,MPI_BYTE,
                0,MPI_COMM[0]);
    bodies = buffer;
    delete[] rcnt;
    delete[] rdsp;
    if( key[0] == 0 ) sort(bodies,buffer);
  }

  void bisection(Bodies &bodies) {
    int const MPI_TYPE = getType(bodies[0].I);
    nprocs[0] = nprocs[1] = SIZE;                               // Initialize number of processes in groups
    offset[0] = offset[1] = 0;                                  // Initialize offset of body in groups
     color[0] =  color[1] =  color[2] = 0;                      // Initialize color of communicators
       key[0] =    key[1] =    key[2] = 0;                      // Initialize key of communicators
    int newSize;
    bigint numLocal = bodies.size();
    bigint numGlobal;
    MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM_WORLD);
    bigint nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];
    binBodies(bodies,0);
    buffer.resize(numLocal);
    sort(bodies,buffer);
    bigint iSplit = nth_element(bodies,nthGlobal);
    int nthLocal = splitBodies(bodies,iSplit);

    for( int l=0; l!=LEVEL; ++l ) {
      bisectionGetComm(l);

      splitDomain(iSplit,l,l%3);

      bisectionAlltoall(bodies,nthLocal,numLocal,newSize);

      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] <= nprocs[1] )
        bisectionScatter(bodies,nthLocal,newSize);
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] >= nprocs[1] )
        bisectionGather(bodies,nthLocal,newSize);

#ifdef DEBUG
      for(int j=0; j!=SIZE; ++j) {
        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100);
        if( RANK == j ) {
          std::cout << "RANK : " << j << std::endl;
          for(int i=0; i!=9; ++i) std::cout << bodies[numLocal/10*i].I << " ";
          std::cout << std::endl;
        }
      }
#endif
      numLocal = newSize;
      MPI_Allreduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,MPI_COMM[0]);
      nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];
      binBodies(bodies,(l+1) % 3);
      buffer.resize(numLocal);
      sort(bodies,buffer);
      iSplit = nth_element(bodies,nthGlobal,MPI_COMM[0]);
      nthLocal = splitBodies(bodies,iSplit);
    }
  }

  void octsection(Bodies &bodies) {
    int byte = sizeof(bodies[0]);
    int level = int(log(SIZE + 1) / M_LN2 / 3);
    BottomUp::setIndex(bodies,level);
    buffer.resize(bodies.size());
    sort(bodies,buffer);
    int *scnt = new int [SIZE];
    int *sdsp = new int [SIZE];
    int *rcnt = new int [SIZE];
    int *rdsp = new int [SIZE];
    for( int i=0; i!=SIZE; ++i ) {
      scnt[i] = 0;
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B->I - ((1 << 3*level) - 1) / 7;
      scnt[i]++;
    }
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD);
    sdsp[0] = rdsp[0] = 0;
    for( int irank=0; irank!=SIZE-1; ++irank ) {
      sdsp[irank+1] = sdsp[irank] + scnt[irank];
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];
    }
    buffer.resize(rdsp[SIZE-1]+rcnt[SIZE-1]);
    for( int irank=0; irank!=SIZE; ++irank ) {
      scnt[irank] *= byte;
      sdsp[irank] *= byte;
      rcnt[irank] *= byte;
      rdsp[irank] *= byte;
    }
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);
    bodies = buffer;
    XMIN[0] = X0-R0;
    XMAX[0] = X0+R0;
    for( int l=0; l!=level; ++l ) {
      for( int d=0; d!=3; ++d ) {
        int i = 3 * l + d;
        XMIN[i+1] = XMIN[i];
        XMAX[i+1] = XMAX[i];
        if( (RANK >> (3 * (level - l - 1) + 2 - d)) % 2 ) {
          XMIN[i+1][2-d] = (XMAX[i][2-d]+XMIN[i][2-d]) / 2;
        } else {
          XMAX[i+1][2-d] = (XMAX[i][2-d]+XMIN[i][2-d]) / 2;
        }
      }
    }
    delete[] scnt;
    delete[] sdsp;
    delete[] rcnt;
    delete[] rdsp;
  }
};

#endif
