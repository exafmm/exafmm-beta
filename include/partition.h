#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "types.h"
#include "sort.h"

class Partition : public MyMPI, public Sort {
private:
  Bodies &bodies;                                               // Bodies to be partitioned
  vect XMIN;                                                    // Minimum position vector of bodies
  vect XMAX;                                                    // Maximum position vector of bodies
public:
  Partition(Bodies &b) : bodies(b) {}                           // Constructor
  ~Partition() {}                                               // Destructor

  void setDomain() {                                            // Set bounds of domain to be partitioned
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN = XMAX = B->pos;                                       // Initialize xmin,xmax
    int const MPI_TYPE = getType(XMIN[0]);                      // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < XMIN[d]) XMIN[d] = B->pos[d];       //   Determine xmin
        else if(B->pos[d] > XMAX[d]) XMAX[d] = B->pos[d];       //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over all bodies
    vect buffer;                                                // Define recv buffer
    MPI_Reduce(XMAX,buffer,3,MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD);// Reduce global maximum
    XMAX = buffer;                                              // Get data from buffer
    MPI_Bcast(XMAX,3,MPI_TYPE,0,MPI_COMM_WORLD);                // Broadcast global maximum
    MPI_Reduce(XMIN,buffer,3,MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD);// Reduce global minimum
    XMIN = buffer;                                              // Get data from buffer
    MPI_Bcast(XMIN,3,MPI_TYPE,0,MPI_COMM_WORLD);                // Broadcast global minimum
  }

  void binPosition(bigint *Ibody, int d) {                      // Turn positions into indices of bins
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      Ibody[B-bodies.begin()] = bigint((B->pos[d] - XMIN[d]) / (XMAX[d] - XMIN[d]) * bodies.size());// Bin bodies
    }                                                           // End loop over all bodies
  }

  template<typename T>
  int getBucket(T *data, int numData, int lOffset, std::vector<T> &send, std::vector<T> &recv) {
    int maxBucket(send.size());                                 // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample(std::min(maxBucket/SIZE,numData));            // Number of local samples
    int const MPI_TYPE = getType(data[0]);                      // Get MPI data type
    std::vector<int> rcnt(SIZE);                                // MPI receive count
    std::vector<int> rdsp(SIZE);                                // MPI receive displacement
    for( int i=0; i!=numSample; ++i ) {                         // Loop over local samples
      int stride = numData/numSample;                           //  Sampling stride
      send[i] = data[lOffset + i * stride];                     //  Put sampled data in send buffer
    }                                                           // End loop over samples
    MPI_Gather(&numSample,1,MPI_INT,                            // Gather size of sample data to rank 0
               &rcnt[0],  1,MPI_INT,
               0,MPI_COMM_WORLD);
    if( RANK == 0 ) {                                           // Only rank 0 operates on gathered info
      numBucket = 0;                                            //  Initialize number of buckets
      for( int irank=0; irank!=SIZE; ++irank ) {                //  Loop over all processes
        rdsp[irank] = numBucket;                                //   Calculate receive displacement
        numBucket += rcnt[irank];                               //   Accumulate receive count to get number of buckets
      }                                                         //  End loop over all processes
      recv.resize(numBucket);                                   //  Resize recv so that end() is valid
    }                                                           // Endif for rank 0
    MPI_Gatherv(&send[0],numSample,MPI_TYPE,                    // Gather sample data to rank 0
                &recv[0],&rcnt[0],&rdsp[0],MPI_TYPE,
                0,MPI_COMM_WORLD);
    if( RANK == 0 ) {                                           // Only rank 0 operates on gathered info
      std::sort(recv.begin(),recv.end());                       //  Sort the bucket data
      numBucket = std::unique(recv.begin(),recv.end()) - recv.begin();// Remove duplicate bucket data
      recv.resize(numBucket);                                   //  Resize recv again
    }                                                           // Endif for rank 0
    MPI_Bcast(&numBucket,1,MPI_INT,0,MPI_COMM_WORLD);           // Broadcast number of buckets
    MPI_Bcast(&recv[0],numBucket,MPI_TYPE,0,MPI_COMM_WORLD);    // Broadcast bucket data
    return numBucket;                                           // Return number of buckets
  }

  template<typename T, typename T2>
  T2 nth_element(T *data, int numData, T2 n) {                  // Find nth element of global data
    int maxBucket(1000);                                        // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int lOffset(0);                                             // Local offset of region being considered
    int const MPI_TYPE = getType(n);                            // Get MPI data type
    std::vector<int>    rcnt(SIZE);                             // MPI receive count
    T2 gOffset(0);                                              // Global offset of region being considered
    std::vector<T>   send(maxBucket);                           // MPI send buffer for data
    std::vector<T>   recv(maxBucket);                           // MPI receive buffer for data
    std::vector<T2> isend(maxBucket);                           // MPI send buffer for index
    std::vector<T2> irecv(maxBucket);                           // MPI receive buffer for index
    std::vector<T2> iredu(maxBucket);                           // Local scan of index buffer
    numBucket = getBucket(data,numData,lOffset,send,recv);      // Get global bucket data

    while( numBucket > 1 ) {                                    // While there are multipole candidates
      int ic(0),nth(0);                                         //  Initialize counters
      std::fill(isend.begin(),isend.end(),0);                   //  Initialize bucket counter
      for( int i=0; i!=numData; ++i ) {                         //  Loop over range of data
        while( data[lOffset + i] > recv[ic] && ic < numBucket-1 ) ++ic;// Set counter to current bucket
        isend[ic]++;                                            //   Increment bucket counter
      }                                                         //  End loop over data
      MPI_Reduce(&isend[0],&irecv[0],numBucket,MPI_TYPE,        //  Reduce bucket counter
                 MPI_SUM,0,MPI_COMM_WORLD);
      if( RANK == 0 ) {                                         //  Only rank 0 operates on reduced data
        iredu[0] = 0;                                           //   Initialize global scan index
        for( int i=0; i!=numBucket-1; ++i )                     //   Loop over buckets
          iredu[i+1] = iredu[i] + irecv[i];                     //    Increment global scan index
        nth = 0;                                                //   Initialize index for bucket containing nth element
        while( n - gOffset > iredu[nth] && nth < numBucket ) ++nth;// Find index for bucket containing nth element
        --nth;                                                  //   Account for overshoot (do while?)
        gOffset += iredu[nth];                                  //   Increment global offset of region being considered
      }                                                         //  Endif for rank 0
      MPI_Bcast(&nth,1,MPI_INT,0,MPI_COMM_WORLD);               //  Broadcast index for bucket containing nth element
      MPI_Bcast(&gOffset,1,MPI_TYPE,0,MPI_COMM_WORLD);          //  Broadcast global offset
      iredu[0] = 0;                                             //  Initialize local scan index
      for( int i=0; i!=numBucket-1; ++i )                       //  Loop over buckets
        iredu[i+1] = iredu[i] + isend[i];                       //   Increment local scan index
      if( nth == numBucket-1 )                                  //  If nth is last bucket
        numData = numData-iredu[nth];                           //   Region of interest is that bucket
      else                                                      //  If nth is not the last bucket
        numData = iredu[nth+1]-iredu[nth];                      //   Region of interest is that bucket
      lOffset += iredu[nth];                                    //  Increment local offset to that bucket
      numBucket = getBucket(data,numData,lOffset,send,recv);    // Get global bucket data
    }                                                           // End while loop
    return recv[0];                                             // Return nth element
  }

  void shiftBodies(Bodies &buffer) {                            // Send bodies to next rank (wrap around)
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
    MPI_Isend(&bodies[0].pos[0],oldSize*bytes,MPI_CHAR,irecv,   // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0].pos[0],newSize*bytes,MPI_CHAR,isend,   // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

};

#endif
