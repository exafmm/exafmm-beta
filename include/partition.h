#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "types.h"
#include "sort.h"

class Partition : public MyMPI, public Sort {
public:
  bigint nth_element(bigint *data, int numData, bigint n) {     // Find nth element of global data
    int maxBucket(1000);                                        // Maximum number of buckets
    int numBucket;                                              // Number of buckets
    int numSample(std::min(maxBucket/SIZE,numData));            // Number of local samples
    int lOffset(0);                                             // Local offset of region being considered
    bigint gOffset(0);                                          // Global offset of region being considered
    std::vector<int>    rcnt(SIZE);                             // MPI recieve count
    std::vector<int>    rdsp(SIZE);                             // MPI recieve displacement
    std::vector<bigint>  send(maxBucket);                       // MPI send buffer for data
    std::vector<bigint>  recv(maxBucket);                       // MPI recieve buffer for data
    std::vector<bigint> isend(maxBucket);                       // MPI send buffer for index
    std::vector<bigint> irecv(maxBucket);                       // MPI recieve buffer for index
    std::vector<bigint> iredu(maxBucket);                       // Local scan of index buffer
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
        rdsp[irank] = numBucket;                                //   Calculate recieve displacement
        numBucket += rcnt[irank];                               //   Accumulate recieve count to get number of buckets
      }                                                         //  End loop over all processes
    }                                                           // Endif for rank 0
    MPI_Gatherv(&send[0],numSample,MPI_LONG,                    // Gather sample data to rank 0
                &recv[0],&rcnt[0],&rdsp[0],MPI_LONG,
                0,MPI_COMM_WORLD);
    if( RANK == 0 ) {                                           // Only rank 0 operates on gathered info
      recv.resize(numBucket);                                   //  Resize recv so that end() is valid
      std::sort(recv.begin(),recv.end());                       //  Sort the bucket data
      numBucket = std::unique(recv.begin(),recv.end()) - recv.begin(); // Remove duplicate bucket data
      recv.resize(numBucket);                                   //  Resize recv again
    }                                                           // Endif for rank 0
    MPI_Bcast(&numBucket,1,MPI_INT,0,MPI_COMM_WORLD);           // Broadcast number of buckets

    while( numBucket > 1 ) {                                    // While there are multipole candidates
      MPI_Bcast(&recv[0],numBucket,MPI_LONG,0,MPI_COMM_WORLD);  //  Broadcast bucket data
      int ic(0),nth(0);                                         //  Initialize counters
      std::fill(isend.begin(),isend.end(),0);                   //  Initialize bucket counter
      for( int i=0; i!=numData; ++i ) {                         //  Loop over range of data
        while( data[lOffset + i] > recv[ic] && ic < numBucket-1 ) ++ic;// Set counter to current bucket
        isend[ic]++;                                            //   Increment bucket counter
      }                                                         //  End loop over data
      MPI_Reduce(&isend[0],&irecv[0],numBucket,MPI_LONG,        //  Reduce bucket counter
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
      MPI_Bcast(&gOffset,1,MPI_LONG,0,MPI_COMM_WORLD);          //  Broadcast global offset
      iredu[0] = 0;                                             //  Initialize local scan index
      for( int i=0; i!=numBucket-1; ++i )                       //  Loop over buckets
        iredu[i+1] = iredu[i] + isend[i];                       //   Increment local scan index
      if( nth == numBucket-1 )                                  //  If nth is last bucket
        numData = numData-iredu[nth];                           //   Region of interest is that bucket
      else                                                      //  If nth is not the last bucket
        numData = iredu[nth+1]-iredu[nth];                      //   Region of interest is that bucket
      lOffset += iredu[nth];                                    //  Increment local offset to that bucket
      numSample = std::min(maxBucket/SIZE,numData);             //  Update number of local samples
      for( int i=0; i!=numSample; ++i ) {                       //  Loop over local samples
        int stride = numData/numSample;                         //   Sampling stride
        send[i] = data[lOffset + i * stride];                   //   Put sampled data in send buffer
      }                                                         //  End loop over samples
      MPI_Gather(&numSample,1,MPI_INT,                          //  Gather sample data to rank 0
                 &rcnt[0],  1,MPI_INT,
                 0,MPI_COMM_WORLD);
      if( RANK == 0 ) {                                         //  Only rank 0 operates on gathered info
        numBucket = 0;                                          //   Initialize number of buckets
        for( int irank=0; irank!=SIZE; ++irank ) {              //   Loop over all processes
          rdsp[irank] = numBucket;                              //    Calculate recieve displacement
          numBucket += rcnt[irank];                             //    Accumulate recieve count to get number of buckets
        }                                                       //   End loop over all processes
      }                                                         //  Endif for rank 0
      recv.resize(numBucket);                                   //  Resize recv
      MPI_Gatherv(&send[0],numSample,MPI_LONG,                  //  Gather sample data to rank 0
                  &recv[0],&rcnt[0],&rdsp[0],MPI_LONG,
                  0,MPI_COMM_WORLD);
      if( RANK == 0 ) {                                         //  Only rank 0 operates on gathered info
        recv.resize(numBucket);                                 //   Resize recv so that end() is valid
        std::sort(recv.begin(),recv.end());                     //   Sort the bucket data
        numBucket = std::unique(recv.begin(),recv.end()) - recv.begin();// Remove duplicate bucket data
        recv.resize(numBucket);                                 //   Resize recv again
      }                                                         //  Endif for rank 0
      MPI_Bcast(&numBucket,1,MPI_INT,0,MPI_COMM_WORLD);         //  Broadcast number of buckets
    }                                                           // End while loop
    MPI_Bcast(&recv[0],1,MPI_INT,0,MPI_COMM_WORLD);             // Broadcast nth element
    return recv[0];                                             // Return nth element
  }
  
};

#endif
