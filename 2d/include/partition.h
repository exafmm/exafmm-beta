#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "logger.h"
#include "types.h"

//! Handles all the partitioning of domains
class Partition : public MyMPI, public Logger {
 protected:
  Bodies sendBodies;                                            //!< Send buffer for bodies
  Bodies recvBodies;                                            //!< Receive buffer for bodies
  int * sendBodyCount;                                          //!< Send count
  int * sendBodyDispl;                                          //!< Send displacement
  int * recvBodyCount;                                          //!< Receive count
  int * recvBodyDispl;                                          //!< Receive displacement

 protected:
//! Exchange send count for bodies
  void alltoall(Bodies &bodies) {
    for (int i=0; i<mpisize; i++) {                             // Loop over ranks
      sendBodyCount[i] = 0;                                     //  Initialize send counts
    }                                                           // End loop over ranks
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      assert(0 <= B->IPROC && B->IPROC < mpisize);
      sendBodyCount[B->IPROC]++;                                //  Fill send count bucket
      B->IPROC = mpirank;                                       //  Tag for sending back to original rank
    }                                                           // End loop over bodies
    MPI_Alltoall(sendBodyCount, 1, MPI_INT,                     // Communicate send count to get receive count
                 recvBodyCount, 1, MPI_INT, MPI_COMM_WORLD);
    sendBodyDispl[0] = recvBodyDispl[0] = 0;                    // Initialize send/receive displacements
    for (int irank=0; irank<mpisize-1; irank++) {               // Loop over ranks
      sendBodyDispl[irank+1] = sendBodyDispl[irank] + sendBodyCount[irank];//  Set send displacement
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

//! Exchange bodies
  void alltoallv(Bodies &bodies) {
    int word = sizeof(bodies[0]) / 4;                           // Word size of body structure
    recvBodies.resize(recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1]);// Resize receive buffer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendBodyCount[irank] *= word;                             //  Multiply send count by word size of data
      sendBodyDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvBodyCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0], sendBodyCount, sendBodyDispl, MPI_INT,// Communicate bodies
                  &recvBodies[0], recvBodyCount, recvBodyDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendBodyCount[irank] /= word;                             //  Divide send count by word size of data
      sendBodyDispl[irank] /= word;                             //  Divide send displacement by word size of data
      recvBodyCount[irank] /= word;                             //  Divide receive count by word size of data
      recvBodyDispl[irank] /= word;                             //  Divide receive displacement by word size of data
    }                                                           // End loop over ranks
  }

 public:
//! Constructor
  Partition() {
    sendBodyCount = new int [mpisize];                          // Allocate send count
    sendBodyDispl = new int [mpisize];                          // Allocate send displacement
    recvBodyCount = new int [mpisize];                          // Allocate receive count
    recvBodyDispl = new int [mpisize];                          // Allocate receive displacement
  }
//! Destructor
  ~Partition() {
    delete[] sendBodyCount;                                     // Deallocate send count
    delete[] sendBodyDispl;                                     // Deallocate send displacement
    delete[] recvBodyCount;                                     // Deallocate receive count
    delete[] recvBodyDispl;                                     // Deallocate receive displacement
  }

//! Send bodies to next rank (round robin)
  void shiftBodies(Bodies &bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int word = sizeof(bodies[0]) / 4;                     // Word size of body structure
    const int isend = (mpirank + 1          ) % mpisize;        // Send to next rank (wrap around)
    const int irecv = (mpirank - 1 + mpisize) % mpisize;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);// Send current number of bodies
    MPI_Irecv(&newSize, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);// Receive new number of bodies
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete

    recvBodies.resize(newSize);                                 // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0], oldSize*word, MPI_INT, irecv,         // Send bodies to next rank
              1, MPI_COMM_WORLD, &sreq);
    MPI_Irecv(&recvBodies[0], newSize*word, MPI_INT, isend,     // Receive bodies from previous rank
              1, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete
    bodies = recvBodies;                                        // Copy bodies from buffer
  }

//! Allgather bodies
  Bodies allgatherBodies(Bodies &bodies) {
    const int word = sizeof(bodies[0]) / 4;                     // Word size of body structure
    sendBodyCount[0] = bodies.size();                           // Determine send count
    MPI_Allgather(sendBodyCount, word, MPI_INT,                 // Allgather number of bodies
                  recvBodyCount, word, MPI_INT, MPI_COMM_WORLD);
    recvBodyDispl[0] = 0;                                       // Initialize receive displacement
    for (int irank=0; irank<mpisize-1; irank++) {               // Loop over ranks
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];// Set receive displacement
    }                                                           // End loop over ranks
    recvBodies.resize(recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1]);// Resize receive buffer
    recvBodies.resize(recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1]);// Resize receive buffer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      recvBodyCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Allgatherv(&bodies[0], sendBodyCount[0]*word, MPI_FLOAT,// Allgather bodies
                   &recvBodies[0], recvBodyCount, recvBodyDispl, MPI_FLOAT, MPI_COMM_WORLD);
    return recvBodies;                                          // Return bodies
  }

//! Allreduce int from all ranks
  int allreduceInt(int send) {
    int recv;                                                   // Receive buffer
    MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);// Communicate values
    return recv;                                                // Return received values
  }

//! Allreduce fvec2 from all ranks
  fvec2 allreduceVec2(fvec2 send) {
    fvec2 recv;                                                 // Receive buffer
    MPI_Allreduce(send, recv, 2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Communicate values
    return recv;                                                // Return received values
  }

//! Allreduce bounds from all ranks
  Bounds allreduceBounds(Bounds local) {
    fvec2 localXmin, localXmax, globalXmin, globalXmax;
    for (int d=0; d<2; d++) {                                   // Loop over dimensions
      localXmin[d] = local.Xmin[d];                             //  Convert Xmin to float
      localXmax[d] = local.Xmax[d];                             //  Convert Xmax to float
    }                                                           // End loop over dimensions
    MPI_Allreduce(localXmin, globalXmin, 2, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);// Reduce domain Xmin
    MPI_Allreduce(localXmax, globalXmax, 2, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmax
    Bounds global;
    for (int d=0; d<2; d++) {                                   // Loop over dimensions
      real_t leeway = (globalXmax[d] - globalXmin[d]) * 1e-6;   //  Adding a bit of leeway to global domain
      global.Xmin[d] = globalXmin[d] - leeway;                  //  Convert Xmin to real_t
      global.Xmax[d] = globalXmax[d] + leeway;                  //  Convert Xmax to real_t
    }                                                           // End loop over dimensions
    return global;                                              // Return global bounds
  }

//! Partition bodies
  Bounds partition(Bodies &bodies, Bounds global) {
    startTimer("Partition");                                    // Start timer
    int size = mpisize;                                         // Initialize MPI size counter
    vec<2,int> Npartition = 1;                                  // Number of partitions in each direction
    int d = 0;                                                  // Initialize dimension counter
    while (size != 1) {                                         // Divide domain while counter is not one
      Npartition[d] <<= 1;                                      //  Divide this dimension
      d = (d+1) % 2;                                            //  Increment dimension
      size >>= 1;                                               //  Right shift the bits of counter
    }                                                           // End while loop for domain subdivision
    vec2 Xpartition;                                            // Size of partitions in each direction
    for (d=0; d<2; d++) {                                       // Loop over dimensions
      Xpartition[d] = (global.Xmax[d] - global.Xmin[d]) / Npartition[d];//  Size of partition in each direction
    }                                                           // End loop over dimensions
    int ix[2];                                                  // Index vector
    ix[0] = mpirank % Npartition[0];                            // x index of partition
    ix[1] = mpirank / Npartition[0];                            // y index
    Bounds local;                                               // Local bounds
    for (d=0; d<2; d++) {                                       // Loop over dimensions
      local.Xmin[d] = global.Xmin[d] + ix[d] * Xpartition[d];   // Xmin of local domain at current rank
      local.Xmax[d] = global.Xmin[d] + (ix[d] + 1) * Xpartition[d];// Xmax of local domain at current rank
    }                                                           // End loop over dimensions
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      for (d=0; d<2; d++) {                                     //  Loop over dimensions
        ix[d] = int((B->X[d] - global.Xmin[d]) / Xpartition[d]);//   Index vector of partition
      }                                                         //  End loop over dimensions
      B->IPROC = ix[0] + Npartition[0] * ix[1];                 //  Set send rank
      assert(0 <= B->IPROC && B->IPROC < mpisize);
      B->ICELL = B->IPROC;                                      //  Do this to sort accroding to IPROC
    }                                                           // End loop over bodies
    stopTimer("Partition");                                     // Stop timer
    return local;
  }

//! Send bodies back to where they came from
  void unpartition(Bodies &bodies) {
    startTimer("Unpartition");                                  // Start timer
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->ICELL = B->IPROC;                                      //  Do this to sortaccroding to IPROC
    }                                                           // End loop over bodies
    stopTimer("Unpartition");                                   // Stop timer
  }

  //! Write send count to file
  void writeSendCount() {
    std::stringstream name;                                     // File name
    name << "send" << std::setfill('0') << std::setw(6)         // Set format
         << mpirank << ".dat";                                  // Create file name for timer
    std::ofstream sendCountFile(name.str().c_str(), std::ios::app); // Open sendCount log file
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendCountFile << std::setw(stringLength) << std::left     //  Set format
		    << sendBodyCount[irank] << std::endl;       //  Print send count
    }                                                           // End loop over ranks
    sendCountFile.close();
  }
};
#endif
