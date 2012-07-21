#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "serialfmm.h"

//! Handles all the partitioning of domains
class Partition : public MyMPI, public SerialFMM {
protected:
  std::vector<vect> XMIN;                                       //!< Minimum position vector of bodies
  std::vector<vect> XMAX;                                       //!< Maximum position vector of bodies
  Bodies sendBodies;                                            //!< Send buffer for bodies
  Bodies recvBodies;                                            //!< Receive buffer for bodies
  int *sendBodyCount;                                           //!< Send count
  int *sendBodyDispl;                                           //!< Send displacement
  int *recvBodyCount;                                           //!< Receive count
  int *recvBodyDispl;                                           //!< Receive displacement

private:
//! Allgather bounds of all partitions
  void allgather(vect xmin, vect xmax) {
    MPI_Datatype MPI_TYPE = getType(xmin[0]);                   // Get MPI data type
    XMIN.resize(MPISIZE);                                       // Xmin of every local domain
    XMAX.resize(MPISIZE);                                       // Xmax of every local domain
    MPI_Allgather(xmin,3,MPI_TYPE,&XMIN[0],3,MPI_TYPE,MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(xmax,3,MPI_TYPE,&XMAX[0],3,MPI_TYPE,MPI_COMM_WORLD);// Gather all domain bounds
  }

//! Set partition of global domain
  void setPartition(Bodies &bodies) {
    startTimer("Partition");                                    // Start timer
    int mpisize = MPISIZE;                                      // Initialize MPI size counter
    vec<3,int> Npartition = 1;                                  // Number of partitions in each direction
    int d = 0;                                                  // Initialize dimension counter
    while( mpisize != 1 ) {                                     // Divide domain  while counter is not one
      Npartition[d] <<= 1;                                      //  Divide this dimension
      d = (d+1) % 3;                                            //  Increment dimension
      mpisize >>= 1;                                            //  Right shift the bits of counter
    }                                                           // End while loop for domain subdivision
    vect Xpartition;                                            // Size of partitions in each direction
    for( d=0; d!=3; ++d ) {                                     // Loop over dimensions
      Xpartition[d] = 2 * R0 / Npartition[d];                   //  Size of partition in each direction
    }                                                           // End loop over dimensions
    int ix = MPIRANK % Npartition[0];                           // x index of partition
    int iy = MPIRANK / Npartition[0] % Npartition[1];           // y index
    int iz = MPIRANK / Npartition[0] / Npartition[1];           // z index
    vect xmin, xmax;                                            // Local domain boundaries at current rank
    xmin[0] = X0[0] - R0 + ix * Xpartition[0];                  // xmin of local domain at current rank
    xmin[1] = X0[1] - R0 + iy * Xpartition[1];                  // ymin
    xmin[2] = X0[2] - R0 + iz * Xpartition[2];                  // zmin
    xmax[0] = X0[0] - R0 + (ix + 1) * Xpartition[0];            // xmax of local domain at current rank
    xmax[1] = X0[1] - R0 + (iy + 1) * Xpartition[1];            // ymax
    xmax[2] = X0[2] - R0 + (iz + 1) * Xpartition[2];            // zmax
    allgather(xmin,xmax);                                       // Allgather bounds of partitions
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      ix = int((B->X[0] + R0 - X0[0]) / Xpartition[0]);         //  x contribution to send rank
      iy = int((B->X[1] + R0 - X0[1]) / Xpartition[1]);         //  y contribution
      iz = int((B->X[2] + R0 - X0[2]) / Xpartition[2]);         //  z contribution
      B->IPROC = ix + Npartition[0] * (iy + iz * Npartition[1]);//  Set send rank
      B->ICELL = B->IPROC;                                      //  Do this to sort accroding to IPROC
    }                                                           // End loop over bodies
    Bodies buffer = bodies;                                     // Sort buffer
    stopTimer("Partition",printNow);                            // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order of ICELL
  }


protected:
//! Exchange send count for bodies
  void alltoall(Bodies &bodies) {
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendBodyCount[i] = 0;                                     //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      sendBodyCount[B->IPROC]++;                                //  Fill send count bucket
      B->IPROC = MPIRANK;                                       //  Tag for sending back to original rank
    }                                                           // End loop over bodies
    MPI_Alltoall(sendBodyCount,1,MPI_INT,                       // Communicate send count to get receive count
                 recvBodyCount,1,MPI_INT,MPI_COMM_WORLD);
    sendBodyDispl[0] = recvBodyDispl[0] = 0;                    // Initialize send/receive displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sendBodyDispl[irank+1] = sendBodyDispl[irank] + sendBodyCount[irank];//  Set send displacement
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

//! Exchange bodies
  void alltoallv(Bodies &bodies) {
    int word = sizeof(bodies[0]) / 4;                           // Word size of body structure
    recvBodies.resize(recvBodyDispl[MPISIZE-1]+recvBodyCount[MPISIZE-1]);// Resize receive buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      sendBodyCount[irank] *= word;                             //  Multiply send count by word size of data
      sendBodyDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvBodyCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],sendBodyCount,sendBodyDispl,MPI_INT,// Communicate bodies
                  &recvBodies[0],recvBodyCount,recvBodyDispl,MPI_INT,MPI_COMM_WORLD);
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      sendBodyCount[irank] /= word;                             //  Divide send count by word size of data
      sendBodyDispl[irank] /= word;                             //  Divide send displacement by word size of data
      recvBodyCount[irank] /= word;                             //  Divide receive count by word size of data
      recvBodyDispl[irank] /= word;                             //  Divide receive displacement by word size of data
    }                                                           // End loop over ranks
  }

public:
//! Constructor
  Partition() : MyMPI(), XMIN(), XMAX(),
                sendBodyCount(NULL), sendBodyDispl(NULL), recvBodyCount(NULL), recvBodyDispl(NULL) {
    sendBodyCount = new int [MPISIZE];                          // Allocate send count
    sendBodyDispl = new int [MPISIZE];                          // Allocate send displacement
    recvBodyCount = new int [MPISIZE];                          // Allocate receive count
    recvBodyDispl = new int [MPISIZE];                          // Allocate receive displacement
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
    const int isend = (MPIRANK + 1          ) % MPISIZE;        // Send to next rank (wrap around)
    const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete

    recvBodies.resize(newSize);                                 // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0],oldSize*word,MPI_INT,irecv,            // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&recvBodies[0],newSize*word,MPI_INT,isend,        // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    bodies = recvBodies;                                        // Copy bodies from buffer
  }

//! Partition bodies
  void partition(Bodies &bodies) {
    setBounds(bodies);                                          // Set global bounds
    X0 = 0;                                                     // Override center of global domain
    R0 = M_PI;                                                  // Override radius of global domain
    setPartition(bodies);                                       // Set partitioning strategy
    startTimer("Partition comm");                               // Start timer
    alltoall(bodies);                                           // Alltoall send count
    alltoallv(bodies);                                          // Alltoallv bodies
    bodies = recvBodies;                                        // Copy receive buffer to bodies
    stopTimer("Partition comm",printNow);                       // Stop timer 
  }

//! Send bodies back to where they came from
  void unpartition(Bodies &bodies) {
    startTimer("Unpartition");                                  // Start timer
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = B->IPROC;                                      //  Do this to sort accroding to IPROC
    }                                                           // End loop over bodies
    Bodies buffer = bodies;                                     // Resize sort buffer
    stopTimer("Unpartition",printNow);                          // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Unpartition comm");                             // Start timer
    alltoall(bodies);                                           // Alltoall send count
    alltoallv(bodies);                                          // Alltoallv bodies
    bodies = recvBodies;                                        // Copy receive buffer to bodies
    stopTimer("Unpartition comm",printNow);                     // Stop timer 
  }
};

#endif
