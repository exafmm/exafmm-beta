#ifndef partition_h
#define partition_h
#include "logger.h"
#include "sort.h"

//! Handles all the partitioning of domains
class Partition {
private:
  int mpirank;                                                  //!< Rank of MPI communicator
  int mpisize;                                                  //!< Size of MPI communicator

public:
  //! Constructor
  Partition() {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);                    // Get rank of current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);                    // Get number of MPI processes
  }

  //! Partition bodies with geometric octsection
  Bounds octsection(Bodies & bodies, Bounds global) {
    logger::startTimer("Partition");                            // Start timer
    int size = mpisize;                                         // Initialize MPI size counter
    vec<3,int> Npartition = 1;                                  // Number of partitions in each direction
    int d = 0;                                                  // Initialize dimension counter
    while (size != 1) {                                         // Divide domain while counter is not one
      Npartition[d] <<= 1;                                      //  Divide this dimension
      d = (d+1) % 3;                                            //  Increment dimension
      size >>= 1;                                               //  Right shift the bits of counter
    }                                                           // End while loop for domain subdivision
    vec3 Xpartition;                                            // Size of partitions in each direction
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
      B->IPROC = iX[0] + Npartition[0] * (iX[1] + iX[2] * Npartition[1]);//  Set send rank
      assert(0 <= B->IPROC && B->IPROC < mpisize);
    }                                                           // End loop over bodies
    logger::stopTimer("Partition");                             // Stop timer
    logger::startTimer("Sort");                                 // Start timer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.iproc(bodies);                                // Sort bodies according to IPROC
    logger::stopTimer("Sort");                                  // Stop timer
    return local;
  }

  //! Send bodies back to where they came from
  void unpartition(Bodies & bodies) {
    logger::startTimer("Sort");                                 // Start timer
    Sort sort;                                                  // Instantiate sort class
    bodies = sort.iproc(bodies);                                // Sort bodies according to IPROC
    logger::stopTimer("Sort");                                  // Stop timer
  }
};
#endif
