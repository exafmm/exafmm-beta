#ifndef base_mpi_h
#define base_mpi_h
#include <mpi.h>
#include <iostream>
#include "types.h"
#include <unistd.h>

namespace exafmm {
  //! Custom MPI utilities
  class BaseMPI {
  private:
    int external;                                               //!< Flag to indicate external MPI_Init/Finalize

  protected:
    const int wait;                                             //!< Waiting time between output of different ranks

  public:
    int mpirank;                                                //!< Rank of MPI communicator
    int mpisize;                                                //!< Size of MPI communicator

  public:
    //! Constructor
    BaseMPI() : external(0), wait(100) {                        // Initialize variables
      int argc(0);                                              // Dummy argument count
      char **argv;                                              // Dummy argument value
      MPI_Initialized(&external);                               // Check if MPI_Init has been called
      if (!external) MPI_Init(&argc, &argv);                    // Initialize MPI communicator
      MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);                  // Get rank of current MPI process
      MPI_Comm_size(MPI_COMM_WORLD, &mpisize);                  // Get number of MPI processes
    }

    //! Destructor
    ~BaseMPI() {
      if (!external) MPI_Finalize();                            // Finalize MPI communicator
    }

    //! Allreduce int type from all ranks
    int allreduceInt(int send) {
      int recv;                                                 // Receive buffer
      MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);// Communicate values
      return recv;                                              // Return received values
    }

    //! Allreduce vec3 type from all ranks
    vec3 allreduceVec3(vec3 send) {
      float fsend[3], frecv[3];                                 // Single precision buffers
      for (int d=0; d<3; d++) fsend[d] = send[d];               // Copy to send buffer
      MPI_Allreduce(fsend, frecv, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);// Communicate values
      vec3 recv;                                                // Receive buffer
      for (int d=0; d<3; d++) recv[d] = frecv[d];               // Copy from recv buffer
      return recv;                                              // Return received values
    }

    //! Allreduce bounds type from all ranks
    Bounds allreduceBounds(Bounds local) {
      float localXmin[3], localXmax[3], globalXmin[3], globalXmax[3];
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	localXmin[d] = local.Xmin[d];                           //  Convert Xmin to float
	localXmax[d] = local.Xmax[d];                           //  Convert Xmax to float
      }                                                         // End loop over dimensions
      MPI_Allreduce(localXmin, globalXmin, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);// Reduce domain Xmin
      MPI_Allreduce(localXmax, globalXmax, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmax
      Bounds global;
      for (int d=0; d<3; d++) {                                 // Loop over dimensions
	real_t leeway = (globalXmax[d] - globalXmin[d]) * 1e-6; //  Adding a bit of leeway to global domain
	global.Xmin[d] = globalXmin[d] - leeway;                //  Convert Xmin to real_t
	global.Xmax[d] = globalXmax[d] + leeway;                //  Convert Xmax to real_t
      }                                                         // End loop over dimensions
      return global;                                            // Return global bounds
    }

    //! Print a scalar value on all ranks
    template<typename T>
    void print(T data) {
      for (int irank=0; irank<mpisize; irank++ ) {              // Loop over ranks
	MPI_Barrier(MPI_COMM_WORLD);                            //  Sync processes
	usleep(wait);                                           //  Wait "wait" milliseconds
	if (mpirank == irank) std::cout << data << " ";         //  If it's my turn print "data"
      }                                                         // End loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              // Sync processes
      usleep(wait);                                             // Wait "wait" milliseconds
      if (mpirank == mpisize-1) std::cout << std::endl;         // New line
    }

    //! Print a scalar value on irank
    template<typename T>
    void print(T data, const int irank) {
      MPI_Barrier(MPI_COMM_WORLD);                              // Sync processes
      usleep(wait);                                             // Wait "wait" milliseconds
      if( mpirank == irank ) std::cout << data;                 // If it's my rank print "data"
    }

    //! Print a vector value on all ranks
    template<typename T>
    void print(T * data, const int begin, const int end) {
      for (int irank=0; irank<mpisize; irank++) {               // Loop over ranks
	MPI_Barrier(MPI_COMM_WORLD);                            //  Sync processes
	usleep(wait);                                           //  Wait "wait" milliseconds
	if (mpirank == irank) {                                 //  If it's my turn to print
	  std::cout << mpirank << " : ";                        //   Print rank
	  for (int i=begin; i<end; i++) {                       //   Loop over data
	    std::cout << data[i] << " ";                        //    Print data[i]
	  }                                                     //   End loop over data
	  std::cout << std::endl;                               //   New line
	}                                                       //  Endif for my turn
      }                                                         // End loop over ranks
    }

    //! Print a vector value on irank
    template<typename T>
    void print(T * data, const int begin, const int end, const int irank) {
      MPI_Barrier(MPI_COMM_WORLD);                              // Sync processes
      usleep(wait);                                             // Wait "wait" milliseconds
      if (mpirank == irank) {                                   // If it's my rank
	std::cout << mpirank << " : ";                          //  Print rank
	for (int i=begin; i<end; i++) {                         //  Loop over data
	  std::cout << data[i] << " ";                          //   Print data[i]
	}                                                       //  End loop over data
	std::cout << std::endl;                                 //  New line
      }                                                         // Endif for my rank
    }
  };
}
#endif
