#ifndef mympi_h
#define mympi_h
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <unistd.h>

//! Custom MPI utilities
class MyMPI {
 private:
  int external;                                                 //!< Flag to indicate external MPI_Init/Finalize

 protected:
  const int wait;                                               //!< Waiting time between output of different ranks

 public:
  int mpirank;                                                  //!< Rank of MPI communicator
  int mpisize;                                                  //!< Size of MPI communicator

 public:
//! Constructor, initialize wait time
  MyMPI() : external(0), wait(100) {                            // Constructor, initialize wait time
    int argc(0);                                                // Dummy argument count
    char **argv;                                                // Dummy argument value
    MPI_Initialized(&external);                                 // Check if MPI_Init has been called
    if (!external) MPI_Init(&argc, &argv);                      // Initialize MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);                    // Get number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);                    // Get rank of current MPI process
  }

//! If n is power of two return true
  bool isPowerOfTwo(const int n) {
    return ((n != 0) && !(n & (n - 1)));                        // Decrement and compare bits
  }

//! Split range and return partial range
  void splitRange(int &begin, int &end, int iSplit, int numSplit) {
    int size = end - begin;                                     // Size of range
    int increment = size / numSplit;                            // Increment of splitting
    int remainder = size % numSplit;                            // Remainder of splitting
    begin += iSplit * increment + std::min(iSplit,remainder);   // Increment the begin counter
    end = begin + increment;                                    // Increment the end counter
    if (remainder > iSplit) end++;                              // Adjust the end counter for remainder
  }

//! Print a scalar value on all ranks
  template<typename T>
  void print(T data) {
    for (int irank=0; irank<mpisize; irank++ ) {                // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
      usleep(wait);                                             //  Wait "wait" milliseconds
      if (mpirank == irank) std::cout << data << " ";           //  If it's my turn print "data"
    }                                                           // End loop over ranks
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(wait);                                               // Wait "wait" milliseconds
    if (mpirank == mpisize-1) std::cout << std::endl;           // New line
  }

//! Print a scalar value on irank
  template<typename T>
  void print(T data, const int irank) {
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(wait);                                               // Wait "wait" milliseconds
    if( mpirank == irank ) std::cout << data;                   // If it's my rank print "data"
  }

//! Print a vector value on all ranks
  template<typename T>
  void print(T * data, const int begin, const int end) {
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
      usleep(wait);                                             //  Wait "wait" milliseconds
      if (mpirank == irank) {                                   //  If it's my turn to print
        std::cout << mpirank << " : ";                          //   Print rank
        for (int i=begin; i<end; i++) {                        //   Loop over data
          std::cout << data[i] << " ";                          //    Print data[i]
        }                                                       //   End loop over data
        std::cout << std::endl;                                 //   New line
      }                                                         //  Endif for my turn
    }                                                           // End loop over ranks
  }

//! Print a vector value on irank
  template<typename T>
  void print(T * data, const int begin, const int end, const int irank) {
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(wait);                                               // Wait "wait" milliseconds
    if (mpirank == irank) {                                     // If it's my rank
      std::cout << mpirank << " : ";                            //  Print rank
      for (int i=begin; i<end; i++) {                          //  Loop over data
        std::cout << data[i] << " ";                            //   Print data[i]
      }                                                         //  End loop over data
      std::cout << std::endl;                                   //  New line
    }                                                           // Endif for my rank
  }
};
#endif
