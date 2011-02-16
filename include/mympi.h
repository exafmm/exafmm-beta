#ifndef mympi_h
#define mympi_h
#include <mpi.h>
#include <typeinfo>
#include "types.h"

class MyMPI {                                                   // My own MPI utilities
private:
  int const WAIT;                                               // Waiting time between output of different ranks
protected:
  int       SIZE;                                               // Number of MPI processes
  int       RANK;                                               // Index of current MPI process
  int       SIZES;                                              // Number of MPI processes for split communicator
  int       RANKS;                                              // Index of current MPI process for split communicator
public:
  MyMPI() : WAIT(100) {                                         // Constructor, initialize WAIT time
    int argc(0);                                                // Dummy argument count
    char **argv;                                                // Dummy argument value
    MPI_Init(&argc,&argv);                                      // Initialize MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD,&SIZE);                        // Get number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);                        // Get index of current MPI process
    MPISIZE = SIZE;                                             // Set global variable MPISIZE
    MPIRANK = RANK;                                             // Set global variable MPIRANK
  }

  ~MyMPI() {                                                    // Destructor
    MPI_Finalize();                                             // Finalize MPI communicator
  }

  int commSize() { return SIZE; }                               // Number of MPI processes
  int commRank() { return RANK; }                               // Index of current MPI process

  bool isPowerOfTwo(int const n) {                              // If n is power of two return true
    return ((n != 0) && !(n & (n - 1)));                        // Decrement and compare bits
  }

  void splitRange(int &begin, int &end, int iSplit, int numSplit) {
    int size = end - begin;
    int increment = size / numSplit;
    int remainder = size % numSplit;
    begin += iSplit * increment + std::min(iSplit,remainder);
    end = begin + increment;
    if( remainder > iSplit ) end++;
  }

  template<typename T>
  int getType(T object) {
    int type;
    if       ( typeid(object) == typeid(char) ) {
      type = MPI_CHAR;
    } else if( typeid(object) == typeid(short) ) {
      type = MPI_SHORT;
    } else if( typeid(object) == typeid(int) ) {
      type = MPI_INT;
    } else if( typeid(object) == typeid(long) ) {
      type = MPI_LONG;
    } else if( typeid(object) == typeid(long long) ) {
      type = MPI_LONG_LONG;
    } else if( typeid(object) == typeid(unsigned char) ) {
      type = MPI_UNSIGNED_CHAR;
    } else if( typeid(object) == typeid(unsigned short) ) {
      type = MPI_UNSIGNED_SHORT;
    } else if( typeid(object) == typeid(unsigned int) ) {
      type = MPI_UNSIGNED;
    } else if( typeid(object) == typeid(unsigned long) ) {
      type = MPI_UNSIGNED_LONG;
    } else if( typeid(object) == typeid(unsigned long long) ) {
      type = MPI_UNSIGNED_LONG_LONG;
    } else if( typeid(object) == typeid(float) ) {
      type = MPI_FLOAT;
    } else if( typeid(object) == typeid(double) ) {
      type = MPI_DOUBLE;
    } else if( typeid(object) == typeid(long double) ) {
      type = MPI_LONG_DOUBLE;
    } else if( typeid(object) == typeid(std::complex<float>) ) {
      type = MPI_COMPLEX;
    } else if( typeid(object) == typeid(std::complex<double>) ) {
      type = MPI_DOUBLE_COMPLEX;
    }
    return type;
  }

  template<typename T>                                          // Detect data type to print
  void print(T data) {                                          // Print a scalar value on all ranks
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
      usleep(WAIT);                                             //  Wait "WAIT" milliseconds
      if( RANK == irank ) std::cout << data << " ";             //  If it's my turn print "data"
    }                                                           // End loop over ranks
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( RANK == 0 ) std::cout << std::endl;                     // New lint
  }

  template<typename T>
  void print(T data, int const irank) {
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( RANK == irank ) std::cout << data;
  }

  template<typename T>
  void print(T *data, int const begin, int const end) {
    for( int irank=0; irank!=SIZE; ++irank ) {
      MPI_Barrier(MPI_COMM_WORLD);
      usleep(WAIT);
      if( RANK == irank ) {
        std::cout << RANK << " : ";
        for( int i=begin; i!=end; ++i ) {
          std::cout << data[i] << " ";
        }
        std::cout << std::endl;
      }
    }
  }

  template<typename T>
  void print(T *data, int const begin, int const end, int const irank) {
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( RANK == irank ) {
      std::cout << RANK << " : ";
      for( int i=begin; i!=end; ++i ) {
        std::cout << data[i] << " ";
      }
      std::cout << std::endl;
    }
  }

};

#endif
