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
#ifndef mympi_h
#define mympi_h
#include <mpi.h>
#include <typeinfo>
#include "types.h"

//! Custom MPI utilities
class MyMPI {
protected:
  const int WAIT;                                               //!< Waiting time between output of different ranks
  int       MPISIZES;                                           //!< Number of MPI processes for split communicator
  int       MPIRANKS;                                           //!< Rank of current MPI process for split communicator
public:
//! Constructor, initialize WAIT time
  MyMPI() : WAIT(100) {                                         // Constructor, initialize WAIT time
    int argc(0);                                                // Dummy argument count
    char **argv;                                                // Dummy argument value
    int ExternalMPI;                                            // Flag for external MPI call
    MPI_Initialized(&ExternalMPI);                              // Check if MPI_Init has been called
    if(!ExternalMPI) MPI_Init(&argc,&argv);                     // Initialize MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD,&MPISIZE);                     // Get number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD,&MPIRANK);                     // Get rank of current MPI process
    DEVICE = MPIRANK % GPUS;                                    // Get GPU device ID from MPI rank
  }

//! Destructor
  ~MyMPI() {
    int ExternalMPI;                                            // Flag for external MPI call
    MPI_Initialized(&ExternalMPI);                              // Check if MPI_Init has been called
    if(!ExternalMPI) MPI_Finalize();                            // Finalize MPI communicator
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
    if( remainder > iSplit ) end++;                             // Adjust the end counter for remainder
  }

//! Get MPI data type
  template<typename T>
  MPI_Datatype getType(T object) {
    MPI_Datatype type = MPI_BYTE;                               // MPI data type
    if       ( typeid(object) == typeid(char) ) {               // If data type is char
      type = MPI_CHAR;                                          //  use MPI_CHAR
    } else if( typeid(object) == typeid(short) ) {              // If data type is short
      type = MPI_SHORT;                                         //  use MPI_SHORT
    } else if( typeid(object) == typeid(int) ) {                // If data type is int
      type = MPI_INT;                                           //  use MPI_INT
    } else if( typeid(object) == typeid(long) ) {               // If data type is long
      type = MPI_LONG;                                          //  use MPI_LONG
    } else if( typeid(object) == typeid(long long) ) {          // If data type is long long
      type = MPI_LONG_LONG;                                     //  use MPI_LONG_LONG
    } else if( typeid(object) == typeid(unsigned char) ) {      // If data type is unsigned char
      type = MPI_UNSIGNED_CHAR;                                 //  use MPI_UNSIGNED_CHAR
    } else if( typeid(object) == typeid(unsigned short) ) {     // If data type is unsigned short
      type = MPI_UNSIGNED_SHORT;                                //  use MPI_UNSIGNED_SHORT
    } else if( typeid(object) == typeid(unsigned int) ) {       // If data type is unsigned int
      type = MPI_UNSIGNED;                                      //  use MPI_UNSIGNED
    } else if( typeid(object) == typeid(unsigned long) ) {      // If data type is unsigned long
      type = MPI_UNSIGNED_LONG;                                 //  use MPI_UNSIGNED_LONG
    } else if( typeid(object) == typeid(unsigned long long) ) { // If data type is unsigned long long
      type = MPI_UNSIGNED_LONG_LONG;                            //  use MPI_UNSIGNED_LONG_LONG
    } else if( typeid(object) == typeid(float) ) {              // If data type is float
      type = MPI_FLOAT;                                         //  use MPI_FLOAT
    } else if( typeid(object) == typeid(double) ) {             // If data type is double
      type = MPI_DOUBLE;                                        //  use MPI_DOUBLE
    } else if( typeid(object) == typeid(long double) ) {        // If data type is long double
      type = MPI_LONG_DOUBLE;                                   //  use MPI_LONG_DOUBLE
    } else if( typeid(object) == typeid(std::complex<float>) ) {// If data type is complex<float>
      type = MPI_COMPLEX;                                       //  use MPI_COMPLEX
    } else if( typeid(object) == typeid(std::complex<double>) ) {// If data type is compelx<double>
      type = MPI_DOUBLE_COMPLEX;                                //  use MPI_DOUBLE_COMPLEX
    }                                                           // Endif for data type
    return type;                                                // Return MPI data type
  }

//! Print a scalar value on all ranks
  template<typename T>
  void print(T data) {
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
      usleep(WAIT);                                             //  Wait "WAIT" milliseconds
      if( MPIRANK == irank ) std::cout << data << " ";          //  If it's my turn print "data"
    }                                                           // End loop over ranks
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( MPIRANK == 0 ) std::cout << std::endl;                  // New line
  }

//! Print a scalar value on irank
  template<typename T>
  void print(T data, const int irank) {
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( MPIRANK == irank ) std::cout << data;                   // If it's my rank print "data"
  }

//! Print a vector value on all ranks
  template<typename T>
  void print(T *data, const int begin, const int end) {
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
      usleep(WAIT);                                             //  Wait "WAIT" milliseconds
      if( MPIRANK == irank ) {                                  //  If it's my turn to print
        std::cout << MPIRANK << " : ";                          //   Print rank
        for( int i=begin; i!=end; ++i ) {                       //   Loop over data
          std::cout << data[i] << " ";                          //    Print data[i]
        }                                                       //   End loop over data
        std::cout << std::endl;                                 //   New line
      }                                                         //  Endif for my turn
    }                                                           // End loop over ranks
  }

//! Print a vector value on irank
  template<typename T>
  void print(T *data, const int begin, const int end, const int irank) {
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( MPIRANK == irank ) {                                    // If it's my rank
      std::cout << MPIRANK << " : ";                            //  Print rank
      for( int i=begin; i!=end; ++i ) {                         //  Loop over data
        std::cout << data[i] << " ";                            //   Print data[i]
      }                                                         //  End loop over data
      std::cout << std::endl;                                   //  New line
    }                                                           // Endif for my rank
  }
};

#endif
