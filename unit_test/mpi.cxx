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
#include <mpi.h>
#include <iostream>
#include <unistd.h>

int main(int argc, char **argv) {
  char hostname[256];                                           // Define hostname
  int size,rank;                                                // Define MPI size and rank
  MPI_Init(&argc,&argv);                                        // Initialize MPI communicator
  MPI_Comm_size(MPI_COMM_WORLD,&size);                          // Get number of MPI processes
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);                          // Get rank of current MPI process
  gethostname(hostname,sizeof(hostname));                       // Get hostname
  for( int irank=0; irank!=size; ++irank ) {                    // Loop over MPI ranks
    MPI_Barrier(MPI_COMM_WORLD);                                //  Synchronize processes
    if( rank == irank ) {                                       //  If loop counter matches MPI rank
      std::cout << hostname << " " << rank << " / " << size << std::endl;// Print hostname, rank, and size
    }                                                           //  Endif for loop counter
    usleep(100);                                                //  Wait for 100 microseconds
  }                                                             // End loop over MPI ranks
  MPI_Finalize();                                               // Finalize MPI communicator
}
