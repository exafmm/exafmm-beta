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
#include "mympi.h"
#include "evaluator.h"

int main() {
  char hostname[256];                                           // Define hostname
  const int numBodies = 1000;                                   // Number of bodies
  IMAGES = 1;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of bodies
  Evaluator<Laplace> FMM;                                       // Instantiate Evaluator class
  MyMPI M;                                                      // Instantiate MyMPI class
  FMM.initialize();                                             // Initialize FMM
  FMM.preCalculation();                                         // Kernel pre-processing
  gethostname(hostname,sizeof(hostname));                       // Get hostname
  if( MPIRANK == 0 ) FMM.printNow = true;                       // Print only if MPIRANK == 0

  FMM.startTimer("Set bodies");                                 // Start timer
  FMM.sphere(bodies);                                           // Initialize bodies on a spherical shell
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer

  FMM.startTimer("Set domain");                                 // Start timer
  FMM.setDomain(bodies);                                        // Set domain size of FMM
  FMM.stopTimer("Set domain",FMM.printNow);                     // Stop timer

  FMM.startTimer("Direct GPU");                                 // Start timer
  FMM.evalP2P(bodies,bodies);                                   // Direct summation on GPU
  FMM.stopTimer("Direct GPU",FMM.printNow);                     // Stop timer

  FMM.startTimer("Direct CPU");                                 // Start timer
  bool onCPU = true;                                            // Bool for CPU run
  Bodies bodies2 = bodies;                                      // Define new bodies vector for direct sum
  FMM.initTarget(bodies2);                                      // Reinitialize target values
  FMM.evalP2P(bodies2,bodies2,onCPU);                           // Direct summation on CPU
  FMM.stopTimer("Direct CPU",FMM.printNow);                     // Stop timer

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;              // Initialize accumulators
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);        // Evaluate error

  for( int irank=0; irank!=MPISIZE; ++irank ) {                 // Loop over MPI ranks
    MPI_Barrier(MPI_COMM_WORLD);                                //  Synchronize processes
    if( MPIRANK == irank ) {                                    //  If loop counter matches MPI rank
      std::cout << hostname << " @ rank : " << MPIRANK << " / " << MPISIZE;// Print hostname MPI rank and size
      std::cout << " @ device : " << DEVICE << " / " << GPUS << std::endl;// Print GPU device ID and device count
      std::cout << std::setw(FMM.stringLength) << std::left     //   Set format
                << "Error" << " : " << std::sqrt(diff1/norm1) << std::endl;// Print relative L2 norm error
    }                                                           //  Endif for loop counter
    usleep(100);                                                //  Wait for 100 microseconds
  }                                                             // End loop over MPI ranks
  FMM.postCalculation();                                        // Kernel post-processing
  FMM.finalize();                                               // Finalize FMM
}
