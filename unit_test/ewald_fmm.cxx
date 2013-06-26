#include "serialfmm.h"

int main() {
  const int numBodies = 1000;                                   // Number of bodies
  const real_t xmax = 100.0;                                    // Size of domain
  const real_t ksize = 44.0;                                    // Ewald wave number
  const real_t alpha = 0.2;                                     // Ewald alpha value
  const real_t sigma = .25 / M_PI;                              // Ewald sigma value
  IMAGES = 8;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrt(4);                                          // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of bodies
  Bodies jbodies;                                               // Define vector of source bodies
  Cells cells, jcells;                                          // Define vector of cells
  SerialFMM<Laplace> FMM;                                       // Instantiate SerialFMM class
  FMM.initialize();                                             // Initialize FMM
  FMM.printNow = true;                                          // Print timings

  FMM.startTimer("Set bodies");                                 // Start timer
  srand48(0);                                                   // Seed for random number generator
  real_t average = 0;                                           // Initialize average charge
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over dimensions
      B->X[d] = drand48() * xmax;                               //   Initialize positions
    }                                                           //  End loop over dimensions
    B->SRC = drand48() / numBodies;                             //  Set charges
    average += B->SRC;                                          //  Accumulate charges
    B->TRG = 0;                                                 //  Initialize target values
    B->IBODY = B-bodies.begin();
  }                                                             // End loop over bodies
  average /= numBodies;                                         // Divide by total to get average
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->SRC -= average;                                          //  Make average charge 0
  }                                                             // End loop over bodies
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer
  FMM.eraseTimer("Set bodies");                                 // Erase entry from timer to avoid timer overlap

  FMM.startTimer("Set domain");                                 // Start timer
  FMM.setDomain(bodies,xmax/2,xmax/2);                          // Set domain size of FMM
  FMM.stopTimer("Set domain",FMM.printNow);                     // Stop timer
  FMM.eraseTimer("Set domain");                                 // Erase entry from timer to avoid timer overlap

  FMM.bottomup(bodies,cells);                                   // Tree construction (bottom up) & upward sweep
  FMM.setEwald(ksize,alpha,sigma);                              // Set Ewald method paramters
  jcells = cells;                                               // Copy cells to jcells
  FMM.Ewald(bodies,cells,jcells);                               // Ewald summation

  Bodies bodies2 = bodies;                                      // Define new bodies vector for direct sum
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->TRG = 0;
  FMM.startTimer("Downward");                                   // Start timer
  FMM.downward(cells,jcells);                                   // Downward sweep
  FMM.stopTimer("Downward",FMM.printNow);                       // Stop timer
  FMM.eraseTimer("Downward");                                   // Erase entry from timer to avoid timer overlap

  real_t diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;            // Initialize accumulators
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2,true);   // Evaluate error
  FMM.printError(diff1,norm1,diff2,norm2);                      // Print the L2 norm error
  FMM.finalize();                                               // Finalize FMM
}
