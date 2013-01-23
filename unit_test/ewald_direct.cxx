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
#include "serialfmm.h"

int main() {
  const int numBodies = 1000;                                   // Number of bodies
  const int numTarget = 100;                                    // Number of target points to be used for error eval
  const real xmax = 100.0;                                      // Size of domain
  const real ksize = 11.0;                                      // Ewald wave number
  const real alpha = 0.1;                                       // Ewald alpha value
  const real sigma = .25 / M_PI;                                // Ewald sigma value
  IMAGES = 3;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrt(4);                                          // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of bodies
  Bodies jbodies;                                               // Define vector of source bodies
  Cells cells, jcells;                                          // Define vector of cells
  SerialFMM<Laplace> FMM;                                       // Instantiate SerialFMM class
  FMM.initialize();                                             // Initialize FMM
  FMM.printNow = true;                                          // Print timings

  FMM.startTimer("Set bodies");                                 // Start timer
  srand48(2);                                                   // Seed for random number generator
  real average = 0;                                             // Initialize average charge
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over dimensions
      B->X[d] = drand48() * xmax;                               //   Initialize positions
    }                                                           //  End loop over dimensions
    B->SRC = drand48();                                         //  Set charges
    average += B->SRC;                                          //  Accumulate charges
    B->TRG = 0;                                                 //  Initialize target values
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

  FMM.startTimer("Direct sum");                                 // Start timer
  jbodies = bodies;                                             // Save jbodies
  FMM.sampleBodies(bodies,numTarget);                           // Shrink target bodies vector to save time
  FMM.buffer = bodies;                                          // Define new bodies vector for direct sum
  FMM.initTarget(FMM.buffer);                                   // Reinitialize target values
  FMM.evalP2P(FMM.buffer,jbodies);                              // Direct summation between buffer and jbodies
  FMM.stopTimer("Direct sum",FMM.printNow);                     // Stop timer
  FMM.eraseTimer("Direct sum");                                 // Erase entry from timer to avoid timer overlap

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;              // Initialize accumulators
  FMM.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2,true);// Evaluate error
  FMM.printError(diff1,norm1,diff2,norm2);                      // Print the L2 norm error
  FMM.finalize();                                               // Finalize FMM
}
