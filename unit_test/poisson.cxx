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
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int k = 4;                                              // Order of polynomial
  const int level = 4;                                          // Depth of tree
  const int numGrid = k * (1 << level);                         // Number of grid points in one direction
  const int numBodies = numGrid * numGrid * numGrid;            // Number of bodies
  const real L = 10;                                            // Size of Gaussian distribution
  std::cout << "N             : " << numBodies << std::endl;    // Print number of bodies
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrt(4);                                          // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of target bodies
  Bodies jbodies;                                               // Define vector of source bodies
  Cells cells, jcells;                                          // Define vector of cells
  SerialFMM<Laplace> FMM;                                       // Instantiate SerialFMM class
  FMM.initialize();                                             // Initialize FMM
  FMM.printNow = true;                                          // Print timings
  assert( NCRIT == k*k*k );                                     // NCRIT must be set to k^3 for this app.

  FMM.startTimer("Set bodies");                                 // Start timer
  real dV = 8. / numBodies;                                     // Volume per body
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    int i = B-bodies.begin();                                   //  Loop counter
    int ix = i / numGrid / numGrid;                             //  x index
    int iy = (i % (numGrid*numGrid)) / numGrid;                 //  y index
    int iz = i % numGrid;                                       //  z index
    B->X[0] = ix * (2. / numGrid) - 1 + 1. / numGrid;           //  x position
    B->X[1] = iy * (2. / numGrid) - 1 + 1. / numGrid;           //  y position
    B->X[2] = iz * (2. / numGrid) - 1 + 1. / numGrid;           //  z position
    B->SRC = -(4 * L * L * norm(B->X) - 6 * L) * exp(-L * norm(B->X)) * dV * .25 / M_PI;// Gaussian source
    B->TRG = 0;                                                 //  Initialize target values
  }                                                             // End loop over bodies
  vec3 X0 = 0;                                                  // Center of root cell
  FMM.setX0(X0);                                                // Set center of root cell
  FMM.setR0(1);                                                 // Set radius of root cell
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer
  FMM.eraseTimer("Set bodies");                                 // Erase entry from timer to avoid timer overlap

  FMM.bottomup(bodies,cells);                                   // Tree construction (bottom up) & upward sweep
  jcells = cells;                                               // Vector of source cells
  FMM.startTimer("Downward");                                   // Start timer
  FMM.downward(cells,jcells);                                   // Downward sweep
  FMM.stopTimer("Downward",FMM.printNow);                       // Stop timer
  FMM.eraseTimer("Downward");                                   // Erase entry from timer to avoid timer overlap

  real diff1 = 0, norm1 = 0;                                    // Initialize accumulators
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    real exact = exp(-L * norm(B->X));                          //  Analytical solution
    diff1 += (B->TRG[0] - exact) * (B->TRG[0] - exact);         //  Difference between analytical solution
    norm1 += exact * exact;                                     //  L2 norm of analytical solution
  }                                                             // End loop over bodies
  std::cout << std::setw(stringLength) << std::left             // Set format
            << "Error" << " : " << std::sqrt(diff1/norm1) << std::endl;// Print L2 norm error

  FMM.finalize();                                               // Finalize FMM
}
