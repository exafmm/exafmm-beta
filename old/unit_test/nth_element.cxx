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
#include "partition.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 1000000;                                // Number of bodies
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of bodies
  Partition<Laplace> FMM;                                       // Instantiate Partition class
  FMM.initialize();                                             // Initialize FMM
  if( MPIRANK == 0 ) FMM.printNow = true;                       // Print only if MPIRANK == 0

  FMM.startTimer("Set bodies");                                 // Start timer
  FMM.cube(bodies,MPIRANK+1);                                   // Initialize bodies in a cube
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer

  FMM.startTimer("Set domain");                                 // Start timer
  FMM.setGlobDomain(bodies);                                    // Set global domain size of FMM
  FMM.stopTimer("Set domain",FMM.printNow);                     // Stop timer

  FMM.setIndex(bodies);                                         // Set index of cells
  FMM.binBodies(bodies,0);                                      // Bin bodies into leaf level cells

  FMM.buffer.resize(bodies.size());                             // Resize sort buffer
  FMM.sortBodies(bodies,FMM.buffer);                            // Sort bodies in ascending order

  FMM.startTimer("Nth element");                                // Start timer
  bigint nthGlobal = numBodies * MPISIZE / 3;                   // Split at nth global element
  bigint iSplit = FMM.nth_element(bodies,nthGlobal);            // Get cell index of nth global element
  int nthLocal = FMM.splitBodies(bodies,iSplit);                // Split bodies based on iSplit
  FMM.stopTimer("Nth element",FMM.printNow);                    // Stop timer
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->ICELL = B-bodies.begin() > nthLocal;                     //  Set cell index according to the split
  }                                                             // End loop over bodies

#ifdef VTK
  if( MPIRANK == 0 ) {                                          // If MPI rank is 0
    int Ncell = 0;                                              //  Initialize number of cells
    vtkPlot vtk;                                                //  Instantiate vtkPlot class
    vtk.setDomain(FMM.getR0(),FMM.getX0());                     //  Set bounding box for VTK
    vtk.setGroupOfPoints(bodies,Ncell);                         //  Set group of points
    vtk.plot(Ncell);                                            //  plot using VTK
  }                                                             // Endif for MPI rank
#endif
  FMM.finalize();                                               // Finalize FMM
}
