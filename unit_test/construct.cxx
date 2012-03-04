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
  const int numBodies = 1000000;                                // Number of bodies
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of bodies
  Cells cells;                                                  // Define vector of cells
  SerialFMM<Laplace> FMM;                                       // Instantiate SerialFMM class
  FMM.initialize();                                             // Initialize FMM
  FMM.printNow = true;                                          // Print timings

  FMM.startTimer("Set bodies");                                 // Start timer
  FMM.cube(bodies);                                             // Initialize bodies in a cube
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer

  FMM.startTimer("Set domain");                                 // Start timer
  FMM.setDomain(bodies);                                        // Set domain size of FMM
  FMM.stopTimer("Set domain",FMM.printNow);                     // Stop timer

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);                                    // Tree construction (top down) & upward sweep
#else
  FMM.bottomup(bodies,cells);                                   // Tree construction (bottom up) & upward sweep
#endif

#ifdef VTK
  int Ncell = 0;                                                // Initialize number of cells
  vtkPlot vtk;                                                  // Instantiate vtkPlot class
  vtk.setDomain(FMM.getR0(),FMM.getX0());                       // Set bounding box for VTK
  vtk.setGroupOfPoints(bodies,Ncell);                           // Set group of points
  vtk.plot(Ncell);                                              // plot using VTK
#endif
  FMM.finalize();                                               // Finalize FMM
}
