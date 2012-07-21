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
#include "parallelfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;                                  // Number of bodies
  const int numTarget = 100;                                    // Number of target points to be used for error eval
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies bodies(numBodies);                                     // Define vector of target bodies
  Bodies jbodies;                                               // Define vector of source bodies
  Cells cells, jcells;                                          // Define vector of cells
  ParallelFMM<Laplace> FMM;                                     // Instantiate ParallelFMM class
  FMM.initialize();                                             // Initialize FMM
  if( MPIRANK == 0 ) FMM.printNow = true;                       // Print only if MPIRANK == 0

  FMM.startTimer("Set bodies");                                 // Start timer
  FMM.cube(bodies,MPIRANK+1);                                   // Initialize bodies in a cube
  Bodies bodies2 = bodies;                                      // Define another vector of target bodies
  FMM.stopTimer("Set bodies",FMM.printNow);                     // Stop timer

  FMM.startTimer("Set domain");                                 // Start timer
  FMM.setGlobDomain(bodies2);                                   // Set global domain size of FMM
  FMM.stopTimer("Set domain",FMM.printNow);                     // Stop timer

#ifndef VTK
 if( IMAGES != 0 ) {                                            // For periodic boundary condition
    FMM.startTimer("Set periodic");                             //  Start timer
    jbodies = FMM.periodicBodies(bodies2);                      //  Copy source bodies for all periodic images
    FMM.stopTimer("Set periodic",FMM.printNow);                 //  Stop timer
    FMM.eraseTimer("Set periodic");                             //  Erase entry from timer to avoid timer overlap
  } else {                                                      // For free field boundary condition
    jbodies = bodies2;                                          //  Copy source bodies
  }                                                             // End if for periodic boundary condition
  FMM.startTimer("Direct sum");                                 // Start timer
  bodies2.resize(numTarget);                                    // Shrink target bodies vector to save time
  FMM.initTarget(bodies2);                                      // Reinitialize target values
  for( int i=0; i!=MPISIZE; ++i ) {                             // Loop over all MPI processes
    FMM.shiftBodies(jbodies);                                   //  Communicate bodies round-robin
    FMM.evalP2P(bodies2,jbodies);                               //  Direct summation between bodies2 and jbodies2
    if(FMM.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;// Print loop counter
  }                                                             // End loop over all MPI processes
  FMM.stopTimer("Direct sum",FMM.printNow);                     // Stop timer
  FMM.eraseTimer("Direct sum");                                 // Erase entry from timer to avoid timer overlap
#endif

  FMM.initTarget(bodies);                                       // Reinitialize target values

  FMM.octsection(bodies);                                       // Partition domain and redistribute bodies

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);                                    // Tree construction (top down) & upward sweep
#else
  FMM.bottomup(bodies,cells);                                   // Tree construction (bottom up) & upward sweep
#endif

  FMM.commBodies(cells);                                        // Send bodies (not receiving yet)

  jbodies = bodies;                                             // Vector of source bodies
  jcells = cells;                                               // Vector of source cells
  FMM.commCells(jbodies,jcells);                                // Communicate cells (receive bodies here)

  FMM.startTimer("Downward");                                   // Start timer
  FMM.downward(cells,jcells);                                   // Downward sweep
  FMM.stopTimer("Downward",FMM.printNow);                       // Stop timer
  FMM.eraseTimer("Downward");                                   // Erase entry from timer to avoid timer overlap

  FMM.startTimer("Unpartition");                                // Start timer
  FMM.unpartition(bodies);                                      // Send bodies back to where they came from
  FMM.stopTimer("Unpartition",FMM.printNow);                    // Stop timer
  FMM.eraseTimer("Unpartition");                                // Erase entry from timer to avoid timer overlap

  FMM.startTimer("Unsort bodies");                              // Start timer
  std::sort(bodies.begin(),bodies.end());                       // Sort bodies back to original order
  FMM.stopTimer("Unsort bodies",FMM.printNow);                  // Stop timer
  FMM.eraseTimer("Unsort bodies");                              // Erase entry from timer to avoid timer overlap
  if(FMM.printNow) FMM.writeTime();                             // Write timings of all events to file
  if(FMM.printNow) FMM.writeTime();                             // Write again to have at least two data sets

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);        // Evaluate error on the reduced set of bodies
  MPI_Datatype MPI_TYPE = FMM.getType(diff1);                   // Get MPI datatype
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Reduce difference in potential
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Reduce norm of potential
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Reduce difference in force
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);// Recude norm of force
  if(FMM.printNow) FMM.printError(diff3,norm3,diff4,norm4);     // Print the L2 norm error

#ifdef VTK
  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) B->ICELL = 0;// Reinitialize cell index
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {        // Loop over source cells
    Body body;                                                  //  Create one body per jcell
    body.ICELL = 1;                                             //  Set cell index to 1
    body.X     = C->X;                                          //  Copy cell position to body
    body.SRC   = 0;                                             //  Set source value to 0
    jbodies.push_back(body);                                    //  Push body into vector
  }                                                             // End loop over source cells

  int Ncell = 0;                                                // Initialize number of cells
  vtkPlot vtk;                                                  // Instantiate vtkPlot class
  if( MPIRANK == 0 ) {                                          // If MPI rank is 0
    vtk.setDomain(FMM.getR0(),FMM.getX0());                     //  Set bounding box for VTK
    vtk.setGroupOfPoints(jbodies,Ncell);                        //  Set group of points
  }                                                             // Endif for MPI rank
  for( int i=1; i!=MPISIZE; ++i ) {                             // Loop over MPI ranks
    FMM.shiftBodies(jbodies);                                   //  Communicate bodies round-robin
    if( MPIRANK == 0 ) {                                        //  If MPI rank is 0
      vtk.setGroupOfPoints(jbodies,Ncell);                      //   Set group of points
    }                                                           //  Endif for MPI rank
  }                                                             // End loop over MPI ranks
  if( MPIRANK == 0 ) {                                          // If MPI rank is 0
    vtk.plot(Ncell);                                            //  plot using VTK
  }                                                             // Endif for MPI rank
#endif
  FMM.finalize();                                               // Finalize FMM
}
