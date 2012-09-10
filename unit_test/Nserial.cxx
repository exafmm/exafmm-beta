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
  int numBodies = 10000;                                        // Number of bodies
  int numTarget = 100;                                          // Number of target points to be used for error eval
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies bodies, jbodies;                                       // Define vector of bodies
  Cells cells, jcells;                                          // Define vector of cells
  SerialFMM<Laplace> FMM;                                       // Instantiate SerialFMM class
  FMM.initialize();                                             // Initialize FMM
  bool printNow = true;                                         // Print timer

  for( int it=0; it!=25; ++it ) {                               // Loop over FMM iterations
    numBodies = int(pow(10,(it+24)/8.0));                       //  Exponentially increase N
    std::cout << "N             : " << numBodies << std::endl;  //  Print N
    bodies.resize(numBodies);                                   //  Resize bodies vector
    FMM.cube(bodies,1,1);                                       //  Initialize bodies in a cube
    FMM.startTimer("FMM");                                      //  Start timer
    FMM.setDomain(bodies);                                      //  Set domain size of FMM
    cells.clear();                                              //  Make sure cells vector is empty
#ifdef TOPDOWN
    FMM.topdown(bodies,cells);                                  //  Tree construction (top down) & upward sweep
#else
    FMM.bottomup(bodies,cells);                                 //  Tree construction (bottom up) & upward sweep
#endif
    jcells = cells;                                             //  Vector of source cells
    FMM.downward(cells,jcells);                                 //  Downward sweep
    FMM.stopTimer("FMM",printNow);                              //  Stop timer
    FMM.eraseTimer("FMM");                                      //  Erase entry from timer to avoid timer overlap

    FMM.startTimer("Direct sum");                               //  Start timer
    FMM.buffer = bodies;                                        //  Define new bodies vector for direct sum
#if 1
    FMM.initTarget(FMM.buffer);                                 //  Reinitialize target values
    if( IMAGES != 0 ) {                                         //  For periodic boundary condition
      jbodies = FMM.periodicBodies(FMM.buffer);                 //   Copy source bodies for all periodic images
    } else {                                                    //  For free field boundary condition
      jbodies = FMM.buffer;                                     //   Source bodies = target bodies
    }                                                           //  End if for periodic boundary condition
    FMM.buffer.resize(numTarget);                               //  Shrink target bodies vector to save time
    FMM.evalP2P(FMM.buffer,jbodies);                            //  Direct summation between buffer and jbodies
    FMM.writeTarget(FMM.buffer);                                //  Write direct summation results to file
#else
    FMM.readTarget(FMM.buffer);                                 //  Read direct summation results from file
#endif
    FMM.stopTimer("Direct sum",printNow);                       //  Stop timer
    FMM.eraseTimer("Direct sum");                               //  Erase entry from timer to avoid timer overlap
    FMM.writeTime();                                            //  Write timings of all events to file
    FMM.resetTimer();                                           //  Erase all events in timer

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;            //  Initialize accumulators
    bodies.resize(numTarget);                                   //  Shrink target bodies vector to save time
    FMM.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2);   //  Evaluate error on the reduced set of bodies
    FMM.printError(diff1,norm1,diff2,norm2);                    //  Print the L2 norm error
  }                                                             // End loop over FMM iteration
  FMM.finalize();                                               // Finalize FMM
}
