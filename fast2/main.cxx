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
#include "dataset.h"
#include "parallelfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, jbodies;
  Cells cells;
  Dataset DATA;
  ParallelFMM<Laplace> FMM;
  bool printNow = MPIRANK == 0;
#if HYBRID
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
  FMM.printNow = MPIRANK == 0;
#if BUILD
  for ( int it=32; it<33; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
#endif
  numBodies = int(pow(10,(it+24)/8.0));
  if( printNow ) std::cout << "N                    : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.cube(bodies,MPIRANK);
  FMM.startTimer("FMM");
#if BOTTOMUP
  FMM.bottomup(bodies,cells);
#else
  FMM.topdown(bodies,cells);
#endif
#if BUILD
#else
  Cells jcells = cells;
  FMM.startPAPI();
  FMM.evaluate(cells,jcells);
  FMM.stopPAPI();
  FMM.stopTimer("FMM",printNow);
  FMM.eraseTimer("FMM");
  FMM.writeTime();
  FMM.resetTimer();
  if(FMM.printNow) FMM.printTreeData(cells);

  Bodies bodies2 = bodies;
#ifdef MANY
  bodies2.resize(100);
#endif
  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic");
    jbodies = FMM.periodicBodies(bodies);
    FMM.stopTimer("Set periodic",printNow);
    FMM.eraseTimer("Set periodic");
  } else {
    jbodies = bodies;
  }
  DATA.initTarget(bodies2);
  FMM.startTimer("Direct sum");
  FMM.direct(bodies2,jbodies);
  FMM.stopTimer("Direct sum",printNow);
  FMM.eraseTimer("Direct sum");

#ifdef MANY
  bodies.resize(100);
#endif
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  if( printNow ) DATA.printError(diff1,norm1,diff2,norm2);
#endif
  }

//  FMM.octsection(bodies);
//  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->ICELL = 0;
//  FMM.commBodies(cells);
#ifdef VTK
  int Ncell = 0;                                                // Initialize number of cells
  vtkPlot vtk;                                                  // Instantiate vtkPlot class
  if( MPIRANK == 0 ) {                                          // If MPI rank is 0
    vtk.setDomain(FMM.getR0(),FMM.getX0());                     //  Set bounding box for VTK
    vtk.setGroupOfPoints(bodies,Ncell);                         //  Set group of points
  }                                                             // Endif for MPI rank
  for( int i=1; i!=MPISIZE; ++i ) {                             // Loop over MPI ranks
    FMM.shiftBodies(bodies);                                    //  Communicate bodies round-robin
    if( MPIRANK == 0 ) {                                        //  If MPI rank is 0
      vtk.setGroupOfPoints(bodies,Ncell);                       //   Set group of points
    }                                                           //  Endif for MPI rank
  }                                                             // End loop over MPI ranks
  if( MPIRANK == 0 ) {                                          // If MPI rank is 0
    vtk.plot(Ncell);                                            //  plot using VTK
  }                                                             // Endif for MPI rank
#endif
}
