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
  const int numBodies = 1000000;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies bodies(numBodies);
  Partition<Laplace> FMM;
  FMM.initialize();
  if( MPIRANK == 0 ) FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  FMM.random(bodies,MPIRANK+1);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setGlobDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  FMM.setIndex(bodies);
  FMM.binBodies(bodies,0);

  FMM.buffer.resize(bodies.size());
  FMM.sortBodies(bodies,FMM.buffer);

  FMM.startTimer("Nth element  ");
  bigint nthGlobal = numBodies * MPISIZE / 3;
  bigint iSplit = FMM.nth_element(bodies,nthGlobal);
  int nthLocal = FMM.splitBodies(bodies,iSplit);
  FMM.stopTimer("Nth element  ",FMM.printNow);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    B->ICELL = B-bodies.begin() > nthLocal;
  }

#ifdef VTK
  if( MPIRANK == 0 ) {
    int Ncell = 0;
    vtkPlot vtk;
    vtk.setDomain(FMM.getR0(),FMM.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
  FMM.finalize();
}
