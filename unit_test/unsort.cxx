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
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies bodies(numBodies);
  Bodies bodies2;
  Bodies jbodies;
  Cells cells, jcells;
  SerialFMM<Laplace> FMM;
  FMM.initialize();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  FMM.random(bodies,1,1);
  bodies2 = bodies;
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies2);
    FMM.stopTimer("Set periodic ",FMM.printNow);
  } else {
    jbodies = bodies2;
  }

  FMM.startTimer("Direct sum   ");
  bodies2.resize(numTarget);
  FMM.evalP2P(bodies2,jbodies);
  FMM.stopTimer("Direct sum   ",FMM.printNow);
  FMM.eraseTimer("Direct sum   ");

  FMM.startTimer("Set domain   ");
  FMM.initTarget(bodies);
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);
#else
  FMM.bottomup(bodies,cells);
#endif

  jcells = cells;
  FMM.startTimer("Downward     ");
  FMM.downward(cells,jcells);
  FMM.stopTimer("Downward     ",FMM.printNow);

  FMM.startTimer("Unsort bodies");
  std::sort(bodies.begin(),bodies.end());
  FMM.stopTimer("Unsort bodies",FMM.printNow);
  FMM.writeTime();
  FMM.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  bodies.resize(numTarget);
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  FMM.printError(diff1,norm1,diff2,norm2);
  FMM.finalize();
}
