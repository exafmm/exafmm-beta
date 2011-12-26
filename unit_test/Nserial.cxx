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
  int numBodies = 10000;
  int numTarget = 100;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  SerialFMM<Laplace> FMM;
  FMM.initialize();

  for( int it=0; it!=25; ++it ) {
    numBodies = int(pow(10,(it+24)/8.0));
    std::cout << "N             : " << numBodies << std::endl;
    bodies.resize(numBodies);
    FMM.random(bodies,1,1);
    FMM.startTimer("FMM          ");
    FMM.setDomain(bodies);
    cells.clear();
#ifdef TOPDOWN
    FMM.topdown(bodies,cells);
#else
    FMM.bottomup(bodies,cells);
#endif
    jcells = cells;
    FMM.downward(cells,jcells);
    FMM.stopTimer("FMM          ",true);
    FMM.eraseTimer("FMM          ");

    FMM.startTimer("Direct sum   ");
    FMM.buffer = bodies;
#if 1
    FMM.initTarget(FMM.buffer);
    if( IMAGES != 0 ) {
      jbodies = FMM.periodicBodies(FMM.buffer);
    } else {
      jbodies = FMM.buffer;
    }
    FMM.buffer.resize(numTarget);
    FMM.evalP2P(FMM.buffer,jbodies);
    FMM.writeTarget(FMM.buffer);
#else
    FMM.readTarget(FMM.buffer);
#endif
    FMM.stopTimer("Direct sum   ",true);
    FMM.eraseTimer("Direct sum   ");
    FMM.writeTime();
    FMM.resetTimer();

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    bodies.resize(numTarget);
    FMM.evalError(bodies,FMM.buffer,diff1,norm1,diff2,norm2);
    FMM.printError(diff1,norm1,diff2,norm2);
  }
  FMM.finalize();
}
