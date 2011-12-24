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
#include "tree.h"

int main() {
  int numBodies = 1000;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies, bodies2;
  Cells cells;
  Dataset DATA;
  SerialFMM FMM;
#if HYBRID
  FMM.timeKernels();
#endif
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
  FMM.printNow = true;
#if BUILD
  for ( int it=32; it<33; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
#endif
  numBodies = int(pow(10,(it+24)/8.0));
  std::cout << "N             : " << numBodies << std::endl;
  bodies.resize(numBodies);
  DATA.random(bodies);
  FMM.startTimer("FMM          ");
#if BOTTOMUP
  FMM.bottomup(bodies,cells);
#else
  FMM.topdown(bodies,cells);
#endif
#if BUILD
#else
#if IneJ
  FMM.evaluate(cells,cells);
#else
  FMM.evaluate(cells);
#endif
  FMM.stopTimer("FMM          ",true);
  FMM.eraseTimer("FMM          ");
  FMM.writeTime();
  FMM.resetTimer();

#ifdef DIRECT
  bodies2 = bodies;
  DATA.initTarget(bodies);
  FMM.startTimer("Direct sum   ");
#if IneJ
  FMM.direct(bodies,bodies);
#else
  FMM.direct(bodies);
#endif
  FMM.stopTimer("Direct sum   ",true);
  FMM.eraseTimer("Direct sum   ");

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  DATA.printError(diff1,norm1,diff2,norm2);
#endif
#endif
  }
}
