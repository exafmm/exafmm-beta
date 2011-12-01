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
#include "evaluator.h"

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Evaluator<Laplace> FMM;
  FMM.initialize();
  FMM.preCalculation();
  FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  FMM.sphere(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies);
    FMM.stopTimer("Set periodic ",FMM.printNow);
  } else {
    jbodies = bodies;
  }

  FMM.startTimer("Direct GPU   ");
  FMM.evalP2P(bodies,jbodies);
  FMM.stopTimer("Direct GPU   ",FMM.printNow);

  FMM.startTimer("Direct CPU   ");
  bool onCPU = true;
  bodies.resize(numTarget);
  Bodies bodies2 = bodies;
  FMM.initTarget(bodies2);
  FMM.evalP2P(bodies2,jbodies,onCPU);
  FMM.stopTimer("Direct CPU   ",FMM.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  FMM.printError(diff1,norm1,diff2,norm2);
  FMM.postCalculation();
  FMM.finalize();
}
