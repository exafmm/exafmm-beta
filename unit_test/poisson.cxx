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
  const int k = 4;
  const int level = 4;
  const int numGrid = k * (1 << level);
  const int numBodies = numGrid * numGrid * numGrid;
  const real L = 10;
  std::cout << "N             : " << numBodies << std::endl;
  IMAGES = 0;
  THETA = 1 / sqrt(4);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells, jcells;
  SerialFMM<Laplace> FMM;
  FMM.initialize();
  FMM.printNow = true;
  assert( NCRIT == k*k*k );

  FMM.startTimer("Set bodies   ");
  real dV = 8. / numBodies;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    int ix = i / numGrid / numGrid;
    int iy = (i % (numGrid*numGrid)) / numGrid;
    int iz = i % numGrid;
    B->X[0] = ix * (2. / numGrid) - 1 + 1. / numGrid;
    B->X[1] = iy * (2. / numGrid) - 1 + 1. / numGrid;
    B->X[2] = iz * (2. / numGrid) - 1 + 1. / numGrid;
    B->SRC = -(4 * L * L * norm(B->X) - 6 * L) * exp(-L * norm(B->X)) * dV * .25 / M_PI;
    B->TRG = 0;
  }
  vect X0 = 0;
  FMM.setX0(X0);
  FMM.setR0(1);
  FMM.stopTimer("Set bodies   ",FMM.printNow);
  FMM.eraseTimer("Set bodies   ");

  FMM.bottomup(bodies,cells);
  jcells = cells;
  FMM.startTimer("Downward     ");
  FMM.downward(cells,jcells);
  FMM.stopTimer("Downward     ",FMM.printNow);
  FMM.eraseTimer("Downward     ");

  real diff1 = 0, norm1 = 0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    real exact = exp(-L * norm(B->X));
    diff1 += (B->TRG[0] - exact) * (B->TRG[0] - exact);
    norm1 += exact * exact;
  }
  std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;

  FMM.finalize();
}
