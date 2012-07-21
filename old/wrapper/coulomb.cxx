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

extern "C" void FMMcalccoulomb(int n, double* x, double* q, double *p, double* f, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  THETA = .5;
  Bodies bodies(n),jbodies;
  Cells cells,jcells;
  ParallelFMM<Laplace> FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    B->SRC  = q[i];
    B->TRG[0] =  p[i];
    B->TRG[1] = -f[3*i+0];
    B->TRG[2] = -f[3*i+1];
    B->TRG[3] = -f[3*i+2];
    B->IBODY = i;
    B->IPROC = MPIRANK;
  }

  FMM.initialize();
  FMM.setGlobDomain(bodies);
  FMM.octsection(bodies);
  FMM.bottomup(bodies,cells);
  FMM.commBodies(cells);
  jbodies = bodies;
  jcells = cells;
  FMM.commCells(jbodies,jcells);

  FMM.downward(cells,jcells);
  FMM.unpartition(bodies);
  std::sort(bodies.begin(),bodies.end());
  FMM.finalize();

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i]     = B->SRC;
    p[i]     =  B->TRG[0];
    f[3*i+0] = -B->TRG[1];
    f[3*i+1] = -B->TRG[2];
    f[3*i+2] = -B->TRG[3];
  }
}
