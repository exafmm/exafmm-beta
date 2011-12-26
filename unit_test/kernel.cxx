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
  const int numBodies = 100;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies ibodies(numBodies);
  Bodies ibodies2(numBodies);
  Bodies jbodies(numBodies);
  Cells  icells;
  Cells  jcells;
  Evaluator<Laplace> FMM;
  FMM.initialize();
  FMM.preCalculation();

  for( int it=0; it!=10; ++it ) {
    real dist = (1 << it) / 2;
    for( B_iter B=ibodies.begin(); B!=ibodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->X[d] = -drand48() - dist;
      }
    }
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        B->X[d] = drand48();
      }
    }
    FMM.initSource(jbodies);
    bool IeqJ = false;
    FMM.initTarget(ibodies,IeqJ);

    Cell cell;
    cell.NDLEAF    = numBodies;
    cell.LEAF     = jbodies.begin();
    cell.X        = 0.5;
    cell.M        = 0;
    cell.ICELL    = 8;
    cell.NCHILD   = 0;
    cell.PARENT   = 1;
    jcells.push_back(cell);
    FMM.evalP2M(jcells);
    cell.X        = 1;
    cell.M        = 0;
    cell.ICELL    = 0;
    cell.NCHILD   = 1;
    cell.CHILD    = 0;
    jcells.push_back(cell);
    FMM.evalM2M(jcells,jcells);
    jcells.erase(jcells.begin());
    cell.X        = -1 - dist;
    cell.M        = 1;
    cell.L        = 0;
    cell.ICELL    = 0;
    icells.push_back(cell);
    FMM.addM2L(jcells.begin());
    FMM.evalM2L(icells);
    cell.NDLEAF    = numBodies;
    cell.LEAF     = ibodies.begin();
    cell.X        = -0.5 - dist;
    cell.L        = 0;
    cell.ICELL    = 1;
    cell.NCHILD   = 0;
    cell.PARENT   = 1;
    icells.insert(icells.begin(),cell);
    FMM.evalL2L(icells);
    icells.pop_back();
    FMM.evalL2P(icells);

    ibodies2 = ibodies;
    FMM.initTarget(ibodies2,IeqJ);
    FMM.evalP2P(ibodies2,jbodies);

    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    FMM.evalError(ibodies,ibodies2,diff1,norm1,diff2,norm2);
    std::cout << "Distance      : " << dist << std::endl;
    FMM.printError(diff1,norm1,diff2,norm2);

    FMM.initTarget(ibodies);
    FMM.addM2P(jcells.begin());
    FMM.evalM2P(icells);
    icells.clear();
    jcells.clear();
    diff1 = norm1 = diff2 = norm2 = 0;
    FMM.evalError(ibodies,ibodies2,diff1,norm1,diff2,norm2);
    FMM.printError(diff1,norm1,diff2,norm2);
  }
  FMM.postCalculation();
  FMM.finalize();
}
