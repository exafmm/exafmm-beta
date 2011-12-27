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
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  const int numBodies = 10000;
  const int numTarget = 100;
  IMAGES = 0;
  THETA = 1 / sqrtf(4);
  Bodies bodies(numBodies);
  Bodies jbodies;
  Cells cells, jcells;
  ParallelFMM<Laplace> FMM;
  FMM.initialize();
  if( MPIRANK == 0 ) FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  FMM.random(bodies,MPIRANK+1);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setGlobDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  FMM.octsection(bodies);

#ifdef TOPDOWN
  FMM.topdown(bodies,cells);
#else
  FMM.bottomup(bodies,cells);
#endif

  FMM.commBodies(cells);

#ifndef VTK
  if( IMAGES != 0 ) {
    FMM.startTimer("Set periodic ");
    jbodies = FMM.periodicBodies(bodies);
    FMM.stopTimer("Set periodic ",FMM.printNow);
    FMM.eraseTimer("Set periodic ");
  } else {
    jbodies = bodies;
  }
  FMM.startTimer("Direct sum   ");
  Bodies bodies2 = bodies;
  bodies2.resize(numTarget);
  FMM.initTarget(bodies2);
  for( int i=0; i!=MPISIZE; ++i ) {
    FMM.shiftBodies(jbodies);
    FMM.evalP2P(bodies2,jbodies);
    if(FMM.printNow) std::cout << "Direct loop   : " << i+1 << "/" << MPISIZE << std::endl;
  }
  FMM.stopTimer("Direct sum   ",FMM.printNow);
#endif

  FMM.resetTimer();
  FMM.initTarget(bodies);
  FMM.evalP2M(cells);
  FMM.evalM2M(cells,cells);
  jcells = cells;
  if( MPISIZE != 1 ) {
    #pragma omp parallel sections num_threads(2)
    {
      #pragma omp section
      {
        FMM.downward(cells,jcells,false);
      }
      #pragma omp section
      {
        FMM.updateBodies();
      }
    }
    jbodies = bodies;
    jcells = cells;
    FMM.commCells(jbodies,jcells);
    FMM.eraseLocalTree(jcells);
  }
  FMM.downward(cells,jcells);
  if(FMM.printNow) FMM.writeTime();
  if(FMM.printNow) FMM.writeTime();

#ifndef VTK
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  bodies.resize(numTarget);
  FMM.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  MPI_Datatype MPI_TYPE = FMM.getType(diff1);
  MPI_Reduce(&diff1,&diff3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm1,&norm3,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&diff2,&diff4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&norm2,&norm4,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(FMM.printNow) FMM.printError(diff3,norm3,diff4,norm4);

#else
  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) B->ICELL = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }

  int Ncell = 0;
  vtkPlot vtk;
  if( MPIRANK == 0 ) {
    vtk.setDomain(FMM.getR0(),FMM.getX0());
    vtk.setGroupOfPoints(jbodies,Ncell);
  }
  for( int i=1; i!=MPISIZE; ++i ) {
    FMM.shiftBodies(jbodies);
    if( MPIRANK == 0 ) {
      vtk.setGroupOfPoints(jbodies,Ncell);
    }
  }
  if( MPIRANK == 0 ) {
    vtk.plot(Ncell);
  }
#endif
  FMM.finalize();
}
