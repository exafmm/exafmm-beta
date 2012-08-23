#include "parallelfmm.h"

extern "C" void FMMcalccoulomb(int n, double* x, double* q, double *p, double* f, int periodicflag) {
  Bodies bodies, jbodies;
  Cells cells, jcells;
  ParallelFMM FMM;
  FMM.NCRIT = 10;
  FMM.NSPAWN = 1000;
  FMM.IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  FMM.THETA = 0.6;
  FMM.printNow = FMM.MPIRANK == 0;

  bodies.resize(n);
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
    B->IBODY  = i;
  }

  FMM.partition(bodies);
  FMM.buildTree(bodies, cells);
  FMM.upwardPass(cells);
  FMM.setLET(cells);
  FMM.commBodies();
  FMM.commCells();
  FMM.evaluate(cells, cells);
  jbodies = bodies;
  for( int irank=1; irank<FMM.MPISIZE; irank++ ) {
    FMM.getLET(jcells, (FMM.MPIRANK + irank) % FMM.MPISIZE);
    FMM.evaluate(cells, jcells);
  }
  FMM.downwardPass(cells);
  FMM.unpartition(bodies);

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
