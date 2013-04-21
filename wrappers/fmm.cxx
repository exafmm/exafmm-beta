#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"

extern "C" void FMM(int n, double* x, double* q, double *p, double* f, int periodicflag) {
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Sort sort;

  const int NCRIT = 16;
  const int NSPAWN = 1000;
  const int IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  const real_t THETA = 0.6;
  const real_t cycle = 2 * M_PI;

  BoundBox boundbox(NSPAWN);
  BuildTree tree(NCRIT,NSPAWN);
  UpDownPass pass(IMAGES,THETA);
  Traversal traversal(NSPAWN,IMAGES);
  LocalEssentialTree LET(IMAGES);
  boundbox.printNow = LET.MPIRANK == 0;
  tree.printNow = LET.MPIRANK == 0;
  pass.printNow = LET.MPIRANK == 0;
  traversal.printNow = LET.MPIRANK == 0;
  LET.printNow = LET.MPIRANK == 0;

  bodies.resize(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
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

  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  Box box = boundbox.bounds2box(localBounds);
  tree.buildTree(bodies, cells, box);
  pass.upwardPass(cells);
  LET.setLET(cells,localBounds,cycle);
  LET.commBodies();
  LET.commCells();
  traversal.dualTreeTraversal(cells, cells, cycle);
  for (int irank=1; irank<LET.MPISIZE; irank++) {
    LET.getLET(jcells,(LET.MPIRANK+irank)%LET.MPISIZE);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }

  pass.downwardPass(cells);
  LET.unpartition(bodies);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->ICELL = B->IBODY;
  }
  Bodies buffer = bodies;
  bodies = sort.sortBodies(buffer);

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
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
