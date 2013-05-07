#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "localessentialtree.h"

extern "C" void FMM(int ni, double* xi, double *pi, int nj, double *xj, double * qj, double nu, double rho) {
  Bodies bodies(ni), jbodies(nj);
  Cells cells, jcells;
  Sort sort;

  const int NCRIT = 16;
  const int NSPAWN = 1000;
  const int IMAGES = 0;
  const real_t THETA = 0.6;
  const real_t cycle = 2 * M_PI;

  BoundBox boundbox(NSPAWN);
  BuildTree tree(NCRIT,NSPAWN);
  UpDownPass pass(THETA);
  Traversal traversal(NSPAWN,IMAGES);
  LocalEssentialTree LET(IMAGES);
  boundbox.verbose = LET.MPIRANK == 0;
  tree.verbose = LET.MPIRANK == 0;
  pass.verbose = LET.MPIRANK == 0;
  traversal.verbose = LET.MPIRANK == 0;
  LET.verbose = LET.MPIRANK == 0;
  traversal.NU = pass.NU = nu;
  traversal.RHO = pass.RHO = rho;

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    B->TRG[0] =  pi[i];
    B->IBODY  = i;
  }

  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0] = xj[3*i+0];
    B->X[1] = xj[3*i+1];
    B->X[2] = xj[3*i+2];
    B->SRC =  qj[i];
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
    p[i]     = B->TRG[0];
  }
}
