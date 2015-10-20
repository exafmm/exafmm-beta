#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include <tbb/task_scheduler_init.h>

real_t cycle;
Bounds localBoundsB;
Bounds localBoundsV;
Bodies bbodies;
Bodies vbodies;
Cells bcells;
Cells vcells;

Args *args;
Logger *logger;
Sort *sort;
BoundBox *boundbox;
BuildTree *tree;
UpDownPass *pass;
Traversal *traversal;
LocalEssentialTree *LET;
tbb::task_scheduler_init * task_init;

void log_initialize() {
  logger->verbose = boundbox->verbose = tree->verbose = pass->verbose
    = traversal->verbose = LET->verbose = args->verbose &= LET->mpirank == 0;
  logger->printTitle("FMM Parameters");
  args->print(logger->stringLength, P);
  logger->printTitle("FMM Profiling");
  logger->startTimer("Total FMM");
  logger->startPAPI();
}

void log_finalize() {
  logger->stopPAPI();
  logger->stopTimer("Total FMM");
  logger->printTitle("Total runtime");
  logger->printTime("Total FMM");
  traversal->printTraversalData();
}

extern "C" void FMM_Init(double eps2, int ncrit, int nworkers,
			 int nb, double * xb, double * yb, double * vb,
			 int nv, double * xv, double * yv, double * vv) {
  const int nspawn = 1000;
  const int images = 0;
  const real_t theta = 0.4;
  args = new Args;
  logger = new Logger;
  sort = new Sort;
  boundbox = new BoundBox(nspawn);
  tree = new BuildTree(ncrit, nspawn);
  pass = new UpDownPass(theta, eps2);
  traversal = new Traversal(nspawn, images, eps2);
  LET = new LocalEssentialTree(images);
  task_init = new tbb::task_scheduler_init(nworkers);

  args->theta = theta;
  args->ncrit = ncrit;
  args->nspawn = nspawn;
  args->images = images;
  args->mutual = 1;
  args->distribution = "external";
#if AUTO
  traversal->timeKernels();
#endif
  bbodies.resize(nb);
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    B->X[0] = xb[i];
    B->X[1] = yb[i];
    B->SRC  = vb[i];
    B->IBODY = i;
  }
  vbodies.resize(nv);
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    B->X[0] = xv[i];
    B->X[1] = yv[i];
    B->SRC  = vv[i];
    B->IBODY = i;
  }
}

extern "C" void FMM_Finalize() {
  boundbox->writeTime(LET->mpirank);
  tree->writeTime(LET->mpirank);
  pass->writeTime(LET->mpirank);
  traversal->writeTime(LET->mpirank);
  LET->writeTime(LET->mpirank);
  LET->writeSendCount();
  logger->resetTimer();
  boundbox->resetTimer();
  tree->resetTimer();
  pass->resetTimer();
  traversal->resetTimer();
  LET->resetTimer();
  delete args;
  delete logger;
  delete sort;
  delete boundbox;
  delete tree;
  delete pass;
  delete traversal;
  delete LET;
  delete task_init;
}

extern "C" void FMM_Partition(int & nb, double * xb, double * yb, double * vb,
                              int & nv, double * xv, double * yv, double * vv) {
  Bounds localBounds = boundbox->getBounds(bbodies);
  localBounds = boundbox->getBounds(vbodies,localBounds);
  Bounds globalBounds = LET->allreduceBounds(localBounds);
  cycle = std::max(globalBounds.Xmax[0] - globalBounds.Xmin[0],
		   globalBounds.Xmax[1] - globalBounds.Xmin[1]);
  localBounds = LET->partition(bbodies,globalBounds);
  bbodies = sort->sortBodies(bbodies);
  bbodies = LET->commBodies(bbodies);
  LET->partition(vbodies,globalBounds);
  vbodies = sort->sortBodies(vbodies);
  vbodies = LET->commBodies(vbodies);
  nb = bbodies.size();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    xb[i] = B->X[0];
    yb[i] = B->X[1];
    vb[i] = B->SRC;
    B->IBODY = i;
  }
  nv = vbodies.size();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    xv[i] = B->X[0];
    yv[i] = B->X[1];
    vv[i] = B->SRC;
    B->IBODY = i;
  }
}

extern "C" void FMM_BuildTree() {
  localBoundsB = boundbox->getBounds(bbodies);
  bcells = tree->buildTree(bbodies, localBoundsB);
  localBoundsV = boundbox->getBounds(vbodies);
  vcells = tree->buildTree(vbodies, localBoundsV);
  boundbox->resetTimer();
  tree->resetTimer();
}

extern "C" void FMM_B2B(double * vi, double * vb, int verbose) {
  args->verbose = verbose;
  args->numBodies = bbodies.size();
  log_initialize();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    B->TRG   = 0;
    B->SRC   = vb[B->IBODY];
  }
  pass->upwardPass(bcells);
  LET->setLET(bcells, localBoundsB, cycle);
  Cells jcells = bcells;
#pragma omp parallel sections
  {
#pragma omp section
    {
      LET->commBodies();
      LET->commCells();
    }
#pragma omp section
    {
      traversal->dualTreeTraversal(bcells, jcells, cycle);
    }
  }
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells, (LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(bcells, jcells, cycle);
  }
  pass->downwardPass(bcells);
  log_finalize();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    vi[B->IBODY] += B->TRG;
  }
}

extern "C" void FMM_V2B(double * vb, double * vv, int verbose) {
  args->verbose = verbose;
  args->numBodies = bbodies.size();
  log_initialize();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    B->TRG   = 0;
    B->SRC   = 1;
  }
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    B->SRC   = vv[B->IBODY];
    B->TRG   = 0;
  }
  pass->upwardPass(bcells);
  pass->upwardPass(vcells);
  LET->setLET(vcells, localBoundsV, cycle);
  Cells jcells = vcells;
#pragma omp parallel sections
  {
#pragma omp section
    {
      LET->commBodies();
      LET->commCells();
    }
#pragma omp section
    {
      traversal->dualTreeTraversal(bcells, jcells, cycle);
    }
  }
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells, (LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(bcells, jcells, cycle);
  }
  pass->downwardPass(bcells);
  log_finalize();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    vb[B->IBODY] += B->TRG;
  }
}

extern "C" void FMM_B2V(double * vv, double * vb, int verbose) {
  args->verbose = verbose;
  args->numBodies = vbodies.size();
  log_initialize();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    B->TRG   = 0;
    B->SRC   = 1;
  }
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    B->SRC   = vb[B->IBODY];
    B->TRG   = 0;
  }
  pass->upwardPass(bcells);
  LET->setLET(bcells, localBoundsB, cycle);
  Cells jcells = bcells;
#pragma omp parallel sections
  {
#pragma omp section
    {
      LET->commBodies();
      LET->commCells();
    }
#pragma omp section
    {
      traversal->dualTreeTraversal(vcells, jcells, cycle);
    }
  }
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells, (LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(vcells, jcells, cycle);
  }
  pass->downwardPass(vcells);
  log_finalize();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    vv[B->IBODY] += B->TRG;
  }
}

extern "C" void FMM_V2V(double * vi, double * vv, int verbose) {
  args->verbose = verbose;
  args->numBodies = vbodies.size();
  log_initialize();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    B->TRG = 0;
    B->SRC = vv[B->IBODY];
  }
  pass->upwardPass(vcells);
  LET->setLET(vcells, localBoundsV, cycle);
  Cells jcells = vcells;
#pragma omp parallel sections
  {
#pragma omp section
    {
      LET->commBodies();
      LET->commCells();
    }
#pragma omp section
    {
      traversal->dualTreeTraversal(vcells, jcells, cycle);
    }
  }
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells, (LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(vcells, jcells, cycle);
  }
  pass->downwardPass(vcells);
  log_finalize();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    vi[B->IBODY] += B->TRG;
  }
}

extern "C" void FMM(int ni, double * xi, double * yi, double * vi,
                    int nj, double * xj, double * yj, double * vj,
                    int verbose) {
  args->verbose = verbose;
  args->numBodies = ni;
  log_initialize();
  Bodies bodies(ni);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0]  = xi[i];
    B->X[1]  = yi[i];
    B->TRG   = 0;
    B->SRC   = 1;
    B->IBODY = i;
  }
  Bodies jbodies(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0]  = xj[i];
    B->X[1]  = yj[i];
    B->SRC   = vj[i];
    B->TRG   = 0;
    B->IBODY = i;
  }

  Bounds localBounds = boundbox->getBounds(bodies);
  Cells cells = tree->buildTree(bodies, localBounds);
  pass->upwardPass(cells);
  localBounds = boundbox->getBounds(jbodies);
  Cells jcells = tree->buildTree(jbodies, localBounds);
  pass->upwardPass(jcells);

  LET->setLET(jcells, localBounds, cycle);
#pragma omp parallel sections
  {
#pragma omp section
    {
      LET->commBodies();
      LET->commCells();
    }
#pragma omp section
    {
      traversal->dualTreeTraversal(cells, jcells, cycle);
    }
  }
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells,(LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(cells, jcells, cycle);
  }
  pass->downwardPass(cells);
  log_finalize();
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    vi[B->IBODY] += B->TRG;
  }
}

extern "C" void Direct(int ni, double * xi, double * yi, double * vi,
		       int nj, double * xj, double * yj, double * vj) {
  Bodies bodies(ni);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = xi[i];
    B->X[1] = yi[i];
    B->TRG  = 0;
    B->SRC  = 1;
  }
  Bodies jbodies(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0] = xj[i];
    B->X[1] = yj[i];
    B->SRC  = vj[i];
  }
  for (int i=0; i<LET->mpisize; i++) {
    if (args->verbose) std::cout << "Direct loop          : " << i+1 << "/" << LET->mpisize << std::endl;
    LET->shiftBodies(jbodies);
    traversal->direct(bodies, jbodies, cycle);
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    vi[i] += B->TRG;
  }
}
