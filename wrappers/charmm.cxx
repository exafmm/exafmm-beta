#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "ewald.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"

Args *args;
Logger *logger;
Sort *sort;
Bounds localBounds;
BoundBox *boundbox;
BuildTree *tree;
UpDownPass *pass;
Traversal *traversal;
LocalEssentialTree *LET;

extern "C" void fmm_init_(int * images) {
  const int ncrit = 32;
  const int nspawn = 1000;
  const real_t theta = 0.4;
  args = new Args;
  logger = new Logger;
  sort = new Sort;
  boundbox = new BoundBox(nspawn);
  tree = new BuildTree(ncrit, nspawn);
  pass = new UpDownPass(theta);
  traversal = new Traversal(nspawn, *images);
  LET = new LocalEssentialTree(*images);

  args->theta = theta;
  args->ncrit = ncrit;
  args->nspawn = nspawn;
  args->images = *images;
  args->mutual = 0;
  args->verbose = 1;
  args->distribution = "external";
  args->verbose &= LET->mpirank == 0;
  if (args->verbose) {
    logger->verbose = true;
    boundbox->verbose = true;
    tree->verbose = true;
    pass->verbose = true;
    traversal->verbose = true;
    LET->verbose = true;
  }
  logger->printTitle("Initial Parameters");
  args->print(logger->stringLength, P, LET->mpirank);
}

extern "C" void fmm_partition_(int * nglobal, int * icpumap, double * x, double * q, double * cycle) {
  logger->printTitle("Partition Profiling");
  int nlocal = 0;
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) nlocal++;
  }
  Bodies bodies(nlocal);
  B_iter B = bodies.begin();
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->IBODY = i;
      B++;
    }
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    if( B->X[0] < -*cycle/2 ) B->X[0] += *cycle;
    if( B->X[1] < -*cycle/2 ) B->X[1] += *cycle;
    if( B->X[2] < -*cycle/2 ) B->X[2] += *cycle;
    if( B->X[0] >  *cycle/2 ) B->X[0] -= *cycle;
    if( B->X[1] >  *cycle/2 ) B->X[1] -= *cycle;
    if( B->X[2] >  *cycle/2 ) B->X[2] -= *cycle;
  }
  localBounds = boundbox->getBounds(bodies);
  Bounds globalBounds = LET->allreduceBounds(localBounds);
  localBounds = LET->partition(bodies,globalBounds);
  bodies = sort->sortBodies(bodies);
  bodies = LET->commBodies(bodies);
  Cells cells = tree->buildTree(bodies, localBounds);
  pass->upwardPass(cells);
  for (int i=0; i<*nglobal; i++) {
    icpumap[i] = 0;
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i]     = B->SRC;
    icpumap[i] = 1;
  }
}

extern "C" void fmm_(int * nglobal, int * icpumap, double * x, double * q, double * p, double * f, double * cycle) {
  int nlocal = 0;
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) nlocal++;
  }
  args->numBodies = nlocal;
  logger->printTitle("FMM Parameters");
  args->print(logger->stringLength, P, LET->mpirank);
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  logger->printTitle("FMM Profiling");
  logger->startTimer("Total FMM");
  logger->startPAPI();
  Bodies bodies(nlocal);
  B_iter B = bodies.begin();
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG[0] = p[i];
      B->TRG[1] = f[3*i+0];
      B->TRG[2] = f[3*i+1];
      B->TRG[3] = f[3*i+2];
      B->IBODY = i;
      B++;
    }
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    if( B->X[0] < -*cycle/2 ) B->X[0] += *cycle;
    if( B->X[1] < -*cycle/2 ) B->X[1] += *cycle;
    if( B->X[2] < -*cycle/2 ) B->X[2] += *cycle;
    if( B->X[0] >  *cycle/2 ) B->X[0] -= *cycle;
    if( B->X[1] >  *cycle/2 ) B->X[1] -= *cycle;
    if( B->X[2] >  *cycle/2 ) B->X[2] -= *cycle;
  }
  Cells cells = tree->buildTree(bodies, localBounds);
  pass->upwardPass(cells);
  LET->setLET(cells, localBounds, *cycle);
  LET->commBodies();
  LET->commCells();
  traversal->dualTreeTraversal(cells, cells, *cycle, args->mutual);
  Cells jcells;
  for (int irank=1; irank<LET->mpisize; irank++) {
    LET->getLET(jcells,(LET->mpirank+irank)%LET->mpisize);
    traversal->dualTreeTraversal(cells, jcells, *cycle);
  }
  pass->downwardPass(cells);
  vec3 localDipole = pass->getDipole(bodies,0);
  vec3 globalDipole = LET->allreduceVec3(localDipole);
  int numBodies = LET->allreduceInt(bodies.size());
  pass->dipoleCorrection(bodies, globalDipole, numBodies, *cycle);
  logger->stopPAPI();
  logger->stopTimer("Total FMM");
  logger->printTitle("Total runtime");
  logger->printTime("Total FMM");

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
}

extern "C" void exclusion_(int * nglobal, int * icpumap, double * x, double * q, double * p, double * f, double * cycle, int * numex, int * natex) {
  logger->printTitle("Exclusion Profiling");
  logger->startTimer("Total Exclusion");
  int ic = 0;
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) {
      for (int jc=0; jc<numex[i]; jc++, ic++) {
	int j = natex[ic]-1;
	float dx = x[3*i+0] - x[3*j+0];
	float dy = x[3*i+1] - x[3*j+1];
	float dz = x[3*i+2] - x[3*j+2];
	if (dx < -*cycle/2) dx += *cycle;
	if (dy < -*cycle/2) dy += *cycle;
	if (dz < -*cycle/2) dz += *cycle;
	if (dx >  *cycle/2) dx -= *cycle;
	if (dy >  *cycle/2) dy -= *cycle;
	if (dz >  *cycle/2) dz -= *cycle;
	float R2 = dx * dx + dy * dy + dz * dz;
	float invR = 1 / std::sqrt(R2);
	if (R2 == 0) invR = 0;
	float invR3 = invR * invR * invR;
        float qinvR3 = q[j] * invR3;
	p[i] -= q[j] * invR;
	f[3*i+0] += dx * qinvR3;
	f[3*i+1] += dy * qinvR3;
	f[3*i+2] += dz * qinvR3;
	p[j] -= q[i] * invR;
	qinvR3 = q[i] * invR3;
	f[3*j+0] -= dx * qinvR3;
	f[3*j+1] -= dy * qinvR3;
	f[3*j+2] -= dz * qinvR3;
      }
    } else {
      ic += numex[i];
    }
  }
  logger->stopTimer("Total Exclusion");
  logger->printTime("Total Exclusion");
}

extern "C" void ewald_(int * nglobal, int * icpumap, double * x, double * q, double * p, double * f,
                       int * ksize, double * alpha, double * sigma, double * cutoff, double * cycle) {
  Ewald * ewald = new Ewald(*ksize, *alpha, *sigma, *cutoff, *cycle);
  if (args->verbose) ewald->verbose = true;
  int nlocal = 0;
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) nlocal++;
  }
  args->numBodies = nlocal;
  logger->printTitle("Ewald Parameters");
  args->print(logger->stringLength, P, LET->mpirank);
  ewald->print(logger->stringLength);
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
  logger->printTitle("Ewald Profiling");
  logger->startTimer("Total Ewald");
  logger->startPAPI();
  Bodies bodies(nlocal);
  B_iter B = bodies.begin();
  for (int i=0; i<*nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG[0] = p[i];
      B->TRG[1] = f[3*i+0];
      B->TRG[2] = f[3*i+1];
      B->TRG[3] = f[3*i+2];
      B->IBODY = i;
      B++;
    }
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    if( B->X[0] < -*cycle/2 ) B->X[0] += *cycle;
    if( B->X[1] < -*cycle/2 ) B->X[1] += *cycle;
    if( B->X[2] < -*cycle/2 ) B->X[2] += *cycle;
    if( B->X[0] >  *cycle/2 ) B->X[0] -= *cycle;
    if( B->X[1] >  *cycle/2 ) B->X[1] -= *cycle;
    if( B->X[2] >  *cycle/2 ) B->X[2] -= *cycle;
  }
  Cells cells = tree->buildTree(bodies, localBounds);
  Bodies jbodies = bodies;
  for (int i=0; i<LET->mpisize; i++) {
    if (args->verbose) std::cout << "Ewald loop           : " << i+1 << "/" << LET->mpisize << std::endl;
    LET->shiftBodies(jbodies);
    localBounds = boundbox->getBounds(jbodies);
    Cells jcells = tree->buildTree(jbodies, localBounds);
    ewald->wavePart(bodies, jbodies);
    ewald->realPart(cells, jcells);
  }
  ewald->selfTerm(bodies);
  logger->stopPAPI();
  logger->stopTimer("Total Ewald");
  logger->printTitle("Total runtime");
  logger->printTime("Total Ewald");

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
  delete ewald;
}
