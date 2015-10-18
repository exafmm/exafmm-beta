#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "ewald.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "van_der_waals.h"
#if EXAFMM_SERIAL
#include "serialfmm.h"
#else
#include "parallelfmm.h"
#endif
#if EXAFMM_MASS
#error Turn off EXAFMM_MASS for this wrapper
#endif
using namespace exafmm;

static const double Celec = 332.0716;

Args * args;
BaseMPI * baseMPI;
BoundBox * boundBox;
BuildTree * localTree, * globalTree;
Partition * partition;
Traversal * traversal;
TreeMPI * treeMPI;
UpDownPass * upDownPass;

#if EXAFMM_SERIAL
SerialFMM * FMM;
#else
ParallelFMM * FMM;
#endif

Bodies buffer;
Bounds localBounds;
Bounds globalBounds;

extern "C" void fmm_init_(int & images, double & theta, int & verbose, int & nglobal) {
  int ncrit = 32;
  const int nspawn = 1000;
  const int threads = 1;
  const bool useRmax = true;
  const bool useRopt = true;
  kernel::eps2 = 0.0;
  kernel::setup();

  args = new Args;
  baseMPI = new BaseMPI;
  boundBox = new BoundBox(nspawn);
  localTree = new BuildTree(ncrit, nspawn);
  globalTree = new BuildTree(1, nspawn);
  partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
  traversal = new Traversal(nspawn, images);
  treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
  upDownPass = new UpDownPass(theta, useRmax, useRopt);
#if EXAFMM_SERIAL
  FMM = new SerialFMM;
#else
  FMM = new ParallelFMM;
#endif

  args->ncrit = ncrit;
  args->distribution = "external";
  args->dual = 1;
  args->graft = 1;
  args->images = images;
  args->mutual = 0;
  args->numBodies = 0;
  args->useRopt = useRopt;
  args->nspawn = nspawn;
  args->theta = theta;
  args->verbose = verbose & (baseMPI->mpirank == 0);
  args->useRmax = useRmax;
  logger::verbose = args->verbose;
  logger::printTitle("Initial Parameters");
  args->print(logger::stringLength, P);
  if (baseMPI->mpirank == 0) {
    std::cout << "PP      : " << EXAFMM_PP << std::endl;
    std::cout << "DP2P    : " << DP2P << std::endl;
    std::cout << "DM2L    : " << DM2L << std::endl;
  //std::cout << "ALPHA_M : " << ALPHA_M << std::endl;  
  //std::cout << "ALPHA_L : " << ALPHA_L << std::endl;
  }
  int numBodies = nglobal / baseMPI->mpisize;
  ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int gatherLevel = 1;
  FMM->allocate(numBodies, maxLevel, args->images);
  FMM->partitioner(gatherLevel);
}

extern "C" void fmm_finalize_() {
  FMM->deallocate();
  delete args;
  delete baseMPI;
  delete boundBox;
  delete localTree;
  delete globalTree;
  delete partition;
  delete traversal;
  delete treeMPI;
  delete upDownPass;
  delete FMM;
}

extern "C" void fmm_partition_(int & nglobal, int * icpumap, double * x, double * q,
			       double * xold, double & cycle) {
  num_threads(args->threads);
  int ix[3] = {0, 0, 0};
  FMM->R0 = 0.5 * cycle / FMM->numPartition[FMM->maxGlobLevel][0];
  for_3d FMM->RGlob[d] = FMM->R0 * FMM->numPartition[FMM->maxGlobLevel][d];
  FMM->getGlobIndex(ix,FMM->MPIRANK,FMM->maxGlobLevel);
  for_3d FMM->X0[d] = 2 * FMM->R0 * (ix[d] + .5);

  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      nlocal++;
    } else {
      icpumap[i] = 0;
    }
  }
  assert(nlocal % 3 == 0);
  args->numBodies = nlocal;
  logger::printTitle("Partition Profiling");
  Bodies bodies(args->numBodies);
  B_iter B = bodies.begin();
  int iwrap = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG = 0;
#if NOWRAP
      if(i % 3 == 0) iwrap = wrap(B->X, cycle);
      else unwrap(B->X, cycle, iwrap);
#else
      iwrap = wrap(B->X, cycle);
#endif
      B->IBODY = i | (iwrap << shift);
      B++;
    }
  }
  int b = 0;
  for (B=bodies.begin(); B!=bodies.end(); B++, b++) {
    FMM->Jbodies[b][0] = B->X[0] + FMM->RGlob[0];
    FMM->Jbodies[b][1] = B->X[1] + FMM->RGlob[1];
    FMM->Jbodies[b][2] = B->X[2] + FMM->RGlob[2];
    FMM->Jbodies[b][3] = B->SRC;
    FMM->Index[b] = B->IBODY;
    FMM->Ibodies[b][0] = 0;
    FMM->Ibodies[b][1] = 0;
    FMM->Ibodies[b][2] = 0;
    FMM->Ibodies[b][3] = 0;
  }
  FMM->numBodies = bodies.size();
  FMM->partitionComm();
  bodies.resize(FMM->numBodies);
  b = 0;
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
    B->X[0] = FMM->Jbodies[b][0] - FMM->RGlob[0];
    B->X[1] = FMM->Jbodies[b][1] - FMM->RGlob[1];
    B->X[2] = FMM->Jbodies[b][2] - FMM->RGlob[2];
    B->SRC = FMM->Jbodies[b][3];
    B->IBODY = FMM->Index[b];
    B->TRG[0] = FMM->Ibodies[b][0];
    B->TRG[1] = FMM->Ibodies[b][1];
    B->TRG[2] = FMM->Ibodies[b][2];
    B->TRG[3] = FMM->Ibodies[b][3];
  }
  localBounds = boundBox->getBounds(bodies);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  for (int i=0; i<nglobal; i++) {
    icpumap[i] = 0;
  }
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY & mask;
    int iwrap = unsigned(B->IBODY) >> shift;
    unwrap(B->X, cycle, iwrap);
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i] = B->SRC;
    xold[3*i+0] = B->TRG[0];
    xold[3*i+1] = B->TRG[1];
    xold[3*i+2] = B->TRG[2];
    icpumap[i] = 1;
  }
}

extern "C" void fmm_coulomb_(int & nglobal, int * icpumap,
			     double * x, double * q, double * p, double * f,
			     double & cycle) {
  num_threads(args->threads);
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      nlocal++;
    } else {
      icpumap[i] = 0;
    }
  }
  assert(nlocal % 3 == 0);
  args->numBodies = nlocal;
  logger::printTitle("FMM Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startPAPI();

  Bodies bodies(args->numBodies);
  B_iter B = bodies.begin();
  int iwrap = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG = 0;
#if NOWRAP
      if(i % 3 == 0) iwrap = wrap(B->X, cycle);
      else unwrap(B->X, cycle, iwrap);
#else
      iwrap = wrap(B->X, cycle);
#endif
      B->IBODY = i | (iwrap << shift);
      B++;
    }
  }
  int b = 0;
  for (B=bodies.begin(); B!=bodies.end(); B++, b++) {
    FMM->Jbodies[b][0] = B->X[0] + FMM->RGlob[0];
    FMM->Jbodies[b][1] = B->X[1] + FMM->RGlob[1];
    FMM->Jbodies[b][2] = B->X[2] + FMM->RGlob[2];
    FMM->Jbodies[b][3] = B->SRC;
    FMM->Index[b] = B->IBODY;
    FMM->Ibodies[b][0] = 0;
    FMM->Ibodies[b][1] = 0;
    FMM->Ibodies[b][2] = 0;
    FMM->Ibodies[b][3] = 0;
  }
  FMM->sortBodies();
  FMM->buildTree();
  FMM->upwardPass();
#if EXAFMM_SERIAL
#else
  FMM->P2PSend();
  FMM->P2PRecv();
  for( int lev=FMM->maxLevel; lev>0; lev-- ) {
    FMM->M2LSend(lev);
    FMM->M2LRecv(lev);
  }
  FMM->rootGather();
  FMM->globM2M();
  FMM->globM2L();
#endif
  FMM->periodicM2L();
#if EXAFMM_SERIAL
#else
  FMM->globL2L();
#endif
  FMM->downwardPass();

  Cells cells = localTree->buildTree(bodies, buffer, localBounds);
  upDownPass->upwardPass(cells);
  treeMPI->allgatherBounds(localBounds);
  treeMPI->setLET(cells, cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->initListCount(cells);
  traversal->initWeight(cells);
  traversal->traverse(cells, cells, cycle, args->dual, args->mutual);
  Cells jcells;
  if (args->graft) {
    treeMPI->linkLET();
    Bodies gbodies = treeMPI->root2body();
    jcells = globalTree->buildTree(gbodies, buffer, globalBounds);
    treeMPI->attachRoot(jcells);
    traversal->traverse(cells, jcells, cycle, args->dual, false);
  } else {
    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      treeMPI->getLET(jcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
      traversal->traverse(cells, jcells, cycle, args->dual, false);
    }
  }
  upDownPass->downwardPass(cells);

#if 1
  b = 0;
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
    B->X[0] = FMM->Jbodies[b][0] - FMM->RGlob[0];
    B->X[1] = FMM->Jbodies[b][1] - FMM->RGlob[1];
    B->X[2] = FMM->Jbodies[b][2] - FMM->RGlob[2];
    B->SRC = FMM->Jbodies[b][3];
    B->IBODY = FMM->Index[b];
    B->TRG[0] = FMM->Ibodies[b][0];
    B->TRG[1] = FMM->Ibodies[b][1];
    B->TRG[2] = FMM->Ibodies[b][2];
    B->TRG[3] = FMM->Ibodies[b][3];
  }
#endif

  vec3 localDipole = upDownPass->getDipole(bodies,0);
  vec3 globalDipole = baseMPI->allreduceVec3(localDipole);
  int numBodies = baseMPI->allreduceInt(bodies.size());
  upDownPass->dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger::stopPAPI();
  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");

  b = 0;
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
    int i = B->IBODY & mask;
    p[i]     += B->TRG[0] * B->SRC * Celec;
    f[3*i+0] += B->TRG[1] * B->SRC * Celec;
    f[3*i+1] += B->TRG[2] * B->SRC * Celec;
    f[3*i+2] += B->TRG[3] * B->SRC * Celec;
  }
  bodies = treeMPI->getRecvBodies();
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY & mask;
    int iwrap = unsigned(B->IBODY) >> shift;
    unwrap(B->X, cycle, iwrap);
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i] = B->SRC;
    assert(icpumap[i] == 0);
    icpumap[i] = 2;
  }
}

extern "C" void ewald_coulomb_(int & nglobal, int * icpumap, double * x, double * q, double * p, double * f,
			       int & ksize, double & alpha, double & sigma, double & cutoff, double & cycle) {
  num_threads(args->threads);
  Ewald * ewald = new Ewald(ksize, alpha, sigma, cutoff, cycle);
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      nlocal++;
    } else {
      icpumap[i] = 0;
    } 
  }
  args->numBodies = nlocal;
  logger::printTitle("Ewald Parameters");
  args->print(logger::stringLength, P);
  ewald->print(logger::stringLength);
  logger::printTitle("Ewald Profiling");
  logger::startTimer("Total Ewald");
  logger::startPAPI();
  Bodies bodies(nlocal);
  B_iter B = bodies.begin();
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG = 0;
      int iwrap = wrap(B->X, cycle);
      B->IBODY = i | (iwrap << shift);
      B++;
    }
  }
  localBounds = boundBox->getBounds(bodies);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  localBounds = partition->octsection(bodies,globalBounds);
  bodies = treeMPI->commBodies(bodies);
  Cells cells = localTree->buildTree(bodies, buffer, localBounds);
  Bodies jbodies = bodies;
  for (int i=0; i<baseMPI->mpisize; i++) {
    if (args->verbose) std::cout << "Ewald loop           : " << i+1 << "/" << baseMPI->mpisize << std::endl;
    treeMPI->shiftBodies(jbodies);
    localBounds = boundBox->getBounds(jbodies);
    Cells jcells = localTree->buildTree(jbodies, buffer, localBounds);
    ewald->wavePart(bodies, jbodies);
    ewald->realPart(cells, jcells);
  }
  ewald->selfTerm(bodies);
  partition->unpartition(bodies);
  bodies = treeMPI->commBodies(bodies);
  logger::stopPAPI();
  logger::stopTimer("Total Ewald");
  logger::printTitle("Total runtime");
  logger::printTime("Total Ewald");

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY & mask;
    p[i]     += B->TRG[0] * B->SRC * Celec;
    f[3*i+0] += B->TRG[1] * B->SRC * Celec;
    f[3*i+1] += B->TRG[2] * B->SRC * Celec;
    f[3*i+2] += B->TRG[3] * B->SRC * Celec;
  }
  delete ewald;
}

extern "C" void direct_coulomb_(int & nglobal, int * icpumap, double * x, double * q, double * p, double * f, double & cycle) {
  num_threads(args->threads);
  logger::startTimer("Direct Coulomb");
  int images = args->images;
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  real_t Xperiodic[3];
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      real_t pp = 0, fx = 0, fy = 0, fz = 0;
      for (int ix=-prange; ix<=prange; ix++) {
	for (int iy=-prange; iy<=prange; iy++) {
	  for (int iz=-prange; iz<=prange; iz++) {
	    Xperiodic[0] = ix * cycle;
	    Xperiodic[1] = iy * cycle;
	    Xperiodic[2] = iz * cycle;
	    for (int j=0; j<nglobal; j++) {
	      vec3 dX;
	      for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d] - Xperiodic[d];
	      real_t R2 = norm(dX);
	      real_t invR = 1 / std::sqrt(R2);
	      if (R2 == 0) invR = 0;
	      real_t invR3 = q[j] * invR * invR * invR;
	      pp += q[j] * invR;
	      fx += dX[0] * invR3;
	      fy += dX[1] * invR3;
	      fz += dX[2] * invR3;
	    }
	  }
	}
      }
      p[i] += pp * q[i] * Celec;
      f[3*i+0] -= fx * q[i] * Celec;
      f[3*i+1] -= fy * q[i] * Celec;
      f[3*i+2] -= fz * q[i] * Celec;
    }
  }
  real_t dipole[3] = {0, 0, 0};
  for (int i=0; i<nglobal; i++) {
    for (int d=0; d<3; d++) dipole[d] += x[3*i+d] * q[i];
  }
  real_t norm = 0;
  for (int d=0; d<3; d++) {
    norm += dipole[d] * dipole[d];
  }
  real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle) * Celec;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      p[i] -= coef * norm / nglobal;
      f[3*i+0] -= coef * dipole[0] * q[i];
      f[3*i+1] -= coef * dipole[1] * q[i];
      f[3*i+2] -= coef * dipole[2] * q[i];
    }
  }
  logger::stopTimer("Direct Coulomb");
}

extern "C" void coulomb_exclusion_(int & nglobal, int * icpumap,
				   double * x, double * q, double * p, double * f,
				   double & cycle, int * numex, int * natex) {
  num_threads(args->threads);
  logger::startTimer("Coulomb Exclusion");
  for (int i=0, ic=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      real_t pp = 0, fx = 0, fy = 0, fz = 0;
      for (int jc=0; jc<numex[i]; jc++, ic++) {
	int j = natex[ic]-1;
	vec3 dX;
	for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
        wrap(dX, cycle);
	real_t R2 = norm(dX);
	real_t invR = 1 / std::sqrt(R2);
	if (R2 == 0) invR = 0;
	real_t invR3 = q[j] * invR * invR * invR;
	pp += q[j] * invR;
	fx += dX[0] * invR3;
	fy += dX[1] * invR3;
	fz += dX[2] * invR3;
      }
      p[i] -= pp * q[i] * Celec;
      f[3*i+0] += fx * q[i] * Celec; 
      f[3*i+1] += fy * q[i] * Celec;
      f[3*i+2] += fz * q[i] * Celec;
    } else {
      ic += numex[i];
    }
  }
  logger::stopTimer("Coulomb Exclusion");
}

extern "C" void fmm_vanderwaals_(int & nglobal, int * icpumap, int * atype,
				 double * x, double * p, double * f,
				 double & cuton, double & cutoff, double & cycle,
				 int & numTypes, double * rscale, double * gscale, double * fgscale) {
  num_threads(args->threads);
  VanDerWaals * VDW = new VanDerWaals(cuton, cutoff, cycle, numTypes, rscale, gscale, fgscale);
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) nlocal++;
  }
  args->numBodies = nlocal;
  logger::printTitle("VdW Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("VdW Profiling");
  logger::startTimer("Total VdW");
  logger::startPAPI();
  Bodies bodies(nlocal);
  B_iter B = bodies.begin();
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = atype[i] - .5;
      B->TRG = 0;
      int iwrap = wrap(B->X, cycle);
      B->IBODY = i | (iwrap << shift);
      B++;
    }
  }
  localBounds = boundBox->getBounds(bodies);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  localBounds = partition->octsection(bodies,globalBounds);
  bodies = treeMPI->commBodies(bodies);
  Cells cells = localTree->buildTree(bodies, buffer, localBounds);
  upDownPass->upwardPass(cells);
  treeMPI->allgatherBounds(localBounds);
  treeMPI->setLET(cells, cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  VDW->evaluate(cells, cells);
  Cells jcells;
  for (int irank=1; irank<baseMPI->mpisize; irank++) {
    treeMPI->getLET(jcells,(baseMPI->mpirank+irank)%baseMPI->mpisize);
    VDW->evaluate(cells, jcells);
  }
  partition->unpartition(bodies);
  bodies = treeMPI->commBodies(bodies);
  logger::stopPAPI();
  logger::stopTimer("Total VdW");
  logger::printTitle("Total runtime");
  logger::printTime("Total VdW");

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY & mask;
    p[i]     += B->TRG[0];
    f[3*i+0] += B->TRG[1];
    f[3*i+1] += B->TRG[2];
    f[3*i+2] += B->TRG[3];
  }
  delete VDW;
}

extern "C" void direct_vanderwaals_(int & nglobal, int * icpumap, int * atype,
				    double * x, double * p, double * f,
				    double & cuton, double & cutoff, double & cycle,
				    int & numTypes, double * rscale, double * gscale, double * fgscale) {
  num_threads(args->threads);
  logger::startTimer("Direct VdW");
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      int atypei = atype[i]-1;
      kreal_t pp = 0, fx = 0, fy = 0, fz = 0;
      for (int j=0; j<nglobal; j++) {
	vec3 dX;
	for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
	wrap(dX, cycle);
	real_t R2 = norm(dX);
	if (R2 != 0) {
	  int atypej = atype[j]-1;
	  real_t rs = rscale[atypei*numTypes+atypej];
	  real_t gs = gscale[atypei*numTypes+atypej];
	  real_t fgs = fgscale[atypei*numTypes+atypej];
	  real_t R2s = R2 * rs;
	  real_t invR2 = 1.0 / R2s;
	  real_t invR6 = invR2 * invR2 * invR2;
	  real_t cuton2 = cuton * cuton;
	  real_t cutoff2 = cutoff * cutoff;
          if (R2 < cutoff2) {
            real_t tmp = 0, dtmp = 0;
            if (cuton2 < R2) {
              real_t tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
              real_t tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
              tmp = invR6 * (invR6 - 1) * tmp2;
              dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
                - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
            } else {
              tmp = invR6 * (invR6 - 1);
              dtmp = invR2 * invR6 * (2 * invR6 - 1);
            }
            dtmp *= fgs;
            pp += gs * tmp;
            fx += dX[0] * dtmp;
            fy += dX[1] * dtmp;
            fz += dX[2] * dtmp;
          }
	}
      }
      p[i] += pp;
      f[3*i+0] -= fx;
      f[3*i+1] -= fy;
      f[3*i+2] -= fz;
    }
  }
  logger::stopTimer("Direct VdW");
}

extern "C" void vanderwaals_exclusion_(int & nglobal, int * icpumap, int * atype,
				       double * x, double * p, double * f,
				       double & cuton, double & cutoff, double & cycle,
				       int & numTypes, double * rscale, double * gscale,
				       double * fgscale, int * numex, int * natex) {
  num_threads(args->threads);
  logger::startTimer("VdW Exclusion");
  for (int i=0, ic=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      int atypei = atype[i]-1;
      for (int jc=0; jc<numex[i]; jc++, ic++) {
	int j = natex[ic]-1;
	vec3 dX;
	for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
        wrap(dX, cycle);
	real_t R2 = norm(dX);
        if (R2 != 0) {
          int atypej = atype[j]-1;
          real_t rs = rscale[atypei*numTypes+atypej];
          real_t gs = gscale[atypei*numTypes+atypej];
          real_t fgs = fgscale[atypei*numTypes+atypej];
          real_t R2s = R2 * rs;
          real_t invR2 = 1.0 / R2s;
          real_t invR6 = invR2 * invR2 * invR2;
          real_t cuton2 = cuton * cuton;
          real_t cutoff2 = cutoff * cutoff;
          if (R2 < cutoff2) {
            real_t tmp = 0, dtmp = 0;
            if (cuton2 < R2) {
              real_t tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
              real_t tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
              tmp = invR6 * (invR6 - 1) * tmp2;
              dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
                - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
            } else {
              tmp = invR6 * (invR6 - 1);
              dtmp = invR2 * invR6 * (2 * invR6 - 1);
            }
            dtmp *= fgs;
            p[i] -= gs * tmp;
            f[3*i+0] += dX[0] * dtmp;
            f[3*i+1] += dX[1] * dtmp;
            f[3*i+2] += dX[2] * dtmp;
          }
        }
      }
    } else {
      ic += numex[i];
    }
  }
  logger::stopTimer("VdW Exclusion");
}
