#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include <fstream>
using namespace exafmm;

vec3 cycles;
Bodies buffer;
Bounds globalBounds;
Bodies bbodies;
Bodies vbodies;
Cells bcells;
Cells vcells;
Args * args;
BaseMPI * baseMPI;
BoundBox * boundBox;
BuildTree * localTree, * globalTree;
Partition * partition;
Traversal * traversal;
TreeMPI * treeMPI;
UpDownPass * upDownPass;
std::vector<std::vector<double> > nearGauss;
std::vector<int> patches;
void log_initialize() {
  args->verbose &= baseMPI->mpirank == 0;
  logger::verbose = args->verbose;
  logger::printTitle("FMM Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startPAPI();
}

void log_finalize() {
  logger::stopPAPI();
  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");
}

extern "C" void FMM_Init(double eps2, double kreal, double kimag, int ncrit, int threads,
			                   int nb, double * xb, double * yb, double * zb, int* patchids, 
                         std::vector<std::vector<double> > nearGaussPoints, int nhdgqp, int nipp, double nearpd, std::vector<double> ws, std::vector<std::vector<double> > ipolator_near) {
  const int nspawn = 1000;
  const int images = 0;
  const double theta = 0.4;
  const bool useRmax = false;
  const bool useRopt = false;
  const bool verbose = false;  
  kernel::eps2 = eps2;
  kernel::wavek = complex_t(kreal, kimag);
  kernel::nhdgqp = nhdgqp;
  kernel::nipp = nipp;
  kernel::nearpd = nearpd;
  nearGauss = nearGaussPoints;

#if EXAFMM_SINGLE
  kernel::ws.resize(ws.size());
  std::copy(ws.begin(), ws.end(), kernel::ws.begin());
  kernel::ipolator_near.resize(ipolator_near.size());
  for (int i = 0; i < ipolator_near.size(); ++i) {
    kernel::ipolator_near[i].resize(ipolator_near[i].size());
    std::copy(ipolator_near[i].begin(), ipolator_near[i].end(), kernel::ipolator_near[i].begin());
  }
#else
  kernel::ws = ws;
  kernel::ipolator_near = ipolator_near;
#endif
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
  num_threads(threads);
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
  bbodies.resize(nb);
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    int patch = patchids[i];
    patches.push_back(patch);
    B->X[0] = xb[i];
    B->X[1] = yb[i];
    B->X[2] = zb[i];
    B->PATCH = patch;
    B->POINT_LOC = i%nipp;
    for (int j = 0; j < nhdgqp; ++j) { 
      B->GAUSS_NEAR[j][0] = nearGaussPoints[patch*nhdgqp+j][0];
      B->GAUSS_NEAR[j][1] = nearGaussPoints[patch*nhdgqp+j][1];
      B->GAUSS_NEAR[j][2] = nearGaussPoints[patch*nhdgqp+j][2]; 
    }
    B->WEIGHT = 1;
  }
  vbodies.resize(nb);
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    int patch = patchids[i];
    B->X[0] = xb[i];
    B->X[1] = yb[i];
    B->X[2] = zb[i];
    B->PATCH = patch;
    B->POINT_LOC = i%nipp;
    for (int j = 0; j < nhdgqp; ++j) { 
      B->GAUSS_NEAR[j][0] = nearGaussPoints[patch*nhdgqp+j][0];
      B->GAUSS_NEAR[j][1] = nearGaussPoints[patch*nhdgqp+j][1];
      B->GAUSS_NEAR[j][2] = nearGaussPoints[patch*nhdgqp+j][2]; 
    }
    B->WEIGHT = 1;
  }
}

extern "C" void FMM_Finalize() {
  delete args;
  delete baseMPI;
  delete boundBox;
  delete localTree;
  delete globalTree;
  delete partition;
  delete traversal;
  delete treeMPI;
  delete upDownPass;
}

extern "C" void FMM_Partition(int & nb, double * xb, double * yb, double * zb) {
 logger::printTitle("Partition Profiling");
  Bounds localBounds = boundBox->getBounds(bbodies);
  localBounds = boundBox->getBounds(vbodies, localBounds);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  cycles = globalBounds.Xmax - globalBounds.Xmin;
  partition->bisection(bbodies, globalBounds);
  bbodies = treeMPI->commBodies(bbodies);
  partition->bisection(vbodies, globalBounds);
  vbodies = treeMPI->commBodies(vbodies);
  treeMPI->allgatherBounds(localBounds);

  nb = bbodies.size();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    xb[i] = B->X[0];
    yb[i] = B->X[1];
    zb[i] = B->X[2];
    B->IBODY = i;
  }
  //nv = vbodies.size();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    B->IBODY = i;
  }
}

extern "C" void FMM_BuildTree() {
  Bounds localBoundsB = boundBox->getBounds(bbodies);
  bcells = localTree->buildTree(bbodies, buffer, localBoundsB);
  Bounds localBoundsV = boundBox->getBounds(vbodies);
  vcells = localTree->buildTree(vbodies, buffer, localBoundsV);
}

extern "C" void FMM_B2B(std::complex<double>* vi, std::complex<double>* vb, std::complex<double>* wb, bool verbose) { 
  args->verbose = verbose;
  log_initialize();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    B->SRC    = 1;   
    B->QWEIGHT = 1;  
    B->TRG    = 0;
    B->ICELL = 0;  
  }
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    B->SRC    = vb[B->IBODY];   
    B->QWEIGHT = wb[B->IBODY];  
    B->TRG    = 0;
    B->ICELL = 0;  
  }
  upDownPass->upwardPass(bcells);
  upDownPass->upwardPass(vcells);
  treeMPI->setLET(vcells, cycles);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->initListCount(bcells);
  traversal->initWeight(bcells);
  traversal->traverse(bcells, vcells, cycles, args->dual, args->mutual);
  if (baseMPI->mpisize > 1) {
    if (args->graft) {
      treeMPI->linkLET();
      Bodies gbodies = treeMPI->root2body();
      Cells jcells = globalTree->buildTree(gbodies, buffer, globalBounds);
      treeMPI->attachRoot(jcells);
      traversal->traverse(bcells, jcells, cycles, args->dual, false);
    } else {
      for (int irank=0; irank<baseMPI->mpisize; irank++) {
  Cells jcells;
  treeMPI->getLET(jcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
  traversal->traverse(bcells, jcells, cycles, args->dual, false);
      }
    }
  }
  upDownPass->downwardPass(bcells);
  if(verbose) {
    localTree->printTreeData(bcells); 
    traversal->printTraversalData();
  }
  log_finalize();

  //size_t s = bbodies.size();
  //for (int i = 0; i < s; ++i) vi[i] = 0;  
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    vi[B->IBODY] += B->TRG[0];
  }
}

extern "C" void Direct(int ni, double * xi, double * yi, double * zi, std::complex<double>* vi,
		                   int nj, double * xj, double * yj, double * zj, std::complex<double>* vj, std::complex<double>* wj) {
  Bodies bodies(ni);
  int nhdgqp = kernel::nhdgqp;
  int nipp =   kernel::nipp ;
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    int patch = patches[i];
    B->X[0] = xi[i];
    B->X[1] = yi[i];
    B->X[2] = zi[i];
    B->TRG  = 0;
    B->SRC  = 1;
    B->QWEIGHT = 1;
    B->PATCH = patch;
    B->POINT_LOC = i%nipp;
    for (int j = 0; j < nhdgqp; ++j) { 
      B->GAUSS_NEAR[j][0] = nearGauss[patch*nhdgqp+j][0];
      B->GAUSS_NEAR[j][1] = nearGauss[patch*nhdgqp+j][1];
      B->GAUSS_NEAR[j][2] = nearGauss[patch*nhdgqp+j][2]; 
    }
  }
  Bodies jbodies(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    int patch = patches[i];
    B->X[0] = xj[i];
    B->X[1] = yj[i];
    B->X[2] = zj[i];
    B->SRC  = vj[i] ;
    B->QWEIGHT = wj[i];    
    B->PATCH = patch;
    B->POINT_LOC = i%nipp;
    for (int j = 0; j < nhdgqp; ++j) { 
      B->GAUSS_NEAR[j][0] = nearGauss[patch*nhdgqp+j][0];
      B->GAUSS_NEAR[j][1] = nearGauss[patch*nhdgqp+j][1];
      B->GAUSS_NEAR[j][2] = nearGauss[patch*nhdgqp+j][2]; 
    }    
  }  
  for (int irank=0; irank<baseMPI->mpisize; irank++) {
    if (args->verbose) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
    treeMPI->shiftBodies(jbodies);
    traversal->direct(bodies, jbodies, cycles);
  }
  size_t s = bbodies.size();
  for (int i = 0; i < s; ++i) vi[i] = 0; 
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
// #if EXAFMM_HELMHOLTZ
//     vi[i] += std::real(B->TRG[0]);
// #else
    vi[i] += B->TRG[0];
//#endif
  }
}


extern "C" void DirectAll(std::complex<double> * vi, std::complex<double>* vb, std::complex<double>* wb, bool verbose) {
  args->verbose = verbose;
  FMM_BuildTree();

  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();       
    B->SRC = 1.0;//vb[i] * wb[i];  
    B->TRG    = 0;
    B->QWEIGHT = 1.0;
  }
  Bodies jbodies(bbodies);  
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->SRC    = vb[i] * wb[i];  
    B->QWEIGHT = 1;//wb[i];  
  }
  for (int irank=0; irank<baseMPI->mpisize; irank++) {
    if (args->verbose) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
    treeMPI->shiftBodies(jbodies);
    traversal->direct(bbodies, jbodies, cycles);
  }
  size_t s = bbodies.size();
  for (int i = 0; i < s; ++i) vi[i] = 0;  
 
  // for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) { 
  //   int i = B-bbodies.begin();
  //   vi[i] = 0;
  // }
 
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    //int i = B-bbodies.begin();
    vi[B->IBODY] += B->TRG[0];
  }
}
