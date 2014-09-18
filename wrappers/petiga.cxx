#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"

real_t cycle;

Args * args;
BaseMPI * baseMPI;
BoundBox * boundBox;
BuildTree * localTree, * globalTree;
Partition * partition;
Traversal * traversal;
TreeMPI * treeMPI;
UpDownPass * upDownPass;

Bodies buffer;
Bounds localBounds;
Bounds globalBounds;
Bodies bbodies;
Bodies vbodies;

extern "C" void FMM_Init(double eps2, int ncrit, int threads,
			 int nb, double * xb, double * yb, double * zb, double * vb,
			 int nv, double * xv, double * yv, double * zv, double * vv) {
  const int nspawn = 1000;
  const int images = 0;
  const real_t theta = 0.4;
  const bool useRmax = true;
  const bool useRopt = true;
  args = new Args;
  baseMPI = new BaseMPI;
  boundBox = new BoundBox(nspawn);
  localTree = new BuildTree(ncrit, nspawn);
  globalTree = new BuildTree(1, nspawn);
  partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
  traversal = new Traversal(nspawn, images, eps2);
  treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
  upDownPass = new UpDownPass(theta, useRmax, useRopt);
  num_threads(threads);

  args->ncrit = ncrit;
  args->nspawn = nspawn;
  args->threads = threads;
  args->images = images;
  args->theta = theta;
  args->mutual = 0;
  args->verbose = 1;
  args->distribution = "external";
  args->verbose &= baseMPI->mpirank == 0;
  logger::verbose = args->verbose;
  bbodies.resize(nb);
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    B->X[0] = xb[i];
    B->X[1] = yb[i];
    B->X[2] = zb[i];
    B->SRC  = vb[i];
  }
  vbodies.resize(nv);
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    B->X[0] = xv[i];
    B->X[1] = yv[i];
    B->X[2] = zv[i];
    B->SRC  = vv[i];
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

extern "C" void FMM_Partition(int & nb, double * xb, double * yb, double * zb, double * vb,
			      int & nv, double * xv, double * yv, double * zv, double * vv) {
  logger::printTitle("Partition Profiling");
  localBounds = boundBox->getBounds(bbodies);
  localBounds = boundBox->getBounds(vbodies,localBounds);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  cycle = max(globalBounds.Xmax - globalBounds.Xmin);
  localBounds = partition->octsection(bbodies,globalBounds);
  bbodies = treeMPI->commBodies(bbodies);
  partition->octsection(vbodies,globalBounds);
  vbodies = treeMPI->commBodies(vbodies);

  nb = bbodies.size();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    xb[i] = B->X[0];
    yb[i] = B->X[1];
    zb[i] = B->X[2];
    vb[i] = B->SRC;
    B->IBODY = i;
  }
  nv = vbodies.size();
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    xv[i] = B->X[0];
    yv[i] = B->X[1];
    zv[i] = B->X[2];
    vv[i] = B->SRC;
    B->IBODY = i;
  }
}

extern "C" void FMM_Laplace(int nb, double * xb, double * yb, double * zb, double * vb,
			    int nv, double * xv, double * yv, double * zv, double * vv) {
  args->numBodies = nb;
  logger::printTitle("FMM Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startPAPI();
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B-bbodies.begin();
    B->SRC    = 1;
    B->TRG    = 0;
    B->TRG[0] = vb[i];
    B->IBODY  = i;
  }
  for (B_iter B=vbodies.begin(); B!=vbodies.end(); B++) {
    int i = B-vbodies.begin();
    B->SRC    = vv[i];
    B->TRG    = 0;
    B->IBODY  = i;
  }
  Cells bcells = localTree->buildTree(bbodies, buffer, localBounds);
  upDownPass->upwardPass(bcells);
  Cells vcells = localTree->buildTree(vbodies, buffer, localBounds);
  upDownPass->upwardPass(vcells);
  treeMPI->allgatherBounds(localBounds);
  treeMPI->setLET(vcells, cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->initWeight(bcells);
  traversal->dualTreeTraversal(bcells, vcells, cycle, args->mutual);
  if (args->graft) {
    treeMPI->linkLET();
    Bodies gbodies = treeMPI->root2body();
    vcells = globalTree->buildTree(gbodies, buffer, globalBounds);
    treeMPI->attachRoot(vcells);
    traversal->dualTreeTraversal(bcells, vcells, cycle, false);
  } else {
    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      treeMPI->getLET(vcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
      traversal->dualTreeTraversal(bcells, vcells, cycle, false);
    }
  }
  upDownPass->downwardPass(bcells);
  logger::stopPAPI();
  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");
  for (B_iter B=bbodies.begin(); B!=bbodies.end(); B++) {
    int i = B->IBODY;
    vb[i] = B->TRG[0];
  }
}

void MPI_Shift(double * var, int &nold, int mpisize, int mpirank) {
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  int nnew;
  MPI_Request sreq, rreq;
  MPI_Isend(&nold, 1, MPI_DOUBLE, irecv, 0, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(&nnew, 1, MPI_DOUBLE, isend, 0, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  double * buf = new double [nnew];
  MPI_Isend(var, nold, MPI_DOUBLE, irecv, 1, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(buf, nnew, MPI_DOUBLE, isend, 1, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  for (int i=0; i<nnew; i++) {
    var[i] = buf[i];
  }
  nold = nnew;
  delete[] buf;
}

extern "C" void Direct_Laplace(int ni, double * xi, double * yi, double * zi, double * vi,
			       int nj, double * xj, double * yj, double * zj, double * vj) {
  const int Nmax = 1000000;
  double * x2 = new double [Nmax];
  double * y2 = new double [Nmax];
  double * z2 = new double [Nmax];
  double * v2 = new double [Nmax];
  for (int i=0; i<nj; i++) {
    x2[i] = xj[i];
    y2[i] = yj[i];
    z2[i] = zj[i];
    v2[i] = vj[i];
  }
  if (baseMPI->mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
  for (int irank=0; irank<baseMPI->mpisize; irank++) {
    if (baseMPI->mpirank == 0) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
    int n2 = nj;
    MPI_Shift(x2, nj, baseMPI->mpisize, baseMPI->mpirank);
    nj = n2;
    MPI_Shift(y2, nj, baseMPI->mpisize, baseMPI->mpirank);
    nj = n2;
    MPI_Shift(z2, nj, baseMPI->mpisize, baseMPI->mpirank);
    nj = n2;
    MPI_Shift(v2, nj, baseMPI->mpisize, baseMPI->mpirank);
    for (int i=0; i<ni; i++) {
      double pp = 0;
      for (int j=0; j<nj; j++) {
	double dx = xi[i] - x2[j];
	double dy = yi[i] - y2[j];
	double dz = zi[i] - z2[j];
	double R2 = dx * dx + dy * dy + dz * dz;
	double invR = 1 / std::sqrt(R2);
	if (R2 == 0) invR = 0;
	pp += v2[j] * invR;
      }
      vi[i] += pp;
    }
  }
  delete[] x2;
  delete[] y2;
  delete[] z2;
  delete[] v2;
}
