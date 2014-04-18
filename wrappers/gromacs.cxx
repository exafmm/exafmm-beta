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
#if MASS
#error Turn off MASS for this wrapper
#endif

Args * args;
BaseMPI * baseMPI;
BoundBox * boundBox;
BuildTree * localTree, * globalTree;
Partition * partition;
Traversal * traversal;
TreeMPI * treeMPI;
UpDownPass * upDownPass;

Bounds localBounds;
Bounds globalBounds;

extern "C" void FMM_Init(int images) {
  const int ncrit = 32;
  const int nspawn = 1000;
  const real_t theta = 0.4;
  const bool useRmax = true;
  const bool useRopt = true;
  args = new Args;
  baseMPI = new BaseMPI;
  boundBox = new BoundBox(nspawn);
  localTree = new BuildTree(ncrit, nspawn);
  globalTree = new BuildTree(1, nspawn);
  partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
  traversal = new Traversal(nspawn, images);
  treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
  upDownPass = new UpDownPass(theta, useRmax, useRopt);

  args->theta = theta;
  args->ncrit = ncrit;
  args->nspawn = nspawn;
  args->images = images;
  args->mutual = 0;
  args->verbose = 1;
  args->distribution = "external";
  args->verbose &= baseMPI->mpirank == 0;
  logger::verbose = args->verbose;
  logger::printTitle("Initial Parameters");
  args->print(logger::stringLength, P);
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

extern "C" void FMM_Partition(int & n, int * index, double * x, double * q, double cycle) {
  logger::printTitle("Partition Profiling");
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    B->SRC = q[i];
    int iwrap = wrap(B->X, cycle);
    B->IBODY = index[i] | (iwrap << shift);
  }
  localBounds = boundBox->getBounds(bodies);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  localBounds = partition->octsection(bodies,globalBounds);
  bodies = treeMPI->commBodies(bodies);
  Cells cells = localTree->buildTree(bodies, localBounds);
  upDownPass->upwardPass(cells);

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    index[i] = B->IBODY & mask;
    int iwrap = unsigned(B->IBODY) >> shift;
    unwrap(B->X, cycle, iwrap);
    x[3*i+0] = B->X[0];
    x[3*i+1] = B->X[1];
    x[3*i+2] = B->X[2];
    q[i]     = B->SRC;
  }
  n = bodies.size();
}

extern "C" void FMM_Coulomb(int n, double * x, double * q, double * p, double * f, double cycle) {
  args->numBodies = n;
  logger::printTitle("FMM Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startPAPI();
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    wrap(B->X, cycle);
    B->SRC = q[i];
    B->TRG[0] = p[i];
    B->TRG[1] = f[3*i+0];
    B->TRG[2] = f[3*i+1];
    B->TRG[3] = f[3*i+2];
    B->IBODY = i;
  }
  Cells cells = localTree->buildTree(bodies, localBounds);
  upDownPass->upwardPass(cells);
  treeMPI->allgatherBounds(localBounds);
  treeMPI->setLET(cells, cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->initWeight(cells);
  traversal->dualTreeTraversal(cells, cells, cycle, args->mutual);
  Cells jcells;
  if (args->graft) {
    treeMPI->linkLET();
    Bodies gbodies = treeMPI->root2body();
    jcells = globalTree->buildTree(gbodies, globalBounds);
    treeMPI->attachRoot(jcells);
    traversal->dualTreeTraversal(cells, jcells, cycle, false);
  } else {
    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      treeMPI->getLET(jcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
      traversal->dualTreeTraversal(cells, jcells, cycle, false);
    }
  }
  upDownPass->downwardPass(cells);
  vec3 localDipole = upDownPass->getDipole(bodies,0);
  vec3 globalDipole = baseMPI->allreduceVec3(localDipole);
  int numBodies = baseMPI->allreduceInt(bodies.size());
  upDownPass->dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger::stopPAPI();
  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
}

extern "C" void Ewald_Coulomb(int n, double * x, double * q, double * p, double * f,
			      int ksize, double alpha, double sigma, double cutoff, double cycle) {
  Ewald * ewald = new Ewald(ksize, alpha, sigma, cutoff, cycle);
  args->numBodies = n;
  logger::printTitle("Ewald Parameters");
  args->print(logger::stringLength, P);
  ewald->print(logger::stringLength);
  logger::printTitle("Ewald Profiling");
  logger::startTimer("Total Ewald");
  logger::startPAPI();
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0];
    B->X[1] = x[3*i+1];
    B->X[2] = x[3*i+2];
    wrap(B->X, cycle);
    B->SRC = q[i];
    B->TRG[0] = p[i];
    B->TRG[1] = f[3*i+0];
    B->TRG[2] = f[3*i+1];
    B->TRG[3] = f[3*i+2];
    B->IBODY = i;
  }
  Cells cells = localTree->buildTree(bodies, localBounds);
  Bodies jbodies = bodies;
  for (int i=0; i<baseMPI->mpisize; i++) {
    if (args->verbose) std::cout << "Ewald loop           : " << i+1 << "/" << baseMPI->mpisize << std::endl;
    treeMPI->shiftBodies(jbodies);
    localBounds = boundBox->getBounds(jbodies);
    Cells jcells = localTree->buildTree(jbodies, localBounds);
    ewald->wavePart(bodies, jbodies);
    ewald->realPart(cells, jcells);
  }
  ewald->selfTerm(bodies);
  logger::stopPAPI();
  logger::stopTimer("Total Ewald");
  logger::printTitle("Total runtime");
  logger::printTime("Total Ewald");
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
  delete ewald;
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

extern "C" void Direct_Coulomb(int Ni, double * x, double * q, double * p, double * f, double cycle) {
  const int Nmax = 1000000;
  int images = args->images;
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  double * x2 = new double [3*Nmax];
  double * q2 = new double [Nmax];
  for (int i=0; i<Ni; i++) {
    x2[3*i+0] = x[3*i+0];
    x2[3*i+1] = x[3*i+1];
    x2[3*i+2] = x[3*i+2];
    q2[i] = q[i];
  }
  double Xperiodic[3];
  int Nj = Ni, Nj3 = 3 * Ni;
  if (baseMPI->mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
  for (int irank=0; irank<baseMPI->mpisize; irank++) {
    if (baseMPI->mpirank == 0) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
    MPI_Shift(x2, Nj3, baseMPI->mpisize, baseMPI->mpirank);
    MPI_Shift(q2, Nj,  baseMPI->mpisize, baseMPI->mpirank);
    for (int i=0; i<Ni; i++) {
      double pp = 0, fx = 0, fy = 0, fz = 0;
      for (int ix=-prange; ix<=prange; ix++) {
        for (int iy=-prange; iy<=prange; iy++) {
          for (int iz=-prange; iz<=prange; iz++) {
            Xperiodic[0] = ix * cycle;
            Xperiodic[1] = iy * cycle;
            Xperiodic[2] = iz * cycle;
            for (int j=0; j<Nj; j++) {
              double dx = x[3*i+0] - x2[3*j+0] - Xperiodic[0];
              double dy = x[3*i+1] - x2[3*j+1] - Xperiodic[1];
              double dz = x[3*i+2] - x2[3*j+2] - Xperiodic[2];
              double R2 = dx * dx + dy * dy + dz * dz;
              double invR = 1 / std::sqrt(R2);
              if (R2 == 0) invR = 0;
              double invR3 = q2[j] * invR * invR * invR;
              pp += q2[j] * invR;
              fx += dx * invR3;
              fy += dy * invR3;
              fz += dz * invR3;
            }
          }
        }
      }
      p[i] += pp;
      f[3*i+0] -= fx;
      f[3*i+1] -= fy;
      f[3*i+2] -= fz;
    }
  }
  double localDipole[3] = {0, 0, 0};
  for (int i=0; i<Ni; i++) {
    for (int d=0; d<3; d++) localDipole[d] += x[3*i+d] * q[i];
  }
  int N;
  MPI_Allreduce(&Ni, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  double globalDipole[3];
  MPI_Allreduce(localDipole, globalDipole, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = 0;
  for (int d=0; d<3; d++) {
    norm += globalDipole[d] * globalDipole[d];
  }
  double coef = 4 * M_PI / (3 * cycle * cycle * cycle);
  for (int i=0; i<Ni; i++) {
    p[i] -= coef * norm / N / q[i];
    f[3*i+0] -= coef * globalDipole[0];
    f[3*i+1] -= coef * globalDipole[1];
    f[3*i+2] -= coef * globalDipole[2];
  }
  delete[] x2;
  delete[] q2;
}
