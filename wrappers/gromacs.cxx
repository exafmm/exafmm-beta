#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree_from_cluster.h"
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
BuildTreeFromCluster * clusterTree;
BuildTree * localTree, * globalTree;
Partition * partition;
Traversal * traversal;
TreeMPI * treeMPI;
UpDownPass * upDownPass;

Bodies buffer;
Bounds localBounds;
Bounds globalBounds;

extern "C" void FMM_Init(int images, int threads, bool verbose) {
  const int ncrit = 32;
  const int nspawn = 1000;
  const real_t theta = 0.5;
  const bool useRmax = false;
  const bool useRopt = false;
  kernel::eps2 = 0.0;
  kernel::setup();

  args = new Args;
  baseMPI = new BaseMPI;
  boundBox = new BoundBox(nspawn);
  clusterTree = new BuildTreeFromCluster();
  localTree = new BuildTree(ncrit, nspawn);
  globalTree = new BuildTree(1, nspawn);
  partition = new Partition(baseMPI->mpirank, baseMPI->mpisize);
  traversal = new Traversal(nspawn, images);
  treeMPI = new TreeMPI(baseMPI->mpirank, baseMPI->mpisize, images);
  upDownPass = new UpDownPass(theta, useRmax, useRopt);

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
  args->threads = threads;
  args->verbose = verbose & (baseMPI->mpirank == 0);
  args->useRmax = useRmax;
  logger::verbose = args->verbose;
  logger::printTitle("Initial Parameters");
  args->print(logger::stringLength, P);
}

extern "C" void FMM_Finalize() {
  delete args;
  delete baseMPI;
  delete boundBox;
  delete clusterTree;
  delete localTree;
  delete globalTree;
  delete partition;
  delete traversal;
  delete treeMPI;
  delete upDownPass;
}

extern "C" void FMM_Partition(int & n, int * ibody, int * icell, float * x, float * q, float cycle) {
  logger::printTitle("Partition Profiling");
  num_threads(args->threads);
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0] - cycle / 2;
    B->X[1] = x[3*i+1] - cycle / 2;
    B->X[2] = x[3*i+2] - cycle / 2;
    B->SRC = q[i];
    int iwrap = wrap(B->X, cycle);
    B->IBODY = ibody[i] | (iwrap << shift);
    B->ICELL = icell[i];
  }
  localBounds = boundBox->getBounds(bodies);
  globalBounds = baseMPI->allreduceBounds(localBounds);
  localBounds = partition->octsection(bodies,globalBounds);
  bodies = treeMPI->commBodies(bodies);
#if Cluster
  Bodies clusters = clusterTree->setClusterCenter(bodies, cycle);
  Cells cells = globalTree->buildTree(clusters, buffer, localBounds);
  clusterTree->attachClusterBodies(bodies, cells, cycle);
#else
  Cells cells = localTree->buildTree(bodies, buffer, localBounds);
#endif
  upDownPass->upwardPass(cells);

#if Cluster
  clusterTree->shiftBackBodies(bodies, cycle);
#endif
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    ibody[i] = B->IBODY & mask;
    icell[i] = B->ICELL;
    int iwrap = unsigned(B->IBODY) >> shift;
    unwrap(B->X, cycle, iwrap);
    x[3*i+0] = B->X[0] + cycle / 2;
    x[3*i+1] = B->X[1] + cycle / 2;
    x[3*i+2] = B->X[2] + cycle / 2;
    q[i]     = B->SRC;
  }
  n = bodies.size();
}

extern "C" void FMM_Coulomb(int n, int * index, float * x, float * q, float * p, float * f, float cycle) {
  num_threads(args->threads);
  args->numBodies = n;
  logger::printTitle("FMM Parameters");
  args->print(logger::stringLength, P);
  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  logger::startPAPI();
  Bodies bodies(n);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = x[3*i+0] - cycle / 2;
    B->X[1] = x[3*i+1] - cycle / 2;
    B->X[2] = x[3*i+2] - cycle / 2;
    wrap(B->X, cycle);
    B->SRC = q[i];
    B->TRG[0] = p[i];
    B->TRG[1] = f[3*i+0];
    B->TRG[2] = f[3*i+1];
    B->TRG[3] = f[3*i+2];
    B->IBODY = i;
    B->ICELL = index[i];
  }
#if Cluster
  Bodies clusters = clusterTree->setClusterCenter(bodies, cycle);
  Cells cells = globalTree->buildTree(clusters, buffer, localBounds);
  clusterTree->attachClusterBodies(bodies, cells, cycle);
#else
  Cells cells = localTree->buildTree(bodies, buffer, localBounds);
#endif
  upDownPass->upwardPass(cells);
  treeMPI->allgatherBounds(localBounds);
  treeMPI->setLET(cells, cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->initListCount(cells);
  traversal->initWeight(cells);
  traversal->traverse(cells, cells, cycle, args->dual, args->mutual);
#if COUNT_LIST
  traversal->writeList(cells, baseMPI->mpirank);
#endif
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
  vec3 localDipole = upDownPass->getDipole(bodies,0);
  vec3 globalDipole = baseMPI->allreduceVec3(localDipole);
  int numBodies = baseMPI->allreduceInt(bodies.size());
  upDownPass->dipoleCorrection(bodies, globalDipole, numBodies, cycle);
  logger::stopPAPI();
  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");
#if Cluster
  clusterTree->shiftBackBodies(bodies, cycle);
#endif
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B->IBODY;
    p[i]     = B->TRG[0];
    f[3*i+0] = B->TRG[1];
    f[3*i+1] = B->TRG[2];
    f[3*i+2] = B->TRG[3];
  }
}

extern "C" void Ewald_Coulomb(int n, float * x, float * q, float * p, float * f,
			      int ksize, float alpha, float sigma, float cutoff, float cycle) {
  num_threads(args->threads);
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

void MPI_Shift(float * var, int &nold, int mpisize, int mpirank) {
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  int nnew;
  MPI_Request sreq, rreq;
  MPI_Isend(&nold, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(&nnew, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  float * buf = new float [nnew];
  MPI_Isend(var, nold, MPI_FLOAT, irecv, 1, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(buf, nnew, MPI_FLOAT, isend, 1, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  for (int i=0; i<nnew; i++) {
    var[i] = buf[i];
  }
  nold = nnew;
  delete[] buf;
}

extern "C" void Direct_Coulomb(int Ni, float * x, float * q, float * p, float * f, float cycle) {
  const int Nmax = 1000000;
  int images = args->images;
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  float * x2 = new float [3*Nmax];
  float * q2 = new float [Nmax];
  for (int i=0; i<Ni; i++) {
    x2[3*i+0] = x[3*i+0];
    x2[3*i+1] = x[3*i+1];
    x2[3*i+2] = x[3*i+2];
    q2[i] = q[i];
  }
  float Xperiodic[3];
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
  float localDipole[3] = {0, 0, 0};
  for (int i=0; i<Ni; i++) {
    for (int d=0; d<3; d++) localDipole[d] += x[3*i+d] * q[i];
  }
  int N;
  MPI_Allreduce(&Ni, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  float globalDipole[3];
  MPI_Allreduce(localDipole, globalDipole, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  double norm = 0;
  for (int d=0; d<3; d++) {
    norm += globalDipole[d] * globalDipole[d];
  }
  float coef = 4 * M_PI / (3 * cycle * cycle * cycle);
  for (int i=0; i<Ni; i++) {
    p[i] -= coef * norm / N / q[i];
    f[3*i+0] -= coef * globalDipole[0];
    f[3*i+1] -= coef * globalDipole[1];
    f[3*i+2] -= coef * globalDipole[2];
  }
  delete[] x2;
  delete[] q2;
}
