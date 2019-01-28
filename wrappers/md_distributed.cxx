#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree_from_cluster.h"
#include "ewald.h"
#include "kernel.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"

namespace exafmm_laplace {
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);

  Args * args;
  BaseMPI * baseMPI;
  BoundBox * boundBox;
  BuildTreeFromCluster * clusterTree;
  BuildTree * localTree, * globalTree;
  Kernel * kernel;
  Partition * partition;
  Traversal * traversal;
  TreeMPI * treeMPI;
  UpDownPass * upDownPass;
  Verify * verify;

  bool isTime;
  bool pass;
  Bodies buffer;
  Bounds localBounds;
  Bounds globalBounds;

  extern "C" void FMM_Init(int images, int threads, double theta, double cutoff, bool verbose, const char * path) {
    const int P = 10;
    const int ncrit = 32;
    const int nspawn = 1000;

    args = new Args;
    baseMPI = new BaseMPI;
    boundBox = new BoundBox;
    clusterTree = new BuildTreeFromCluster();
    kernel = new Kernel(P, eps2, wavek);
    localTree = new BuildTree(ncrit);
    globalTree = new BuildTree(1);
    partition = new Partition(*baseMPI);
    traversal = new Traversal(*kernel, theta, nspawn, images, path);
    treeMPI = new TreeMPI(*kernel, *baseMPI, theta, images);
    upDownPass = new UpDownPass(*kernel);
    verify = new Verify(path);

    args->accuracy = 1;
    args->ncrit = ncrit;
    args->cutoff = cutoff;
    args->distribution = "external";
    args->dual = 1;
    args->graft = 1;
    args->images = images;
    args->numBodies = 0;
    args->nspawn = nspawn;
    args->path = path;
    args->theta = theta;
    args->threads = threads;
    args->verbose = verbose & (baseMPI->mpirank == 0);
    verify->verbose = args->verbose;
    logger::verbose = args->verbose;
    logger::path = args->path;
    logger::printTitle("Initial Parameters");
    args->print(logger::stringLength);
    pass = true;
    isTime = false;
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

  extern "C" void Set_Index(int * ni, int nimax, int * res_index, double * x, double * q, double * v, double * cycle) {
    num_threads(args->threads);
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    const int shift = 29;
    const int mask = ~(0x7U << shift);
    Bodies bodies(*ni);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B-bodies.begin();
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG[0] = v[3*i+0];
      B->TRG[1] = v[3*i+1];
      B->TRG[2] = v[3*i+2];
      int iwrap = wrap(B->X, cycles);
      B->IBODY = i | (iwrap << shift);
      B->ICELL = res_index[i];
    }
    localBounds = boundBox->getBounds(bodies);
    Cells cells = localTree->buildTree(bodies, buffer, localBounds);
    upDownPass->upwardPass(cells);
    int id = 0;
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {
      int ic = 0;
      for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++) {
        if ((drand48() > 0.5 && ic > 2) || ic > 5) ic = 0;
        if (ic == 0) {
          B->ICELL = id;
          id++;
        } else {
          B->ICELL = -ic;
        }
        ic++;
      }
    }
    *ni = bodies.size();
    if (*ni < nimax) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
        int i = B-bodies.begin();
        res_index[i] = B->ICELL;
        int iwrap = unsigned(B->IBODY) >> shift;
        unwrap(B->X, cycles, iwrap);
        x[3*i+0] = B->X[0];
        x[3*i+1] = B->X[1];
        x[3*i+2] = B->X[2];
        q[i]     = B->SRC;
        v[3*i+0] = B->TRG[0];
        v[3*i+1] = B->TRG[1];
        v[3*i+2] = B->TRG[2];
      }
    }
  }

  extern "C" void FMM_Partition(int * ni, int nimax, int * res_index, double * x, double * q, double * v, double * cycle) {
    num_threads(args->threads);
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    logger::printTitle("Partition Profiling");
    const int shift = 29;
    const int mask = ~(0x7U << shift);
    Bodies bodies(*ni);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B-bodies.begin();
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      B->SRC = q[i];
      B->TRG[0] = v[3*i+0];
      B->TRG[1] = v[3*i+1];
      B->TRG[2] = v[3*i+2];
      int iwrap = wrap(B->X, cycles);
      B->IBODY = i | (iwrap << shift);
      B->ICELL = res_index[i];
    }
    localBounds = boundBox->getBounds(bodies);
    globalBounds = baseMPI->allreduceBounds(localBounds);
    localBounds = partition->octsection(bodies,globalBounds);
    bodies = treeMPI->commBodies(bodies);
    Cells cells = localTree->buildTree(bodies, buffer, localBounds);
    upDownPass->upwardPass(cells);

    *ni = bodies.size();
    if (*ni < nimax) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
        int i = B-bodies.begin();
        res_index[i] = B->ICELL;
        int iwrap = unsigned(B->IBODY) >> shift;
        unwrap(B->X, cycles, iwrap);
        x[3*i+0] = B->X[0];
        x[3*i+1] = B->X[1];
        x[3*i+2] = B->X[2];
        q[i]     = B->SRC;
        v[3*i+0] = B->TRG[0];
        v[3*i+1] = B->TRG[1];
        v[3*i+2] = B->TRG[2];
      }
    }
  }

  extern "C" void FMM_FMM(int ni, int * nj, int * res_index, double * x, double * q, double * p, double * f, double * cycle) {
    num_threads(args->threads);
    double cutoff = args->cutoff;
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    const int shift = 29;
    const int mask = ~(0x7U << shift);
    args->numBodies = ni;
    logger::printTitle("FMM Parameters");
    args->print(logger::stringLength);
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    Bodies bodies(ni);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B-bodies.begin();
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      wrap(B->X, cycles);
      B->SRC = q[i];
      B->TRG[0] = p[i];
      B->TRG[1] = f[3*i+0];
      B->TRG[2] = f[3*i+1];
      B->TRG[3] = f[3*i+2];
      int iwrap = wrap(B->X, cycles);
      B->IBODY = i | (iwrap << shift);
      B->ICELL = res_index[i];
    }
    Cells cells = localTree->buildTree(bodies, buffer, localBounds);
    upDownPass->upwardPass(cells);
    treeMPI->allgatherBounds(localBounds);
    treeMPI->setLET(cells, cycles);
    treeMPI->commBodies();
    treeMPI->commCells();
    traversal->initListCount(cells);
    traversal->initWeight(cells);
    traversal->traverse(cells, cells, cycles, args->dual);
#if EXAFMM_COUNT_LIST
    traversal->writeList(cells, baseMPI->mpirank);
#endif
    Cells jcells;
    if (baseMPI->mpisize > 1) {
      if (args->graft) {
        treeMPI->linkLET();
        Bodies gbodies = treeMPI->root2body();
        jcells = globalTree->buildTree(gbodies, buffer, globalBounds);
        treeMPI->attachRoot(jcells);
        traversal->traverse(cells, jcells, cycles, args->dual);
      } else {
        for (int irank=0; irank<baseMPI->mpisize; irank++) {
          treeMPI->getLET(jcells, (baseMPI->mpirank+irank)%baseMPI->mpisize);
          traversal->traverse(cells, jcells, cycles, args->dual);
        }
      }
    }
    upDownPass->downwardPass(cells);
    vec3 localDipole = upDownPass->getDipole(bodies,0);
    vec3 globalDipole = baseMPI->allreduceVec3(localDipole);
    int numBodies = baseMPI->allreduceInt(bodies.size());
    upDownPass->dipoleCorrection(bodies, globalDipole, numBodies, cycles);
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
    logger::printTitle("Total runtime");
    logger::printTime("Total FMM");
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B - bodies.begin();
      p[i]     = B->TRG[0];
      f[3*i+0] = B->TRG[1];
      f[3*i+1] = B->TRG[2];
      f[3*i+2] = B->TRG[3];
    }
    Bodies jbodies = treeMPI->getRecvBodies();
    jbodies.insert(jbodies.begin(), bodies.begin(), bodies.end());
    bodies.clear();
    vec3 Xmin = localBounds.Xmin - cutoff;
    vec3 Xmax = localBounds.Xmax + cutoff;
    ivec3 iX;
    real_t X[3];
    for (iX[0]=-1; iX[0]<=1; iX[0]++) {
      for (iX[1]=-1; iX[1]<=1; iX[1]++) {
        for (iX[2]=-1; iX[2]<=1; iX[2]++) {
          if (norm(iX) != 0) {
            for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
              for (int d=0; d<3; d++) X[d] = B->X[d] + iX[d] * cycles[d];
              if (Xmin[0] < X[0] && X[0] < Xmax[0] &&
                  Xmin[1] < X[1] && X[1] < Xmax[1] &&
                  Xmin[2] < X[2] && X[2] < Xmax[2]) {
                bodies.push_back(*B);
                for (int d=0; d<3; d++) {
                  bodies.back().X[d] = X[d];
                }
              }
            }
          }
        }
      }
    }
    jbodies = treeMPI->getRecvBodies();
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      if (Xmin[0] < B->X[0] && B->X[0] < Xmax[0] &&
          Xmin[1] < B->X[1] && B->X[1] < Xmax[1] &&
          Xmin[2] < B->X[2] && B->X[2] < Xmax[2]) {
        bodies.push_back(*B);
      }
    }
    int njmax = *nj;
    *nj = ni + bodies.size();
    if (*nj < njmax) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
        int i = B - bodies.begin() + ni;
        res_index[i] = B->ICELL;
        x[3*i+0] = B->X[0];
        x[3*i+1] = B->X[1];
        x[3*i+2] = B->X[2];
        q[i] = B->SRC;
      }
    }
  }

  extern "C" void FMM_Ewald(int ni, double * x, double * q, double * p, double * f,
                            int ksize, double alpha, double sigma, double cutoff, double * cycle) {
    num_threads(args->threads);
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    Ewald * ewald = new Ewald(ksize, alpha, sigma, cutoff, cycles);
    args->numBodies = ni;
    logger::printTitle("Ewald Parameters");
    args->print(logger::stringLength);
    ewald->print(logger::stringLength);
    logger::printTitle("Ewald Profiling");
    logger::startTimer("Total Ewald");
    logger::startPAPI();
    Bodies bodies(ni);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B-bodies.begin();
      B->X[0] = x[3*i+0];
      B->X[1] = x[3*i+1];
      B->X[2] = x[3*i+2];
      wrap(B->X, cycles);
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
      Bounds jlocalBounds = boundBox->getBounds(jbodies);
      Cells jcells = localTree->buildTree(jbodies, buffer, jlocalBounds);
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
    MPI_Isend(&nold, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);
    MPI_Irecv(&nnew, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);
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

  extern "C" void Dipole_Correction(int ni, double * x, double * q, double * p, double * f, double * cycle) {
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    float localDipole[3] = {0, 0, 0};
    for (int i=0; i<ni; i++) {
      for (int d=0; d<3; d++) localDipole[d] += x[3*i+d] * q[i];
    }
    int N;
    MPI_Allreduce(&ni, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    float globalDipole[3];
    MPI_Allreduce(localDipole, globalDipole, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    double norm = 0;
    for (int d=0; d<3; d++) {
      norm += globalDipole[d] * globalDipole[d];
    }
    float coef = 4 * M_PI / (3 * cycles[0] * cycles[1] * cycles[2]);
    for (int i=0; i<ni; i++) {
      p[i] -= coef * norm / N / q[i];
      f[3*i+0] -= coef * globalDipole[0];
      f[3*i+1] -= coef * globalDipole[1];
      f[3*i+2] -= coef * globalDipole[2];
    }
  }

  extern "C" void FMM_Cutoff(int ni, double * x, double * q, double * p, double * f, double cutoff, double * cycle) {
    vec3 cycles;
    for (int d=0; d<3; d++) cycles[d] = cycle[d];
    const int Nmax = 1000000;
    const double cutoff2 = cutoff * cutoff;
    int prange = int(cutoff/min(cycles)*0.999999) + 1;
    double * x2 = new double [3*Nmax];
    double * q2 = new double [Nmax];
    for (int i=0; i<ni; i++) {
      x2[3*i+0] = x[3*i+0];
      x2[3*i+1] = x[3*i+1];
      x2[3*i+2] = x[3*i+2];
      q2[i] = q[i];
    }
    float Xperiodic[3];
    int nj = ni, nj3 = 3 * ni;
    if (baseMPI->mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
    for (int irank=0; irank<baseMPI->mpisize; irank++) {
      if (baseMPI->mpirank == 0) std::cout << "Direct loop          : " << irank+1 << "/" << baseMPI->mpisize << std::endl;
      MPI_Shift(x2, nj3, baseMPI->mpisize, baseMPI->mpirank);
      MPI_Shift(q2, nj,  baseMPI->mpisize, baseMPI->mpirank);
      for (int i=0; i<ni; i++) {
        double pp = 0, fx = 0, fy = 0, fz = 0;
        for (int ix=-prange; ix<=prange; ix++) {
          for (int iy=-prange; iy<=prange; iy++) {
            for (int iz=-prange; iz<=prange; iz++) {
              Xperiodic[0] = ix * cycles[0];
              Xperiodic[1] = iy * cycles[1];
              Xperiodic[2] = iz * cycles[2];
              for (int j=0; j<nj; j++) {
                double dx = x[3*i+0] - x2[3*j+0] - Xperiodic[0];
                double dy = x[3*i+1] - x2[3*j+1] - Xperiodic[1];
                double dz = x[3*i+2] - x2[3*j+2] - Xperiodic[2];
                double R2 = dx * dx + dy * dy + dz * dz;
                if (R2 < cutoff2) {
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
        }
        p[i] += pp;
        f[3*i+0] -= fx;
        f[3*i+1] -= fy;
        f[3*i+2] -= fz;
      }
    }
    delete[] x2;
    delete[] q2;
  }

  extern "C" void FMM_Verify_Accuracy(int &t, double potRel, double accRel) {
    isTime = false;
    logger::printTitle("Accuracy regression");
    if (!baseMPI->mpirank) {
      pass = verify->regression(args->getKey(baseMPI->mpisize), isTime, t, potRel, accRel);
    }
    MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pass) {
      if (verify->verbose) std::cout << "passed accuracy regression at t: " << t << std::endl;
      t = -1;
    }
  }

  extern "C" bool FMM_Only_Accuracy() {
    return args->accuracy == 1;
  }

  extern "C" void FMM_Verify_Time(int &t, double totalFMM) {
    isTime = true;
    logger::printTitle("Time regression");
    double totalFMMGlob;
    MPI_Reduce(&totalFMM, &totalFMMGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    totalFMMGlob /= baseMPI->mpisize;
    if (!baseMPI->mpirank) {
      pass = verify->regression(args->getKey(baseMPI->mpisize), isTime, t, totalFMMGlob);
    }
    MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pass) {
      if (verify->verbose) std::cout << "passed time regression at t: " << t << std::endl;
      t = -1;
    }
  }

  extern "C" void FMM_Verify_End() {
    if (!pass) {
      if (verify->verbose) {
        if(!isTime) std::cout << "failed accuracy regression" << std::endl;
        else std::cout << "failed time regression" << std::endl;
      }
      abort();
    }
  }
}
