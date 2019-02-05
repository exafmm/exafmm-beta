#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "ewald.h"
#include "kernel.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "van_der_waals.h"
#include "verify.h"

namespace exafmm_laplace {
  static const double Celec = 332.0716;
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);

  Args * args;
  BaseMPI * baseMPI;
  BoundBox * boundBox;
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

  extern "C" void fmm_init_(int & expansions, int & images, double & theta, int & verbose, int &, char * path, size_t *) {
    const int P = expansions;
    const int ncrit = 16;
    const int nspawn = 1000;

    args = new Args;
    baseMPI = new BaseMPI;
    boundBox = new BoundBox;
    localTree = new BuildTree(ncrit);
    globalTree = new BuildTree(1);
    kernel = new Kernel(P, eps2, wavek);
    partition = new Partition(*baseMPI);
    traversal = new Traversal(*kernel, theta, nspawn, images, path);
    treeMPI = new TreeMPI(*kernel, *baseMPI, theta, images);
    upDownPass = new UpDownPass(*kernel);
    verify = new Verify(path);

    args->ncrit = ncrit;
    args->accuracy = 1;
    args->distribution = "external";
    args->dual = 1;
    args->graft = 0;
    args->images = images;
    args->numBodies = 0;
    args->path = path;
    args->nspawn = nspawn;
    args->theta = theta;
    args->verbose = verbose & (baseMPI->mpirank == 0);
    verify->verbose = baseMPI->mpirank == 0;
    logger::verbose = args->verbose;
    logger::path = args->path;
    logger::printTitle("Initial Parameters");
    args->print(logger::stringLength);
    pass = true;
    isTime = false;
  }

  extern "C" void fmm_finalize_() {
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

  extern "C" void fmm_partition_(int & nglobal, int * icpumap, double * x, double * q,
                                 double * xold, double & cycle) {
    vec3 cycles = cycle;
    logger::printTitle("Partition Profiling");
    const int shift = 29;
    const int mask = ~(0x7U << shift);
    int nlocal = 0;
    for (int i=0; i<nglobal; i++) {
      if (icpumap[i] == 1) nlocal++;
    }
    Bodies bodies(nlocal);
    B_iter B = bodies.begin();
    for (int i=0; i<nglobal; i++) {
      if (icpumap[i] == 1) {
        B->X[0] = x[3*i+0];
        B->X[1] = x[3*i+1];
        B->X[2] = x[3*i+2];
        B->SRC = q[i];
        B->TRG[0] = xold[3*i+0];
        B->TRG[1] = xold[3*i+1];
        B->TRG[2] = xold[3*i+2];
        int iwrap = wrap(B->X, cycles);
        B->IBODY = i | (iwrap << shift);
        B++;
      }
    }
    localBounds = boundBox->getBounds(bodies);
    globalBounds = baseMPI->allreduceBounds(localBounds);
    localBounds = partition->octsection(bodies,globalBounds);
    bodies = treeMPI->commBodies(bodies);
    for (int i=0; i<nglobal; i++) {
      icpumap[i] = 0;
    }
    for (B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B->IBODY & mask;
      int iwrap = unsigned(B->IBODY) >> shift;
      unwrap(B->X, cycles, iwrap);
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
    vec3 cycles = cycle;
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
    logger::printTitle("FMM Parameters");
    args->print(logger::stringLength);
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    Bodies bodies(nlocal);
    B_iter B = bodies.begin();
    for (int i=0; i<nglobal; i++) {
      if (icpumap[i] == 1) {
        B->X[0] = x[3*i+0];
        B->X[1] = x[3*i+1];
        B->X[2] = x[3*i+2];
        B->SRC = q[i] == 0 ? EPS : q[i];
        B->TRG = 0;
        int iwrap = wrap(B->X, cycles);
        B->IBODY = i | (iwrap << shift);
        B++;
      }
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

    for (B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B->IBODY & mask;
      p[i]     += B->TRG[0] * B->SRC * Celec;
      f[3*i+0] += B->TRG[1] * B->SRC * Celec;
      f[3*i+1] += B->TRG[2] * B->SRC * Celec;
      f[3*i+2] += B->TRG[3] * B->SRC * Celec;
    }
    bodies = treeMPI->getRecvBodies();
    for (B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B->IBODY & mask;
      int iwrap = unsigned(B->IBODY) >> shift;
      unwrap(B->X, cycles, iwrap);
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
    vec3 cycles = cycle;
    Ewald * ewald = new Ewald(ksize, alpha, sigma, cutoff, cycles);
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
    args->print(logger::stringLength);
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
        int iwrap = wrap(B->X, cycles);
        B->IBODY = i | (iwrap << shift);
        B++;
      }
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

    for (B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B->IBODY & mask;
      p[i]     += B->TRG[0] * B->SRC * Celec;
      f[3*i+0] += B->TRG[1] * B->SRC * Celec;
      f[3*i+1] += B->TRG[2] * B->SRC * Celec;
      f[3*i+2] += B->TRG[3] * B->SRC * Celec;
    }
    treeMPI->allgatherBounds(localBounds);
    treeMPI->setLET(cells, cycles);
    treeMPI->commBodies();
    bodies = treeMPI->getRecvBodies();
    for (B=bodies.begin(); B!=bodies.end(); B++) {
      int i = B->IBODY & mask;
      int iwrap = unsigned(B->IBODY) >> shift;
      unwrap(B->X, cycles, iwrap);
      x[3*i+0] = B->X[0];
      x[3*i+1] = B->X[1];
      x[3*i+2] = B->X[2];
      q[i] = B->SRC;
      assert(icpumap[i] == 0);
      icpumap[i] = 2;
    }
    delete ewald;
  }

  extern "C" void direct_coulomb_(int & nglobal, int * icpumap, double * x, double * q, double * p, double * f, double & cycle) {
    vec3 cycles = cycle;
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
              Xperiodic[0] = ix * cycles[0];
              Xperiodic[1] = iy * cycles[1];
              Xperiodic[2] = iz * cycles[2];
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
    real_t coef = 4 * M_PI / (3 * cycles[0] * cycles[1] * cycles[2]) * Celec;
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
    vec3 cycles = cycle;
    logger::startTimer("Coulomb Exclusion");
    for (int i=0, ic=0; i<nglobal; i++) {
      if (icpumap[i] == 1) {
        real_t pp = 0, fx = 0, fy = 0, fz = 0;
        for (int jc=0; jc<numex[i]; jc++, ic++) {
          int j = natex[ic]-1;
          vec3 dX;
          for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
          wrap(dX, cycles);
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
    vec3 cycles = cycle;
    VanDerWaals * VDW = new VanDerWaals(cuton, cutoff, cycles, numTypes, rscale, gscale, fgscale);
    const int shift = 29;
    const int mask = ~(0x7U << shift);
    int nlocal = 0;
    for (int i=0; i<nglobal; i++) {
      if (icpumap[i] == 1) nlocal++;
    }
    args->numBodies = nlocal;
    logger::printTitle("VdW Parameters");
    args->print(logger::stringLength);
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
        int iwrap = wrap(B->X, cycles);
        B->IBODY = i | (iwrap << shift);
        B++;
      }
    }

    Cells cells = localTree->buildTree(bodies, buffer, localBounds);
    upDownPass->upwardPass(cells);
    treeMPI->allgatherBounds(localBounds);
    treeMPI->setLET(cells, cycles);
    treeMPI->commBodies();
    treeMPI->commCells();
    VDW->evaluate(cells, cells);
    Cells jcells;
    for (int irank=1; irank<baseMPI->mpisize; irank++) {
      treeMPI->getLET(jcells,(baseMPI->mpirank+irank)%baseMPI->mpisize);
      VDW->evaluate(cells, jcells);
    }
    logger::stopPAPI();
    logger::stopTimer("Total VdW");
    logger::printTitle("Total runtime");
    logger::printTime("Total VdW");

    for (B=bodies.begin(); B!=bodies.end(); B++) {
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
    vec3 cycles = cycle;
    logger::startTimer("Direct VdW");
    for (int i=0; i<nglobal; i++) {
      if (icpumap[i] == 1) {
        int atypei = atype[i]-1;
        kreal_t pp = 0, fx = 0, fy = 0, fz = 0;
        for (int j=0; j<nglobal; j++) {
          vec3 dX;
          for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
          wrap(dX, cycles);
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
    vec3 cycles = cycle;
    logger::startTimer("VdW Exclusion");
    for (int i=0, ic=0; i<nglobal; i++) {
      if (icpumap[i] == 1) {
        int atypei = atype[i]-1;
        for (int jc=0; jc<numex[i]; jc++, ic++) {
          int j = natex[ic]-1;
          vec3 dX;
          for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
          wrap(dX, cycles);
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

  extern "C" void fmm_verify_accuracy_(int & t, double & potRel, double & accRel) {
    isTime = false;
    if (!baseMPI->mpirank) {
      std::cout << "key: " << args->getKey(baseMPI->mpisize) << std::endl;
      pass = verify->regression(args->getKey(baseMPI->mpisize), isTime, t-1, potRel, accRel);
    }
    MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pass) {
      if (verify->verbose) std::cout << "passed accuracy regression at t: " << t-1 << std::endl;
      t = -1;
    }
  }

  extern "C" void fmm_only_accuracy_(int & accuracy) {
    accuracy = args->accuracy;
  }

  extern "C" void fmm_verify_time_(int & t, double & totalFMM) {
    isTime = true;
    double totalFMMGlob;
    MPI_Reduce(&totalFMM, &totalFMMGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    totalFMMGlob /= baseMPI->mpisize;
    if (!baseMPI->mpirank) {
      pass = verify->regression(args->getKey(baseMPI->mpisize), isTime, t-1, totalFMMGlob);
    }
    MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    if (pass) {
      if (verify->verbose) std::cout << "passed time regression at t: " << t-1 << std::endl;
      t = -1;
    }
  }

  extern "C" void fmm_verify_end_() {
    if (!pass) {
      if (verify->verbose) {
        if(!isTime) std::cout << "failed accuracy regression" << std::endl;
        else std::cout << "failed time regression" << std::endl;
      }
      abort();
    }
  }
}
