#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "logger.h"
#include "partition.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
using namespace exafmm;
vec3 KernelBase::Xperiodic = 0;
real_t KernelBase::eps2 = 0.0;
complex_t KernelBase::wavek = complex_t(10.,1.) / real_t(2 * M_PI);

template<typename Kernel>
void fmm(Args args) {
  typedef typename Kernel::Bodies Bodies;                       //!< Vector of bodies
  typedef typename Kernel::Cells Cells;                         //!< Vector of cells
  typedef typename Kernel::B_iter B_iter;                       //!< Iterator of body vector
  typedef typename Kernel::C_iter C_iter;                       //!< Iterator of cell vector

  const vec3 cycle = 2 * M_PI;
  BaseMPI baseMPI;
  Bodies bodies, bodies2, jbodies, gbodies, buffer;
  BoundBox<Kernel> boundBox(args.nspawn);
  Bounds localBounds, globalBounds;
  BuildTree<Kernel> localTree(args.ncrit, args.nspawn);
  BuildTree<Kernel> globalTree(1, args.nspawn);
  Cells cells, jcells, gcells;
  Dataset<Kernel> data;
  Partition<Kernel> partition(baseMPI.mpirank, baseMPI.mpisize);
  TreeMPI<Kernel> treeMPI(baseMPI.mpirank, baseMPI.mpisize, args.images);
  Traversal<Kernel> traversal(args.nspawn, args.images, args.path);
  UpDownPass<Kernel> upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify<Kernel> verify(args.path);
  num_threads(args.threads);

  Kernel::init();
  //args.numBodies /= baseMPI.mpisize;
  args.verbose &= baseMPI.mpirank == 0;
  verify.verbose = args.verbose;
  logger::verbose = args.verbose;
  logger::path = args.path;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  buffer.reserve(bodies.size());
  if (args.IneJ) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    jbodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank+baseMPI.mpisize, baseMPI.mpisize);
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
  }
  if (args.mass) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->SRC = 1. / bodies.size();
    }
  }
  bool pass = true;
  bool isTime = false;
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    int numIteration = 1;
    if (isTime) numIteration = 10;
    for (int it=0; it<numIteration; it++) {
      std::stringstream title;
      title << "Time average loop " << it;
      logger::printTitle(title.str());
      data.initTarget(bodies);
      localBounds = boundBox.getBounds(bodies);
      if (args.IneJ) {
        localBounds = boundBox.getBounds(jbodies, localBounds);
      }
      globalBounds = baseMPI.allreduceBounds(localBounds);
      partition.bisection(bodies, globalBounds);
      bodies = treeMPI.commBodies(bodies);
      if (args.IneJ) {
        partition.bisection(jbodies, globalBounds);
        jbodies = treeMPI.commBodies(jbodies);
      }
      localBounds = boundBox.getBounds(bodies);
      cells = localTree.buildTree(bodies, buffer, localBounds);
      localBounds = boundBox.getBounds(cells, localBounds);
      upDownPass.upwardPass(cells);
      if (args.IneJ) {
        localBounds = boundBox.getBounds(jbodies);
        jcells = localTree.buildTree(jbodies, buffer, localBounds);
        localBounds = boundBox.getBounds(jcells, localBounds);
        upDownPass.upwardPass(jcells);
      }

#if 1 // Set to 0 for debugging by shifting bodies and reconstructing tree
      treeMPI.allgatherBounds(localBounds);
      if (args.IneJ) {
        treeMPI.setLET(jcells, cycle);
      } else {
        treeMPI.setLET(cells, cycle);
      }
#pragma omp parallel sections
      {
#pragma omp section
        {
          treeMPI.commBodies();
          treeMPI.commCells();
        }
#pragma omp section
        {
          traversal.initListCount(cells);
          traversal.initWeight(cells);
          if (args.IneJ) {
            traversal.traverse(cells, jcells, cycle, args.dual, false);
          } else {
            traversal.traverse(cells, cells, cycle, args.dual, args.mutual);
            jbodies = bodies;
          }
        }
      }
      if (baseMPI.mpisize > 1) {
        if (args.graft) {
          treeMPI.linkLET();
          gbodies = treeMPI.root2body();
          jcells = globalTree.buildTree(gbodies, buffer, globalBounds);
          treeMPI.attachRoot(jcells);
          traversal.traverse(cells, jcells, cycle, args.dual, false);
        } else {
          for (int irank=0; irank<baseMPI.mpisize; irank++) {
            treeMPI.getLET(jcells, (baseMPI.mpirank+irank)%baseMPI.mpisize);
            traversal.traverse(cells, jcells, cycle, args.dual, false);
          }
        }
      }
#else
      if (!args.IneJ) jbodies = bodies;
      for (int irank=0; irank<baseMPI.mpisize; irank++) {
        treeMPI.shiftBodies(jbodies);
        jcells.clear();
        localBounds = boundBox.getBounds(jbodies);
        jcells = localTree.buildTree(jbodies, buffer, localBounds);
        upDownPass.upwardPass(jcells);
        traversal.traverse(cells, jcells, cycle, args.dual, false);
      }
#endif
      upDownPass.downwardPass(cells);
    }
    logger::stopPAPI();
    double totalFMM = logger::stopTimer("Total FMM", 0);
    totalFMM /= numIteration;
    logger::resetTimer("Total FMM");
    if (args.write) {
      logger::writeTime(baseMPI.mpirank);
    }
    traversal.writeList(cells, baseMPI.mpirank);

    if (!isTime) {
      logger::printTitle("MPI direct sum");
      const int numTargets = 100;
      buffer = bodies;
      data.sampleBodies(bodies, numTargets);
      bodies2 = bodies;
      data.initTarget(bodies);
      logger::startTimer("Total Direct");
      for (int i=0; i<baseMPI.mpisize; i++) {
        if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << baseMPI.mpisize << std::endl;
        treeMPI.shiftBodies(jbodies);
        traversal.direct(bodies, jbodies, cycle);
      }
      traversal.normalize(bodies);
      logger::printTitle("Total runtime");
      logger::printTime("Total FMM");
      logger::stopTimer("Total Direct");
      logger::resetTimer("Total Direct");
      double potDif = verify.getDifScalar(bodies, bodies2);
      double potNrm = verify.getNrmScalar(bodies);
      double accDif = verify.getDifVector(bodies, bodies2);
      double accNrm = verify.getNrmVector(bodies);
      double potDifGlob, potNrmGlob, accDifGlob, accNrmGlob;
      MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&accDif, &accDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&accNrm, &accNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      double potRel = std::sqrt(potDifGlob/potNrmGlob);
      double accRel = std::sqrt(accDifGlob/accNrmGlob);
      logger::printTitle("FMM vs. direct");
      verify.print("Rel. L2 Error (pot)",potRel);
      verify.print("Rel. L2 Error (acc)",accRel);
      localTree.printTreeData(cells);
      traversal.printTraversalData();
      logger::printPAPI();
      bodies = buffer;
      if (!baseMPI.mpirank) {
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, potRel, accRel);
      }
      MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
      if (pass) {
        if (verify.verbose) std::cout << "passed accuracy regression at t: " << t << std::endl;
        if (args.accuracy) break;
        t = -1;
        isTime = true;
      }
    } else {
      double totalFMMGlob;
      MPI_Reduce(&totalFMM, &totalFMMGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      totalFMMGlob /= baseMPI.mpisize;
      if (!baseMPI.mpirank) {
        pass = verify.regression(args.getKey(baseMPI.mpisize), isTime, t, totalFMMGlob);
      }
      MPI_Bcast(&pass, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
      if (pass) {
        if (verify.verbose) std::cout << "passed time regression at t: " << t << std::endl;
        break;
      }
    }
    data.initTarget(bodies);
  }
  if (!pass) {
    if (verify.verbose) {
      if(!isTime) std::cout << "failed accuracy regression" << std::endl;
      else std::cout << "failed time regression" << std::endl;
    }
    abort();
  }
  Kernel::finalize();
}

template<int P>
struct CallFMM {
  static inline void LaplaceCartesianCPU_P1(Args args) {
    if(args.P == P) fmm<LaplaceCartesianCPU<P,1> >(args);
    CallFMM<P-1>::LaplaceCartesianCPU_P1(args);
  }
  static inline void LaplaceCartesianCPU_P0(Args args) {
    if(args.P == P) fmm<LaplaceCartesianCPU<P,0> >(args);
    CallFMM<P-1>::LaplaceCartesianCPU_P0(args);
  }
  static inline void LaplaceSphericalCPU_P(Args args) {
    if(args.P == P) fmm<LaplaceSphericalCPU<P> >(args);
    CallFMM<P-1>::LaplaceSphericalCPU_P(args);
  }
  static inline void HelmholtzSphericalCPU_P(Args args) {
    if(args.P == P) fmm<HelmholtzSphericalCPU<P> >(args);
    CallFMM<P-1>::HelmholtzSphericalCPU_P(args);
  }
  static inline void BiotSavartSphericalCPU_P(Args args) {
    if(args.P == P) fmm<BiotSavartSphericalCPU<P> >(args);
    CallFMM<P-1>::BiotSavartSphericalCPU_P(args);
  }
};

template<>
struct CallFMM<Pmin-1> {
  static inline void LaplaceCartesianCPU_P1(Args args) {
    if(args.P < Pmin || Pmax < args.P) {
      fprintf(stderr,"Pmin <= P <= Pmax\n");
      abort();
    }
  }
  static inline void LaplaceCartesianCPU_P0(Args args) {
    if(args.P < Pmin || Pmax < args.P) {
      fprintf(stderr,"Pmin <= P <= Pmax\n");
      abort();
    }
  }
  static inline void LaplaceSphericalCPU_P(Args args) {
    if(args.P < Pmin || Pmax < args.P) {
      fprintf(stderr,"Pmin <= P <= Pmax\n");
      abort();
    }
  }
  static inline void HelmholtzSphericalCPU_P(Args args) {
    if(args.P < Pmin || 2*Pmax < args.P) {
      fprintf(stderr,"Pmin <= P <= 2*Pmax\n");
      abort();
    }
  }
  static inline void BiotSavartSphericalCPU_P(Args args) {
    if(args.P < Pmin || Pmax < args.P) {
      fprintf(stderr,"Pmin <= P <= Pmax\n");
      abort();
    }
  }
};

int main(int argc, char ** argv) {
  Args args(argc, argv);                                        // Argument parser class
  switch (args.equation[0]) {                                   // Case switch for equation
  case 'L':                                                     // Laplace equation
    switch (args.basis[0]) {                                    //  Case switch for basis
    case 'C':                                                   //  Cartesian basis
      if (args.mass)                                            //   If all charges are positive
        CallFMM<Pmax>::LaplaceCartesianCPU_P1(args);            //    Call Laplace Cartesian kernel for mass
      else                                                      //   Elseif charges are both positive and negative
        CallFMM<Pmax>::LaplaceCartesianCPU_P0(args);            //    Call Laplace Cartesian kernel for charge
      break;                                                    //  Break Cartesian basis
    case 'S':                                                   //  Spherical basis
      CallFMM<Pmax>::LaplaceSphericalCPU_P(args);               //   Call Laplace Spherical kernel
      break;                                                    //  Break Spherical basis
    default:                                                    //  No matching case
      fprintf(stderr,"No matching basis\n");                    //   Print error message
      abort();                                                  //   Abort execution
    }                                                           //  End case switch for basis
    break;                                                      // Break Laplace equation
  case 'H':                                                     // Helmholtz equation
    CallFMM<2*Pmax>::HelmholtzSphericalCPU_P(args);             //  Call Helmholtz Spherical kernel
    break;                                                      // Break Helmholtz equation
  case 'B':                                                     // Biot-Savart equation
    CallFMM<Pmax>::BiotSavartSphericalCPU_P(args);              //  Call Biot-Savart Spherical kernel
    break;                                                      // Break Biot-Savart equation
  default:                                                      // No matching case
    fprintf(stderr,"No matching equation\n");                   //  Print error message
    abort();                                                    //  Abort execution
  }                                                             // End case switch for equation
  return 0;
}
