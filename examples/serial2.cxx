#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "utils.h"
#include "verify.h"
#if VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Bodies bodies, bodies2, bodies3, jbodies;
  BoundBox boundBox(args.nspawn);
  Bounds bounds;
  BuildTree buildTree(args.ncrit, args.nspawn);
  Cells cells, jcells;
  Dataset data;
  Traversal traversal(args.nspawn, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;
  num_threads(args.threads);

  int N = args.numBodies;
  int maxlev = 10;
  int maxheight = 2;
  int nbins = (1 << maxlev);
  float* X = (float*)malloc(3 * N * sizeof(float));
  float* Y = (float*)malloc(11 * N * sizeof(float));
  float* Z = (float*)malloc(11 * N * sizeof(float));
  uint* codes = (uint*)malloc(3 * N * sizeof(uint));
  unsigned long int* mcodes = (unsigned long int*)malloc(N * sizeof(unsigned long int));
  uint long* scodes = (uint long*)malloc(N * sizeof(uint long));
  uint* pointIds = (uint*)malloc(N * sizeof(uint));
  uint* index = (uint*)malloc(N * sizeof(uint));
  uint long* bins = (uint long*)malloc(N * sizeof(uint long*));
  int* levels = (int*)malloc(N * sizeof(int));

  const real_t cycle = 2 * M_PI;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    logger::startDAG();
    bounds = boundBox.getBounds(bodies);
#if 0
    cells = buildTree.buildTree(bodies, bounds);
#else
    logger::startTimer("Copy in");
    B_iter B = bodies.begin();
    for (int i=0; i<N; i++, B++) {
      index[i] = i;
      for (int d=0; d<3; d++) {
	X[3*i+d] = B->X[d];
      }
      Y[11*i+ 0] = B->X[0];
      Y[11*i+ 1] = B->X[1];
      Y[11*i+ 2] = B->X[2];
      Y[11*i+ 3] = B->SRC;
      Y[11*i+ 4] = B->IBODY;
      Y[11*i+ 5] = B->IRANK;
      Y[11*i+ 6] = B->WEIGHT;
      Y[11*i+ 7] = B->TRG[0];
      Y[11*i+ 8] = B->TRG[1];
      Y[11*i+ 9] = B->TRG[2];
      Y[11*i+10] = B->TRG[3];
    }
    logger::stopTimer("Copy in");
    logger::startTimer("Quantization");
    compute_quantization_codes_T(codes, X, N, nbins);
    logger::stopTimer("Quantization"); 
    logger::startTimer("Morton");
    morton_encoding_T(mcodes, codes, N);
    logger::stopTimer("Morton");
    logger::startTimer("Radix sort");
    bin_sort_radix6(mcodes, scodes, pointIds, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));
    logger::stopTimer("Radix sort");
    logger::startTimer("Permutation");
    rearrange_dataTL(Z, Y, pointIds, N);
    logger::stopTimer("Permutation");
    logger::startTimer("Copy out");
    B = bodies.begin();
    for (int i=0; i<N; i++, B++) {
      B->X[0]   = Z[11*i+ 0];
      B->X[1]   = Z[11*i+ 1];
      B->X[2]   = Z[11*i+ 2];
      B->SRC    = Z[11*i+ 3];
      B->IBODY  = Z[11*i+ 4];
      B->IRANK  = Z[11*i+ 5];
      B->WEIGHT = Z[11*i+ 6];
      B->TRG[0] = Z[11*i+ 7];
      B->TRG[1] = Z[11*i+ 8];
      B->TRG[2] = Z[11*i+ 9];
      B->TRG[3] = Z[11*i+10];
      std::cout << i << " " << B->X << std::endl;
    }
    logger::stopTimer("Copy out");
    /*
    int ibin = -1;
    for (int i=0; i<N; i++) {
      if(bins[i] != ibin) {
        if (ibin<100) std::cout << i << " " << bins[i] << std::endl;
	ibin = bins[i];
      }
    }
    */
#endif
#if 0
    upDownPass.upwardPass(cells);
    traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
    jbodies = bodies;
    upDownPass.downwardPass(cells);
    logger::printTitle("Total runtime");
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
    logger::resetTimer("Total FMM");
#if WRITE_TIME
    logger::writeTime();
#endif
    const int numTargets = 100;
    bodies3 = bodies;
    data.sampleBodies(bodies, numTargets);
    bodies2 = bodies;
    data.initTarget(bodies);
    logger::startTimer("Total Direct");
    traversal.direct(bodies, jbodies, cycle);
    traversal.normalize(bodies);
    logger::stopTimer("Total Direct");
    double potDif = verify.getDifScalar(bodies, bodies2);
    double potNrm = verify.getNrmScalar(bodies);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    logger::printTitle("FMM vs. direct");
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
    buildTree.printTreeData(cells);
    traversal.printTraversalData();
    logger::printPAPI();
    logger::stopDAG();
    bodies = bodies3;
    data.initTarget(bodies);
#endif
  }
  logger::writeDAG();
#if VTK
  B_iter B = bodies.begin();
  for (int i=0; i<N; i++, B++) {
    B->IBODY = bins[i];
  }
  vtk3DPlot vtk;
  vtk.setBounds(M_PI,0);
  vtk.setGroupOfPoints(bodies);
  vtk.plot();
#endif
  return 0;
}
