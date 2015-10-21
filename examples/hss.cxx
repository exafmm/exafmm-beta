#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "verify.h"
#include "StrumpackDensePackage.hpp"
using namespace exafmm;

/* Laplace, cartesian coordinates example */

void elements(void *, int *, int *, double *, int *);

int main(int argc, char ** argv) {
  const real_t cycle = 2 * M_PI;
  Args args(argc, argv);
  Bodies bodies, bodies2, jbodies, buffer;
  BoundBox boundBox(args.nspawn);
  Bounds bounds;
  BuildTree buildTree(args.ncrit, args.nspawn);
  Cells cells, jcells;
  Dataset data;
  Traversal traversal(args.nspawn, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;
  num_threads(args.threads);

  int mpirank, mpisize;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpisize);

  kernel::eps2 = 0.0;
#if EXAFMM_HELMHOLTZ
  kernel::wavek = complex_t(10.,1.) / real_t(2 * M_PI);
#endif
  kernel::setup();
  logger::verbose = args.verbose && !mpirank;
  logger::printTitle("FMM Parameters");
  if(!mpirank)args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  buffer.reserve(bodies.size());
  if (args.IneJ) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    jbodies = data.initBodies(args.numBodies, args.distribution, 1);
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
  }
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    logger::startDAG();
    bounds = boundBox.getBounds(bodies);
    if (args.IneJ) {
      bounds = boundBox.getBounds(jbodies, bounds);
    }
    cells = buildTree.buildTree(bodies, buffer, bounds);
    upDownPass.upwardPass(cells);
    traversal.initListCount(cells);
    traversal.initWeight(cells);
    if (args.IneJ) {
      jcells = buildTree.buildTree(jbodies, buffer, bounds);
      upDownPass.upwardPass(jcells);
      traversal.traverse(cells, jcells, cycle, args.dual, false);
    } else {
      traversal.traverse(cells, cells, cycle, args.dual, args.mutual);
      jbodies = bodies;
    }
    upDownPass.downwardPass(cells);
    logger::printTitle("Total runtime");
    logger::stopDAG();
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
    logger::resetTimer("Total FMM");
    if (args.write) {
      logger::writeTime();
    }
    traversal.writeList(cells, 0);
    const int numTargets = 100;
    buffer = bodies;
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
    bodies = buffer;
    data.initTarget(bodies);
  }
  if (args.getMatrix) {
    /* BLACS 2D grid, as square as possible */
    int blacsctxt;
    blacs_get_(&IZERO, &IZERO, &blacsctxt);
    int rowsize = floor(sqrt((float) mpisize));
    int colsize = mpisize / rowsize;
    blacs_gridinit_(&blacsctxt, "R", &rowsize, &colsize);
    int rowrank, colrank;
    blacs_gridinfo_(&blacsctxt, &rowsize, &colsize, &rowrank, &colrank);

    int numBodies = bodies.size();
    int numBlocks = 64;

    /* Initialize the solver */
    StrumpackDensePackage<double,double> sdp(MPI_COMM_WORLD);
    sdp.verbose = true;
    sdp.use_HSS = true;
    sdp.tol_HSS = 1e-4;
    sdp.levels_HSS = 5;

    double *A;
    int descA[BLACSCTXTSIZE], descRS[BLACSCTXTSIZE];
    int nrand = std::min(4*(int)floor(sqrt(numBodies)), 1000);
    sdp.min_rand_HSS = nrand;
    sdp.max_rand_HSS = nrand;

    Bodies **pbodies = new Bodies*[2];
    pbodies[0] = &bodies;
    pbodies[1] = &bodies;
    sdp.obj =(void*) pbodies;

    int seed = time(NULL);
    MPI_Bcast((void*)&seed, IONE, MPI_INTEGER, IZERO, MPI_COMM_WORLD);
    srand(seed);

    if(!mpirank) std::cout << "Sampling..." << std::endl;
    double tstart = MPI_Wtime();
    assert(mpirank < rowsize * colsize);
    int locr = numroc_(&numBodies, &numBlocks, &rowrank, &IZERO, &rowsize);
    int locc = numroc_(&nrand, &numBlocks, &colrank, &IZERO, &colsize);
    double * R = new double [locr*locc]();
    double * S = new double [locr*locc]();

    /* BLAS3 sampling.
     * In Rglob, processes store only the columns
     * of the random vectors that are relevant to them.
     * We enforce the same sequence of calls to the
     * random number generator on every process.
     * Then, for a given blocksize nbA, they
     * assemble a block B of nbA rows of A, and they
     * compute S([nbA rows],:)=B*Rglob.
     */

    double * Rglob = new double[numBodies*locc];
    for (int j=1; j<=nrand; j++) {
      if (colrank==indxg2p_(&j, &numBlocks, &colrank, &IZERO, &colsize)) {
	int locj = indxg2l_(&j, &numBlocks, &colrank, &IZERO, &colsize);
	for (int i=0; i<numBodies; i++)
	  Rglob[i+numBodies*(locj-1)] = rand() / (double)RAND_MAX;
      } else {
	for (int i=0; i<numBodies; i++)
	  rand();
      }
    } 

    int nbA = 128;
    int r = 0;
    while (r<locr) {
      int nrows = std::min(nbA,locr-r);
      A = new double [nrows*numBodies];
      for (int j=0; j<numBodies; j++) {
	B_iter Bj = jbodies.begin()+j; 
	for (int i=0; i<nrows; i++) {
	  int locri = r + i + 1;
	  int globri = indxl2g_(&locri, &numBlocks, &rowrank, &IZERO, &rowsize);
	  B_iter Bi = bodies.begin() + globri - 1; 
	  vec3 dX = Bi->X - Bj->X;
	  real_t R2 = norm(dX) + kernel::eps2;
	  A[i+nrows*j] = R2 == 0 ? 0.0 : 1.0 / sqrt(R2);
	}
      }
      gemm('N', 'N', nrows, locc, numBodies, 1.0, A, nrows, Rglob, numBodies, 0.0, &S[r], locr);
      delete[] A;
      r += nbA;
    }
    A=NULL;

    /* Compress the random vectors */
    for (int j=0; j<locc; j++) {
      for (int i=1; i<=numBodies; i++) {
	if (rowrank==indxg2p_(&i, &numBlocks, &rowrank, &IZERO, &rowsize)) {
	  int loci = indxg2l_(&i, &numBlocks, &rowrank, &IZERO, &rowsize);
	  R[loci-1+locr*j] = Rglob[i-1+numBodies*j];
	}
      }
    }
    delete[] Rglob;

    int ierr;
    int dummy=std::max(1,locr);
    descinit_(descA, &numBodies, &numBodies, &numBlocks, &numBlocks, &IZERO, &IZERO, &blacsctxt, &dummy, &ierr);
    descinit_(descRS, &numBodies, &nrand, &numBlocks, &numBlocks, &IZERO, &IZERO, &blacsctxt, &dummy, &ierr);

    double tend = MPI_Wtime();
    if(!mpirank) std::cout << "Sampling time: " << tend-tstart << " seconds" << std::endl;

    sdp.compress(R, R, S, S, descRS, elements);

    double *X, *B, *Btrue;
    int descXB[BLACSCTXTSIZE];
    if (mpirank < rowsize * colsize) {
      int locr=numroc_(&numBodies, &numBlocks, &rowrank, &IZERO, &rowsize);
      int locc=numroc_(&IONE, &numBlocks, &colrank, &IZERO, &colsize);
      X = new double [locr * locc];
      B = new double [locr * locc];
      Btrue = new double [locr * locc];
      for (int i=0; i<locr*locc; i++)
	X[i] = rand() / (double)RAND_MAX;
      int ierr;
      int dummy = std::max(1,locr);
      descinit_(descXB, &numBodies, &IONE, &numBlocks, &numBlocks, &IZERO, &IZERO, &blacsctxt, &dummy, &ierr);
    } else {
      X = NULL;
      B = NULL;
      Btrue = NULL;
      descset_(descXB, &numBodies, &IONE, &numBlocks, &numBlocks, &IZERO, &IZERO, &INONE, &IONE);
    }

    sdp.product('N', 1.0, A, descA, X, descXB, 0.0, B, descXB);
    sdp.print_statistics();

    if (mpirank < rowsize * colsize) {
      double err=plange('F', numBodies, 1, Btrue, IONE, IONE, descXB, (double*)NULL);
      int locr=numroc_(&numBodies, &numBlocks, &rowrank, &IZERO, &rowsize);
      int locc=numroc_(&IONE, &numBlocks, &colrank, &IZERO, &colsize);
      for (int i=0; i<locr*locc; i++)
	Btrue[i] -= B[i];
      err=plange('F', numBodies, 1, Btrue, IONE, IONE, descXB, (double*)NULL) / err;
      if(!mpirank) std::cout << "Product quality = " << err << std::endl;
    }

    delete[] R;
    delete[] S;
    delete[] X;
    delete[] B;
    delete[] Btrue;
  }
  logger::writeDAG();
  MPI_Finalize();
  return 0;
}

void elements(void * obj, int *I, int *J, double *B, int *descB) {
  if(B==NULL)
    return;

  int blacsctxt=descB[1];
  int mB=descB[2];
  int nB=descB[3];
  int mb=descB[4];
  int nb=descB[5];
  int rsrc=descB[6];
  int csrc=descB[7];

  Bodies **pbodies=(Bodies **)obj;
  Bodies *bodies=pbodies[0];
  Bodies *jbodies=pbodies[1];

  int rowsize, colsize;
  int rowrank, colrank;
  blacs_gridinfo_(&blacsctxt,&rowsize,&colsize,&rowrank,&colrank);

  int locr=numroc_(&mB,&mb,&rowrank,&rsrc,&rowsize);
  int locc=numroc_(&nB,&nb,&colrank,&csrc,&colsize);

  for(int i=1;i<=locr;i++) {
    for(int j=1;j<=locc;j++) {
      int ii=indxl2g_(&i,&mb,&rowrank,&csrc,&rowsize);
      int jj=indxl2g_(&j,&nb,&colrank,&csrc,&colsize);
      int iii=I[ii-1];
      int jjj=J[jj-1];
      B_iter Bi=bodies->begin()+iii-1;
      B_iter Bj=jbodies->begin()+jjj-1;
      vec3 dX=Bi->X-Bj->X;
      real_t R2=norm(dX)+kernel::eps2;
      B[locr*(j-1)+(i-1)]=R2==0?0.0:1.0/sqrt(R2); 
    }
  }
}

