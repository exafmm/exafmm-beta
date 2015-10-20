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

  int myid, np;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Comm_size(MPI_COMM_WORLD,&np);

  kernel::eps2 = 0.0;
#if EXAFMM_HELMHOLTZ
  kernel::wavek = complex_t(10.,1.) / real_t(2 * M_PI);
#endif
  kernel::setup();
  logger::verbose = args.verbose && !myid;
  logger::printTitle("FMM Parameters");
  if(!myid)args.print(logger::stringLength, P);
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
    int ctxt;
    int nprow, npcol;
    int myrow, mycol;
    nprow=floor(sqrt((float)np));
    npcol=np/nprow;
    blacs_get_(&IZERO,&IZERO,&ctxt);
    blacs_gridinit_(&ctxt,"R",&nprow,&npcol);
    blacs_gridinfo_(&ctxt,&nprow,&npcol,&myrow,&mycol);

    int n=bodies.size();
    int nb=64;

    /* Initialize the solver */
    StrumpackDensePackage<double,double> sdp(MPI_COMM_WORLD);
    sdp.verbose=true;
    sdp.use_HSS=true;
    sdp.tol_HSS=1e-4;
    sdp.levels_HSS=5;

    double *A, *R, *S;
    int descA[BLACSCTXTSIZE], descRS[BLACSCTXTSIZE];
    /* Matrix-free version.
     * The number of random vectors is fixed.
     */ 
    int nrand=std::min(4*(int)floor(sqrt(n)),1000);
    sdp.min_rand_HSS=nrand;
    sdp.max_rand_HSS=nrand;

    Bodies **pbodies=new Bodies*[2];
    pbodies[0]=&bodies;
    pbodies[1]=&bodies;
    sdp.obj=(void*)pbodies;

    int seed=time(NULL);
    MPI_Bcast((void*)&seed,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD);
    srand(seed);

    if(!myid) std::cout << "Sampling..." << std::endl;
    double tstart=MPI_Wtime();
    if(myid<nprow*npcol) {
      int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
      int locc=numroc_(&nrand,&nb,&mycol,&IZERO,&npcol);
      R=new double[locr*locc]();
      S=new double[locr*locc]();

      /* BLAS3 sampling.
       * In Rglob, processes store only the columns
       * of the random vectors that are relevant to them.
       * We enforce the same sequence of calls to the
       * random number generator on every process.
       * Then, for a given blocksize nbA, they
       * assemble a block B of nbA rows of A, and they
       * compute S([nbA rows],:)=B*Rglob.
       */

      double *Rglob=new double[n*locc];
      for(int j=1;j<=nrand;j++) {
	if(mycol==indxg2p_(&j,&nb,&mycol,&IZERO,&npcol)) {
	  /* I own this column */
	  int locj=indxg2l_(&j,&nb,&mycol,&IZERO,&npcol);
	  for(int i=0;i<n;i++)
	    Rglob[i+n*(locj-1)]=rand()/(double)RAND_MAX;
	} else {
	  /* n calls to rand() to enforce the same calling
	   * sequence on every process.
	   */ 
	  for(int i=0;i<n;i++)
	    rand();
	}
      } 

      int nbA=128;
      int r=0;
      while(r<locr) {
	int nrows=std::min(nbA,locr-r);

	/* Compute the nrows rows of A that correspond
	 * to rows r:r+nrows-1 of the local array S.
	 */
	A=new double[nrows*n];
	for(int j=0;j<n;j++) {
	  B_iter Bj=jbodies.begin()+j; 
	  for(int i=0;i<nrows;i++) {
	    int locri=r+i+1;
	    int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	    B_iter Bi=bodies.begin()+globri-1; 
	    vec3 dX=Bi->X-Bj->X;
	    real_t R2=norm(dX)+kernel::eps2;
	    A[i+nrows*j]=R2==0?0.0:1.0/sqrt(R2);
	  }
	}

	/* Compute nrows of the sample */
	gemm('N','N',nrows,locc,n,1.0,A,nrows,Rglob,n,0.0,&S[r],locr);

	delete[] A;
	r+=nbA;
      }
      A=NULL;

      /* Compress the random vectors */
      for(int j=0;j<locc;j++) {
	for(int i=1;i<=n;i++) {
	  if(myrow==indxg2p_(&i,&nb,&myrow,&IZERO,&nprow)) {
	    /* I own this row */
	    int loci=indxg2l_(&i,&nb,&myrow,&IZERO,&nprow);
	    R[loci-1+locr*j]=Rglob[i-1+n*j];
	  }
	}
      }
      delete[] Rglob;

      int ierr;
      int dummy=std::max(1,locr);
      descinit_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
      descinit_(descRS,&n,&nrand,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
    } else {
      R=NULL;
      S=NULL;
      descset_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
      descset_(descRS,&n,&nrand,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
    }
    double tend=MPI_Wtime();
    if(!myid) std::cout << "Sampling time: " << tend-tstart << " seconds" << std::endl;

    sdp.compress(R,R,S,S,descRS,elements);

    double *X, *B, *Btrue;
    int descXB[BLACSCTXTSIZE];
    if(myid<nprow*npcol) {
      int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
      int locc=numroc_(&IONE,&nb,&mycol,&IZERO,&npcol);
      X=new double[locr*locc];
      B=new double[locr*locc];
      Btrue=new double[locr*locc];
      for(int i=0;i<locr*locc;i++)
        X[i]=rand()/(double)RAND_MAX;
      int ierr;
      int dummy=std::max(1,locr);
      descinit_(descXB,&n,&IONE,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
    } else {
      X=NULL;
      B=NULL;
      Btrue=NULL;
      descset_(descXB,&n,&IONE,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
    }

    sdp.product('N',1.0,A,descA,X,descXB,0.0,B,descXB);
    sdp.print_statistics();

    if(myid<nprow*npcol) {
      double err=plange('F',n,1,Btrue,IONE,IONE,descXB,(double*)NULL);
      int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
      int locc=numroc_(&IONE,&nb,&mycol,&IZERO,&npcol);
      for(int i=0;i<locr*locc;i++)
	Btrue[i]-=B[i];
      err=plange('F',n,1,Btrue,IONE,IONE,descXB,(double*)NULL)/err;
      if(!myid) std::cout << "Product quality = " << err << std::endl;
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

  int ctxt=descB[1];
  int mB=descB[2];
  int nB=descB[3];
  int mb=descB[4];
  int nb=descB[5];
  int rsrc=descB[6];
  int csrc=descB[7];

  Bodies **pbodies=(Bodies **)obj;
  Bodies *bodies=pbodies[0];
  Bodies *jbodies=pbodies[1];

  int nprow, npcol;
  int myrow, mycol;
  blacs_gridinfo_(&ctxt,&nprow,&npcol,&myrow,&mycol);

  int locr=numroc_(&mB,&mb,&myrow,&rsrc,&nprow);
  int locc=numroc_(&nB,&nb,&mycol,&csrc,&npcol);

  for(int i=1;i<=locr;i++) {
    for(int j=1;j<=locc;j++) {
      int ii=indxl2g_(&i,&mb,&myrow,&csrc,&nprow);
      int jj=indxl2g_(&j,&nb,&mycol,&csrc,&npcol);
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

