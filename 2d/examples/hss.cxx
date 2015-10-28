#include "localessentialtree.h"
#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "sort.h"
#include "traversal.h"
#include "updownpass.h"
#include "StrumpackDensePackage.hpp"

/* Laplace, cartesian coordinates example, 2D geometry
 *
 * Run with mpirun -np {#processes} ./laplace --numBodies {matrixsize}
 *
 */

void elements(void *, int *, int *, double *, int *);

const real_t eps2 = 0.0;

double elemops = 0.0;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  Dataset data;
  Logger logger;
  Sort sort;

  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(args.nspawn);
  BuildTree tree(args.ncrit,args.nspawn);
  UpDownPass pass(args.theta,eps2);
  Traversal traversal(args.nspawn,args.images,eps2);
  LocalEssentialTree LET(args.images);
  args.numBodies /= LET.mpisize;
  args.verbose = false;
  args.verbose &= LET.mpirank == 0;
  if (args.verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
    LET.verbose = true;
  }
  int myid = LET.mpirank;
  int np = LET.mpisize;

  logger.printTitle("FMM Parameters");
  args.print(logger.stringLength,P);
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();
  Bodies bodies = data.initBodies(args.numBodies, args.distribution, LET.mpirank, LET.mpisize);
  Bounds localBounds = boundbox.getBounds(bodies);
  Bounds globalBounds = LET.allreduceBounds(localBounds);
  localBounds = LET.partition(bodies,globalBounds);
  bodies = sort.sortBodies(bodies);
  bodies = LET.commBodies(bodies);
  Cells cells = tree.buildTree(bodies, localBounds);
  pass.upwardPass(cells);
  LET.setLET(cells,localBounds,cycle);
  LET.commBodies();
  LET.commCells();
  traversal.dualTreeTraversal(cells, cells, cycle, args.mutual);
  Bodies jbodies = bodies;
  Cells jcells;
  for (int irank=1; irank<LET.mpisize; irank++) {
    LET.getLET(jcells,(LET.mpirank+irank)%LET.mpisize);
    traversal.dualTreeTraversal(cells, jcells, cycle);
  }
  pass.downwardPass(cells);
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  logger.printTitle("MPI direct sum");
  Bodies buffer = bodies;
  data.sampleBodies(bodies, args.numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger.startTimer("Total Direct");
  for (int i=0; i<LET.mpisize; i++) {
    if (args.verbose) std::cout << "Direct loop          : " << i+1 << "/" << LET.mpisize << std::endl;
    LET.shiftBodies(jbodies);
    traversal.direct(bodies, jbodies, cycle);
  }
  traversal.normalize(bodies);
  logger.printTitle("Total runtime");
  logger.printTime("Total FMM");
  logger.stopTimer("Total Direct");
  boundbox.writeTime(LET.mpirank);
  tree.writeTime(LET.mpirank);
  pass.writeTime(LET.mpirank);
  traversal.writeTime(LET.mpirank);
  LET.writeTime(LET.mpirank);
  LET.writeSendCount();
  boundbox.resetTimer();
  tree.resetTimer();
  pass.resetTimer();
  traversal.resetTimer();
  LET.resetTimer();
  logger.resetTimer();
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  data.evalError(bodies2, bodies, diff1, norm1);
  MPI_Reduce(&diff1, &diff2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  logger.printTitle("FMM vs. direct");
  logger.printError(diff2, norm2);
  tree.printTreeData(cells);
  traversal.printTraversalData();
  Bodies gbodies = LET.allgatherBodies(buffer);
  bodies = gbodies;
  jbodies = gbodies;

  /* BLACS 2D grid, as square as possible */
  int ctxt;
  int nprow, npcol;
  int myrow, mycol;
  nprow=floor(sqrt((float)np));
  npcol=np/nprow;
  blacs_get_(&IZERO,&IZERO,&ctxt);
  blacs_gridinit_(&ctxt,"R",&nprow,&npcol);
  blacs_gridinfo_(&ctxt,&nprow,&npcol,&myrow,&mycol);

  int n=args.numBodies * LET.mpisize;
  int nb=64;

  /* Initialize the solver */
  StrumpackDensePackage<double,double> sdp(MPI_COMM_WORLD);
  sdp.verbose=true;
  sdp.use_HSS=true;
  sdp.tol_HSS=1e-4;
  sdp.block_HSS=64;

  double *A, *R, *S;
  int descA[BLACSCTXTSIZE], descRS[BLACSCTXTSIZE];
  int nrand=60;
  sdp.min_rand_HSS=nrand;
  sdp.max_rand_HSS=nrand;

  Bodies **pbodies=new Bodies*[2];
  pbodies[0]=&bodies;
  pbodies[1]=&jbodies;
  sdp.obj=(void*)pbodies;

  int seed=time(NULL);
  MPI_Bcast((void*)&seed,IONE,MPI_INTEGER,IZERO,MPI_COMM_WORLD);
  srand(seed);

  if(!myid) std::cout << "Sampling with " << nrand << " random vectors..." << std::endl;
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

#if 0
    int nbA=128;
    int r=0;
    while(r<locr) {
      int nrows=std::min(nbA,locr-r);
      A=new double[nrows*n];
      for(int j=0;j<n;j++) {
	B_iter Bj=jbodies.begin()+j;
	for(int i=0;i<nrows;i++) {
	  int locri=r+i+1;
	  int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	  B_iter Bi=bodies.begin()+globri-1;
	  vec2 dX=Bi->X-Bj->X;
	  real_t R2=norm(dX)+eps2;
	  A[i+nrows*j]=R2==0?0.0:-log(sqrt(R2));
	}
      }
      gemm('N','N',nrows,locc,n,1.0,A,nrows,Rglob,n,0.0,&S[r],locr);
      delete[] A;
      r+=nbA;
    }
    A=NULL;
#else
    for(int i=0;i<locr;i++) {
      int locri=i+1;
      int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
      bodies[i]=jbodies[globri-1];
      bodies[i].SRC=1;
      bodies[i].IBODY=i;
    }
    bodies.resize(locr);
    for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
      Bj->IBODY = Bj - jbodies.begin();
    }
    localBounds = boundbox.getBounds(bodies);
    localBounds = boundbox.getBounds(jbodies,localBounds);
    cells = tree.buildTree(bodies, localBounds);
    jcells = tree.buildTree(jbodies, localBounds);
    for(int k=1;k<=nrand;k++) {
      if(mycol==indxg2p_(&k,&nb,&mycol,&IZERO,&npcol)) {
	int lock=indxg2l_(&k,&nb,&mycol,&IZERO,&npcol);
	for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
	  int j = Bj-jbodies.begin();
	  Bj->SRC = Rglob[j+n*(lock-1)];
	}
	for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
	  Bi->TRG = 0;
	}
#if 0
	traversal.direct(bodies, jbodies, cycle);
#else
	pass.upwardPass(cells);
	pass.upwardPass(jcells);
	traversal.dualTreeTraversal(cells, jcells, cycle, false);
	pass.downwardPass(cells);
#endif
	for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
	  int i = Bi->IBODY;
	  S[i+locr*(lock-1)] = Bi->TRG;
	}
      }
    }
    bodies = gbodies;
    jbodies = gbodies;
#endif

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
  if(!myid) std::cout << "Sampling time: " << tend-tstart << " seconds" << std::endl << std::endl;

  sdp.compress(R,R,S,S,descRS,elements);

  double minops;
  MPI_Allreduce((void*)&elemops,(void*)&minops,IONE,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  minops/=1e9;
  double maxops;
  MPI_Allreduce((void*)&elemops,(void*)&maxops,IONE,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  maxops/=1e9;
  double sumops;
  MPI_Allreduce((void*)&elemops,(void*)&sumops,IONE,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  sumops/=1e9;

  if(!myid)
    std::cout << "Flops in kernel evaluation (x1e9): min=" << minops << "; max=" << maxops << "; total=" << sumops << std::endl << std::endl;

  /* Matrix-vector product */
  int nrhs=1;
  double *X, *B, *Btrue;
  int descXB[BLACSCTXTSIZE];
  if(myid<nprow*npcol) {
    int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    int locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
    X=new double[locr*locc];
    B=new double[locr*locc];
    Btrue=new double[locr*locc];
    for(int i=0;i<locr*locc;i++)
      X[i]=rand()/(double)RAND_MAX;
    int ierr;
    int dummy=std::max(1,locr);
    descinit_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
  } else {
    X=NULL;
    B=NULL;
    Btrue=NULL;
    descset_(descXB,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  }

  sdp.product('N',1.0,A,descA,X,descXB,0.0,B,descXB);
  sdp.print_statistics();

  /* Error check w/o a matrix: perform a direct sum */
  if(!myid) std::cout << std::endl << "Direct sum..." << std::endl;
  tstart=MPI_Wtime();

  /* BLAS2 direct sum.
   * Xglob is a duplicated version of X.
   * Then, for a given blocksize nbA, the processes
   * that have a piece of X (locr*locc>0) assemble
   * a block B of nbA rows of A, and they
   * compute Btrue([nbA rows],:)=B*Xglob.
   *
   * Construction of Xglob: first X is gathered onto
   * id 0, then it is broadcasted to all the processes.
   */
  double *Xglob=new double[n*nrhs];

  /* We initialize a context with only id 0 */
  int ctxtcent;
  blacs_get_(&IZERO,&IZERO,&ctxtcent);
  blacs_gridinit_(&ctxtcent,"R",&IONE,&IONE);

  /* We initialize a context with all the processes */
  int ctxtglob;
  blacs_get_(&IZERO,&IZERO,&ctxtglob);
  blacs_gridinit_(&ctxtglob,"R",&IONE,&np);

  /* Redistribute X onto id 0 (Xglob) */
  int descXglob[BLACSCTXTSIZE], ierr;
  if(!myid)
    descinit_(descXglob,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&ctxtcent,&n,&ierr);
  else
    descset_(descXglob,&n,&nrhs,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
  pgemr2d(n,nrhs,X,IONE,IONE,descXB,Xglob,IONE,IONE,descXglob,ctxtglob);

  /* Broadcast Xglob */
  MPI_Bcast((void*)Xglob,n*nrhs,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* Direct sum Btrue=A*Xglob */
  int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
  int locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
  if(locr*locc) {
#if 0
    int nbA=128;
    int r=0;
    int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    while(r<locr) {
      int nrows=std::min(nbA,locr-r);
      A=new double[nrows*n];
      for(int j=0;j<n;j++) {
	B_iter Bj=jbodies.begin()+j;
	for(int i=0;i<nrows;i++) {
	  int locri=r+i+1;
	  int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	  B_iter Bi=bodies.begin()+globri-1;
	  vec2 dX=Bi->X-Bj->X;
	  real_t R2=norm(dX)+eps2;
	  A[i+nrows*j]=R2==0?0.0:-log(sqrt(R2));
	}
      }
      gemm('N','N',nrows,nrhs,n,1.0,A,nrows,Xglob,n,0.0,&Btrue[r],locr);
      delete[] A;
      r+=nbA;
    }
    A=NULL;
#else
    for(int i=0;i<locr;i++) {
      int locri=i+1;
      int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
      bodies[i]=jbodies[globri-1];
      bodies[i].SRC=1;
      bodies[i].IBODY=i;
    }
    bodies.resize(locr);
    for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
      Bj->IBODY = Bj - jbodies.begin();
    }
    localBounds = boundbox.getBounds(bodies);
    localBounds = boundbox.getBounds(jbodies,localBounds);
    cells = tree.buildTree(bodies, localBounds);
    jcells = tree.buildTree(jbodies, localBounds);
    for(int k=1;k<=nrhs;k++) {
      if(mycol==indxg2p_(&k,&nb,&mycol,&IZERO,&npcol)) {
	int lock=indxg2l_(&k,&nb,&mycol,&IZERO,&npcol);
	for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
	  int j = Bj-jbodies.begin();
	  Bj->SRC = Xglob[j+n*(lock-1)];
	}
	for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
	  Bi->TRG = 0;
	}
#if 0
        traversal.direct(bodies, jbodies, cycle);
#else
        pass.upwardPass(cells);
        pass.upwardPass(jcells);
        traversal.dualTreeTraversal(cells, jcells, cycle, false);
        pass.downwardPass(cells);
#endif
	for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
	  int i = Bi->IBODY;
	  Btrue[i+locr*(lock-1)] = Bi->TRG;
	}
      }
    }
    bodies = gbodies;
    jbodies = gbodies;
#endif
  }
  delete[] Xglob;
  tend=MPI_Wtime();
  if(!myid) std::cout << "Direct sum time: " << tend-tstart << " seconds" << std::endl << std::endl;

  if(myid<nprow*npcol) {
    /* Compare Btrue against B */
    double err=plange('F',n,nrhs,Btrue,IONE,IONE,descXB,(double*)NULL);
    int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
    int locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
    for(int i=0;i<locr*locc;i++)
      Btrue[i]-=B[i];
    err=plange('F',n,nrhs,Btrue,IONE,IONE,descXB,(double*)NULL)/err;
    if(!myid) std::cout << "Product quality = " << err << std::endl;
  }

  delete[] R;
  delete[] S;
  delete[] X;
  delete[] B;
  delete[] Btrue;

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
      vec2 dX=Bi->X-Bj->X;
      real_t R2=norm(dX)+eps2;
      B[locr*(j-1)+(i-1)]=R2==0?0.0:-log(sqrt(R2));
    }
  }

  elemops+=7*locr*locc;
}

