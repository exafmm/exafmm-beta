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
#include "StrumpackDensePackage.hpp"
using namespace exafmm_helmholtz;
/* Helmholtz, spherical coordinates example, 3D geometry.
 *
 * Run with mpirun -np {#processes} ./helmholtz -DgGmovx -n {matrixsize}
 *
 */

void elements(void *, int *, int *, dcomplex *, int *);

double elemops = 0.0;

const complex_t I1(0.0,1.0);

int main(int argc, char ** argv) {
  const real_t cycle = 2 * M_PI;
  const real_t eps2 = 0.0;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);
  Args args(argc, argv);
  BaseMPI baseMPI;
  Bodies bodies, bodies2, jbodies, gbodies, buffer;
  BoundBox boundBox;
  Bounds localBounds, globalBounds;
  BuildTree localTree(args.ncrit);
  BuildTree globalTree(1);
  Cells cells, jcells, gcells;
  Dataset data;
  Kernel kernel(args.P, eps2, wavek);
  Partition partition(baseMPI);
  Traversal traversal(kernel, args.theta, args.nspawn, args.images, args.path);
  TreeMPI treeMPI(kernel, baseMPI, args.theta, args.images);
  UpDownPass upDownPass(kernel);
  Verify verify(args.path);
  num_threads(args.threads);

  int myid = baseMPI.mpirank;
  int np = baseMPI.mpisize;

  args.numBodies /= baseMPI.mpisize;
  args.verbose &= baseMPI.mpirank == 0;
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength);
  bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
  buffer.reserve(bodies.size());
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    localBounds = boundBox.getBounds(bodies);
    globalBounds = baseMPI.allreduceBounds(localBounds);
    partition.bisection(bodies, globalBounds);
    bodies = treeMPI.commBodies(bodies);
    localBounds = boundBox.getBounds(bodies);
    cells = localTree.buildTree(bodies, buffer, localBounds);
    localBounds = boundBox.getBounds(cells, localBounds);
    upDownPass.upwardPass(cells);
    jbodies = bodies;
    for (int irank=0; irank<baseMPI.mpisize; irank++) {
      treeMPI.shiftBodies(jbodies);
      jcells.clear();
      localBounds = boundBox.getBounds(jbodies);
      jcells = localTree.buildTree(jbodies, buffer, localBounds);
      upDownPass.upwardPass(jcells);
      traversal.traverse(cells, jcells, cycle, args.dual);
    }
    upDownPass.downwardPass(cells);

    logger::stopPAPI();
    logger::stopTimer("Total FMM", 0);
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
    logger::printTitle("FMM vs. direct");
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
    localTree.printTreeData(cells);
    traversal.printTraversalData();
    logger::printPAPI();
    bodies = buffer;
    data.initTarget(bodies);
    logger::resetTimer("Total FMM");
    if (args.write) {
      logger::writeTime(baseMPI.mpirank);
    }
    traversal.writeList(cells, baseMPI.mpirank);
  }

  if (args.getMatrix) {
    gbodies = treeMPI.allgatherBodies(bodies);
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

    int n=args.numBodies * baseMPI.mpisize;
    int nb=64;

    /* Initialize the solver */
    StrumpackDensePackage<dcomplex,double> sdp(MPI_COMM_WORLD);
    sdp.verbose=true;
    sdp.use_HSS=true;
    sdp.tol_HSS=1e-4;
    sdp.levels_HSS=5;

    dcomplex *A, *Rr, *Rc, *Sr, *Sc;
    int descA[BLACSCTXTSIZE], descRS[BLACSCTXTSIZE];
    int nrand=std::min(8*(int)floor(sqrt(n)),1000);
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
      Rr=new dcomplex[locr*locc]();
      Rc=new dcomplex[locr*locc]();
      Sr=new dcomplex[locr*locc]();
      Sc=new dcomplex[locr*locc]();

      /* BLAS3 sampling.
       * In Rglob, processes store only the columns
       * of the random vectors that are relevant to them.
       * We enforce the same sequence of calls to the
       * random number generator on every process.
       * Then, for a given blocksize nbA, they
       * assemble a block B of nbA rows of A, and they
       * compute Sr([nbA rows],:)=B*Rglob.
       *
       * For Sc, same thing, but with rows of A^H,
       * i.e, columns of A with switched sign of
       * imaginary part.
       */

      dcomplex *Rglob=new dcomplex[n*locc];
      for(int j=1;j<=nrand;j++) {
	if(mycol==indxg2p_(&j,&nb,&mycol,&IZERO,&npcol)) {
	  /* I own this column */
	  int locj=indxg2l_(&j,&nb,&mycol,&IZERO,&npcol);
	  for(int i=0;i<n;i++)
	    Rglob[i+n*(locj-1)]=dcomplex(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
	} else {
	  /* n calls to rand() to enforce the same calling
	   * sequence on every process.
	   */
	  for(int i=0;i<n;i++) {
	    rand();
	    rand();
	  }
	}
      }

#if 0
      int nbA=128;
      int r=0;
      while(r<locr) {
	int nrows=std::min(nbA,locr-r);
	A=new dcomplex[nrows*n];
	for(int j=0;j<n;j++) {
	  B_iter Bj=jbodies.begin()+j;
	  for(int i=0;i<nrows;i++) {
	    int locri=r+i+1;
	    int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	    B_iter Bi=bodies.begin()+globri-1;
	    vec3 dX=Bi->X-Bj->X;
	    real_t R2=norm(dX)+kernel.eps2;
	    real_t R=sqrt(R2);
	    A[i+nrows*j]=R2==0?0.0:exp(I1*(kernel.wavek)*R)/R;
	  }
	}
	gemm('N','N',nrows,locc,n,dcomplex(1.0),A,nrows,Rglob,n,dcomplex(0.0),&Sr[r],locr);
	delete[] A;
	r+=nbA;
      }
      A=NULL;
#else
      for(int i=0;i<locr;i++) {
        int locri=i+1;
        int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
        bodies[i]=gbodies[globri-1];
        bodies[i].SRC=1;
	bodies[i].IBODY=i;
      }
      bodies.resize(locr);
      for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
        Bj->IBODY = Bj - jbodies.begin();
      }
      localBounds = boundBox.getBounds(bodies);
      localBounds = boundBox.getBounds(jbodies,localBounds);
      for(int k=1;k<=nrand;k++) {
        if(mycol==indxg2p_(&k,&nb,&mycol,&IZERO,&npcol)) {
          int lock=indxg2l_(&k,&nb,&mycol,&IZERO,&npcol);
          for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
            int j = Bj->IBODY;
            Bj->SRC = Rglob[j+n*(lock-1)];
          }
          for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
            Bi->TRG = 0;
          }
#if 0
          traversal.direct(bodies, jbodies, cycle);
#else
	  cells = localTree.buildTree(bodies, buffer, localBounds);
	  jcells = localTree.buildTree(jbodies, buffer, localBounds);
          upDownPass.upwardPass(cells);
          upDownPass.upwardPass(jcells);
          traversal.traverse(cells, jcells, cycle, args.dual);
          upDownPass.downwardPass(cells);
#endif
          for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
            int i = Bi->IBODY;
            Sr[i+locr*(lock-1)] = Bi->TRG[0];
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
	    Rr[loci-1+locr*j]=Rglob[i-1+n*j];
	  }
	}
      }
      delete[] Rglob;

      Rglob=new dcomplex[n*locc];
      for(int j=1;j<=nrand;j++) {
	if(mycol==indxg2p_(&j,&nb,&mycol,&IZERO,&npcol)) {
	  /* I own this column */
	  int locj=indxg2l_(&j,&nb,&mycol,&IZERO,&npcol);
	  for(int i=0;i<n;i++)
	    Rglob[i+n*(locj-1)]=dcomplex(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
	} else {
	  /* n calls to rand() to enforce the same calling
	   * sequence on every process.
	   */
	  for(int i=0;i<n;i++) {
	    rand();
	    rand();
	  }
	}
      }

#if 0
      int nbA = 128;
      int r=0;
      while(r<locr) {
	int nrows=std::min(nbA,locr-r);
	A=new dcomplex[nrows*n];
	for(int j=0;j<n;j++) {
	  B_iter Bj=jbodies.begin()+j;
	  for(int i=0;i<nrows;i++) {
	    int locri=r+i+1;
	    int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	    B_iter Bi=bodies.begin()+globri-1;
	    vec3 dX=Bi->X-Bj->X;
	    real_t R2=norm(dX)+kernel.eps2;
	    real_t R=sqrt(R2);
	    A[i+nrows*j]=R2==0?0.0:exp(I1*(kernel.wavek)*R)/R;
	    A[i+nrows*j]=std::conj(A[i+nrows*j]);
	  }
	}
	gemm('N','N',nrows,locc,n,dcomplex(1.0),A,nrows,Rglob,n,dcomplex(0.0),&Sc[r],locr);
	delete[] A;
	r+=nbA;
      }
      A=NULL;
#else
      kernel.wavek = complex_t(-std::real(kernel.wavek),std::imag(kernel.wavek));
      for(int i=0;i<locr;i++) {
        int locri=i+1;
        int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
        bodies[i]=gbodies[globri-1];
        bodies[i].SRC=1;
	bodies[i].IBODY=i;
      }
      bodies.resize(locr);
      for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
        Bj->IBODY = Bj - jbodies.begin();
      }
      localBounds = boundBox.getBounds(bodies);
      localBounds = boundBox.getBounds(jbodies,localBounds);
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
          traversal.direct(bodies2, jbodies, cycle);
#else
          cells = localTree.buildTree(bodies, buffer, localBounds);
          jcells = localTree.buildTree(jbodies, buffer, localBounds);
          upDownPass.upwardPass(cells);
          upDownPass.upwardPass(jcells);
          traversal.traverse(cells, jcells, cycle, args.dual);
          upDownPass.downwardPass(cells);
#endif
          for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
            int i = Bi->IBODY;
            Sc[i+locr*(lock-1)] = Bi->TRG[0];
          }
        }
      }
      bodies = gbodies;
      jbodies = gbodies;
      kernel.wavek = complex_t(-std::real(kernel.wavek),std::imag(kernel.wavek));
#endif

      /* Compress the random vectors */
      for(int j=0;j<locc;j++) {
	for(int i=1;i<=n;i++) {
	  if(myrow==indxg2p_(&i,&nb,&myrow,&IZERO,&nprow)) {
	    /* I own this row */
	    int loci=indxg2l_(&i,&nb,&myrow,&IZERO,&nprow);
	    Rc[loci-1+locr*j]=Rglob[i-1+n*j];
	  }
	}
      }
      delete[] Rglob;

      int ierr;
      int dummy=std::max(1,locr);
      descinit_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
      descinit_(descRS,&n,&nrand,&nb,&nb,&IZERO,&IZERO,&ctxt,&dummy,&ierr);
    } else {
      Rr=NULL;
      Rc=NULL;
      Sr=NULL;
      Sc=NULL;
      descset_(descA,&n,&n,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
      descset_(descRS,&n,&nrand,&nb,&nb,&IZERO,&IZERO,&INONE,&IONE);
    }
    double tend=MPI_Wtime();
    if(!myid) std::cout << "Sampling time: " << tend-tstart << " seconds" << std::endl << std::endl;

    sdp.compress(Rr,Rc,Sr,Sc,descRS,elements);

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
    dcomplex *X, *B, *Btrue;
    int descXB[BLACSCTXTSIZE];
    srand(time(NULL)+myid); /* Now processes can follow different random sequences */
    if(myid<nprow*npcol) {
      int locr=numroc_(&n,&nb,&myrow,&IZERO,&nprow);
      int locc=numroc_(&nrhs,&nb,&mycol,&IZERO,&npcol);
      X=new dcomplex[locr*locc];
      B=new dcomplex[locr*locc];
      Btrue=new dcomplex[locr*locc];
      for(int i=0;i<locr*locc;i++)
        X[i]=dcomplex(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
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
    dcomplex *Xglob=new dcomplex[n*nrhs];

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
    MPI_Bcast((void*)Xglob,n*nrhs,MPI_C_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

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
	A=new dcomplex[nrows*n];
	for(int j=0;j<n;j++) {
	  B_iter Bj=jbodies.begin()+j;
	  for(int i=0;i<nrows;i++) {
	    int locri=r+i+1;
	    int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
	    B_iter Bi=bodies.begin()+globri-1;
	    vec3 dX=Bi->X-Bj->X-kernel.Xperiodic;
	    real_t R2=norm(dX)+kernel.eps2;
	    real_t R=sqrt(R2);
	    A[i+nrows*j]=R2==0?0.0:exp(I1*(kernel.wavek)*R)/R;
	  }
	}
	gemm('N','N',nrows,nrhs,n,dcomplex(1.0),A,nrows,Xglob,n,dcomplex(0.0),&Btrue[r],locr);
	delete[] A;
	r+=nbA;
      }
      A=NULL;
#else
      for(int i=0;i<locr;i++) {
        int locri=i+1;
        int globri=indxl2g_(&locri,&nb,&myrow,&IZERO,&nprow);
        bodies[i]=gbodies[globri-1];
        bodies[i].SRC=1;
	bodies[i].IBODY=i;
      }
      bodies.resize(locr);
      for(B_iter Bj=jbodies.begin(); Bj!=jbodies.end(); Bj++) {
        Bj->IBODY = Bj - jbodies.begin();
      }
      localBounds = boundBox.getBounds(bodies);
      localBounds = boundBox.getBounds(jbodies,localBounds);
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
          cells = localTree.buildTree(bodies, buffer, localBounds);
          jcells = localTree.buildTree(jbodies, buffer, localBounds);
          upDownPass.upwardPass(cells);
          upDownPass.upwardPass(jcells);
          traversal.traverse(cells, jcells, cycle, args.dual);
          upDownPass.downwardPass(cells);
#endif
          for(B_iter Bi=bodies.begin(); Bi!=bodies.end(); Bi++) {
            int i = Bi->IBODY;
            Btrue[i+locr*(lock-1)] = Bi->TRG[0];
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

    delete[] Rr;
    delete[] Rc;
    delete[] Sr;
    delete[] Sc;
    delete[] X;
    delete[] B;
    delete[] Btrue;

  }
  return 0;
}

void elements(void * obj, int *I, int *J, dcomplex *B, int *descB) {
  if(B==NULL)
    return;
  const complex_t wavek = complex_t(10.,1.) / real_t(2 * M_PI);
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
      real_t R2=norm(dX);
      real_t R=sqrt(R2);
      B[locr*(j-1)+(i-1)]=R2==0?0.0:exp(I1*wavek*R)/R;
    }
  }

  elemops+=16*locr*locc;
}
