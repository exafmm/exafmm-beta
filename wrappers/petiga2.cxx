#include "tree_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "up_down_pass.h"
#include "verify.h"

Args *args;
Logger *logger;
Bounds localBounds;
BoundBox *boundbox;
BuildTree *build;
UpDownPass *pass;
Traversal *traversal;
TreeMPI *treeMPI;

extern "C" void FMM_Init() {
  const int ncrit = 16;
  const int nspawn = 1000;
  const int images = 0;
  const real_t theta = 0.4;
  args = new Args;
  logger = new Logger;
  boundbox = new BoundBox(nspawn);
  build = new BuildTree(ncrit, nspawn);
  pass = new UpDownPass(theta);
  traversal = new Traversal(nspawn, images);
  treeMPI = new TreeMPI(images);
  args->numBodies /= treeMPI->mpisize;
  args->verbose &= treeMPI->mpirank == 0;
  if (args->verbose) {
    logger->verbose = true;
    boundbox->verbose = true;
    build->verbose = true;
    pass->verbose = true;
    traversal->verbose = true;
    treeMPI->verbose = true;
  }
  logger->printTitle("Initial Parameters");
  args->print(logger->stringLength, P, treeMPI->mpirank);
}

extern "C" void FMM_Finalize() {
  delete args;
  delete logger;
  delete boundbox;
  delete build;
  delete pass;
  delete traversal;
  delete treeMPI;
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

extern "C" void Direct_Laplace(int ni, double * xi, double * yi, double * zi, double * vi,
                               int nj, double * xj, double * yj, double * zj, double * vj) {
  const int Nmax = 2000000;
  double * x2 = new double [Nmax];
  double * y2 = new double [Nmax];
  double * z2 = new double [Nmax];
  double * v2 = new double [Nmax];
  for (int i=0; i<nj; i++) {
    x2[i] = xj[i];
    y2[i] = yj[i];
    z2[i] = zj[i];
    v2[i] = vj[i];
  }
  if (treeMPI->mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
  for (int irank=0; irank<treeMPI->mpisize; irank++) {
    if (treeMPI->mpirank == 0) std::cout << "Direct loop          : " << irank+1 << "/" << treeMPI->mpisize << std::endl;
    int n2 = nj;
    MPI_Shift(x2, nj, treeMPI->mpisize, treeMPI->mpirank);
    nj = n2;
    MPI_Shift(y2, nj, treeMPI->mpisize, treeMPI->mpirank);
    nj = n2;
    MPI_Shift(z2, nj, treeMPI->mpisize, treeMPI->mpirank);
    nj = n2;
    MPI_Shift(v2, nj, treeMPI->mpisize, treeMPI->mpirank);
    for (int i=0; i<ni; i++) {
      double pp = 0;
      for (int j=0; j<nj; j++) {
        double dx = xi[i] - x2[j];
        double dy = yi[i] - y2[j];
        double dz = zi[i] - z2[j];
        double R2 = dx * dx + dy * dy + dz * dz;
        double invR = 1 / std::sqrt(R2);
        if (R2 == 0) invR = 0;
        pp += v2[j] * invR;
      }
      vi[i] += pp;
    }
  }
  delete[] x2;
  delete[] y2;
  delete[] z2;
  delete[] v2;
}

int main(int argc, char ** argv) {
  Dataset data;
  Verify verify;
  const int Nmax = 2000000;
  int ni = 125000;
  int nj = 125000;
  int stringLength = 20;
  double * xi = new double [Nmax];
  double * yi = new double [Nmax];
  double * zi = new double [Nmax];
  double * vi = new double [Nmax];
  double * xj = new double [Nmax];
  double * yj = new double [Nmax];
  double * zj = new double [Nmax];
  double * vj = new double [Nmax];
  double * v2 = new double [Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  srand48(mpirank);
  for (int i=0; i<ni; i++) {
    xi[i] = drand48() - .5;
    yi[i] = drand48() - .5;
    zi[i] = drand48() - .5;
    vi[i] = 0;
  }
  for (int i=0; i<nj; i++) {
    xj[i] = drand48() - .5;
    yj[i] = drand48() - .5;
    zj[i] = drand48() - .5;
    vj[i] = drand48() - .5;
  }

  FMM_Init();
  verify.verbose = treeMPI->mpirank == 0;
  Bodies bodies(ni);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0] = xi[i];
    B->X[1] = yi[i];
    B->X[2] = zi[i];
    B->SRC  = vi[i];
  }
  Bodies jbodies(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0] = xj[i];
    B->X[1] = yj[i];
    B->X[2] = zj[i];
    B->SRC  = vj[i];
  }
  localBounds = boundbox->getBounds(bodies);
  localBounds = boundbox->getBounds(jbodies,localBounds);
  Bounds globalBounds = treeMPI->allreduceBounds(localBounds);
  localBounds = treeMPI->partition(bodies,globalBounds);
  bodies = treeMPI->commBodies(bodies);
  treeMPI->partition(jbodies,globalBounds);
  jbodies = treeMPI->commBodies(jbodies);
  Cells cells = build->buildTree(bodies, localBounds);
  pass->upwardPass(cells);
  Cells jcells = build->buildTree(jbodies, localBounds);
  pass->upwardPass(jcells);

  ni = bodies.size();
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    xi[i] = B->X[0];
    yi[i] = B->X[1];
    zi[i] = B->X[2];
    vi[i] = B->SRC;
  }
  nj = jbodies.size();
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    xj[i] = B->X[0];
    yj[i] = B->X[1];
    zj[i] = B->X[2];
    vj[i] = B->SRC;
  }

  logger->printTitle("FMM Parameters");
  args->print(logger->stringLength, P, treeMPI->mpirank);
  logger->printTitle("FMM Profiling");
  logger->startTimer("Total FMM");
  logger->startPAPI();
  const real_t cycle = 2 * M_PI;
  bodies.resize(ni);
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    B->X[0]   = xi[i];
    B->X[1]   = yi[i];
    B->X[2]   = zi[i];
    B->SRC    = 1;
    B->TRG[0] = vi[i];
    B->TRG[1] = 0;
    B->TRG[2] = 0;
    B->TRG[3] = 0;
    B->IBODY = i;
  }
  jbodies.resize(nj);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    int i = B-jbodies.begin();
    B->X[0]   = xj[i];
    B->X[1]   = yj[i];
    B->X[2]   = zj[i];
    B->SRC    = vj[i];
    B->TRG    = 0;
    B->IBODY = i;
  }
  cells = build->buildTree(bodies, localBounds);
  pass->upwardPass(cells);
  jcells = build->buildTree(jbodies, localBounds);
  pass->upwardPass(jcells);
  treeMPI->setLET(jcells,localBounds,cycle);
  treeMPI->commBodies();
  treeMPI->commCells();
  traversal->dualTreeTraversal(cells, jcells, cycle);
  for (int irank=1; irank<treeMPI->mpisize; irank++) {
    treeMPI->getLET(jcells,(treeMPI->mpirank+irank)%treeMPI->mpisize);
    traversal->dualTreeTraversal(cells, jcells, cycle);
  }
  pass->downwardPass(cells);

  logger->stopPAPI();
  logger->stopTimer("Total FMM");
  logger->printTitle("MPI direct sum");
  data.sampleBodies(bodies, args->numTargets);
  Bodies bodies2 = bodies;
  data.initTarget(bodies);
  logger->startTimer("Total Direct");
  for (int i=0; i<treeMPI->mpisize; i++) {
    if (args->verbose) std::cout << "Direct loop          : " << i+1 << "/" << treeMPI->mpisize << std::endl;
    treeMPI->shiftBodies(jbodies);
    traversal->direct(bodies, jbodies, cycle);
  }
  traversal->normalize(bodies);
  logger->printTitle("Total runtime");
  logger->printTime("Total FMM");
  logger->stopTimer("Total Direct");
  boundbox->writeTime(treeMPI->mpirank);
  build->writeTime(treeMPI->mpirank);
  pass->writeTime(treeMPI->mpirank);
  traversal->writeTime(treeMPI->mpirank);
  treeMPI->writeTime(treeMPI->mpirank);
  double potDif = verify.getDifScalar(bodies, bodies2);
  double potNrm = verify.getNrmScalar(bodies);
  double accDif = verify.getDifVector(bodies, bodies2);
  double accNrm = verify.getNrmVector(bodies);
  double potDifGlob, potNrmGlob, accDifGlob, accNrmGlob;
  MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accDif, &accDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accNrm, &accNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  logger->printTitle("FMM vs. direct");
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDifGlob/potNrmGlob));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDifGlob/accNrmGlob));
  build->printTreeData(cells);
  traversal->printTraversalData();
  logger->printPAPI();

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    int i = B-bodies.begin();
    xi[i] = B->X[0];
    yi[i] = B->X[1];
    zi[i] = B->X[2];
    vi[i] = B->TRG[0];
  }
  ni = bodies.size();
  for (int i=0; i<ni; i++) {
    v2[i] = 0;
  }
  Direct_Laplace(ni, xi, yi, zi, v2, nj, xj, yj, zj, vj);
  potDif = potNrm = 0;
  for (int i=0; i<ni; i++) {
    potDif += (vi[i] - v2[i]) * (vi[i] - v2[i]);
    potNrm += v2[i] * v2[i];
  }
  //double potDifGlob, potNrmGlob;
  MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpirank == 0) {
    std::cout << "--- FMM vs. Direct ---------------" << std::endl;
    std::cout << std::setw(stringLength) << std::left << std::scientific
              << "Rel. L2 Error (pot)" << " : " << std::sqrt(potDifGlob/potNrmGlob) << std::endl;
  }


  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] vi;
  delete[] xj;
  delete[] yj;
  delete[] zj;
  delete[] vj;
  delete[] v2;
  FMM_Finalize();
  MPI_Finalize();
  return 0;
}
