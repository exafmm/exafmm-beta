#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

extern "C" void FMM_Init(double eps2, int ncrit, int nworkers,
			 int nb, double * xb, double * yb, double * vb,
			 int nv, double * xv, double * yv, double * vv);
extern "C" void FMM_Finalize();
extern "C" void FMM_Partition(int & nb, double * xb, double * yb, double * vb,
                              int & nv, double * xv, double * yv, double * vv);
extern "C" void FMM_BuildTree();
extern "C" void FMM_B2B(double * vi, double * vb, int verbose);
extern "C" void FMM_V2B(double * vb, double * vv, int verbose);
extern "C" void FMM_B2V(double * vv, double * vb, int verbose);
extern "C" void FMM_V2V(double * vi, double * vv, int verbose);
extern "C" void Direct(int ni, double * xi, double * yi, double * vi,
		       int nj, double * xj, double * yj, double * vj);

void Validate(int n, double * vi, double * vd, int verbose) {
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  for (int i=0; i<n; i++) {
    diff1 += (vi[i] - vd[i]) * (vi[i] - vd[i]);
    norm1 += vd[i] * vd[i];
  }
  MPI_Reduce(&diff1, &diff2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (verbose) {
    std::cout << "--- FMM vs. direct ---------------" << std::endl;
    std::cout << std::setw(20) << std::left << std::scientific
              << "Rel. L2 Error" << " : " << std::sqrt(diff2/norm2) << std::endl;
  }
}

int main(int argc, char ** argv) {
  const int Nmax = 10000000;
  double * xb = new double [Nmax];
  double * yb = new double [Nmax];
  double * vb = new double [Nmax];
  double * xv = new double [Nmax];
  double * yv = new double [Nmax];
  double * vv = new double [Nmax];
  double * vd = new double [Nmax];
  double * vi = new double [Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int nb = 1000000 / mpisize;
  int nv = 1000000 / mpisize;
  if(mpirank < 2) nb = 0;

  srand48(mpirank);
  for (int i=0; i<nb; i++) {
    xb[i] = drand48();
    yb[i] = drand48();
  }
  for (int i=0; i<nv; i++) {
    xv[i] = drand48();
    yv[i] = drand48();
  }
  FMM_Init(0.0, 64, 16, nb, xb, yb, vb, nv, xv, yv, vv);
  FMM_Partition(nb, xb, yb, vb, nv, xv, yv, vv);
  FMM_BuildTree();

  for (int i=0; i<nb; i++) {
    vb[i] = drand48() - .5;
    vi[i] = 0;
    vd[i] = 0;
  }
  FMM_B2B(vi, vb, 1);
  Direct(100, xb, yb, vd, nb, xb, yb, vb);
  Validate(100, vi, vd, mpirank == 0);

  for (int i=0; i<nb; i++) {
    vb[i] = 0;
    vd[i] = 0;
  }
  for (int i=0; i<nv; i++) {
    vv[i] = drand48() - .5;
  }
  FMM_V2B(vb, vv, 1);
  Direct(100, xb, yb, vd, nv, xv, yv, vv);
  Validate(100, vb, vd, mpirank == 0);

  for (int i=0; i<nb; i++) {
    vb[i] = drand48() - .5;
  }
  for (int i=0; i<nv; i++) {
    vv[i] = 0;
    vd[i] = 0;
  }
  FMM_B2V(vv, vb, 1);
  Direct(100, xv, yv, vd, nb, xb, yb, vb);
  Validate(100, vv, vd, mpirank == 0);

  for (int i=0; i<nv; i++) {
    vv[i] = drand48() - .5;
    vi[i] = 0;
    vd[i] = 0;
  }
  FMM_V2V(vi, vv, 1);
  Direct(100, xv, yv, vd, nv, xv, yv, vv);
  Validate(100, vi, vd, mpirank == 0);

  FMM_Finalize();
  MPI_Finalize();
  delete[] xb;
  delete[] yb;
  delete[] vb;
  delete[] xv;
  delete[] yv;
  delete[] vv;
  delete[] vd;
  delete[] vi;
}
