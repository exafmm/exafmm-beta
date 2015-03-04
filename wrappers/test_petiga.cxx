#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

extern "C" void FMM_Init(double eps2, int ncrit, int threads,
                         int nb, double * xb, double * yb, double * zb, double * vb,
                         int nv, double * xv, double * yv, double * zv, double * vv);
extern "C" void FMM_Finalize();
extern "C" void FMM_Partition(int & nb, double * xb, double * yb, double * zb, double * vb,
			      int & nv, double * xv, double * yv, double * zv, double * vv);
extern "C" void FMM_BuildTree();
extern "C" void FMM_B2B(double * vi, double * vb, bool verbose);
extern "C" void FMM_V2B(double * vb, double * vv, bool verbose);
extern "C" void FMM_B2V(double * vv, double * vb, bool verbose);
extern "C" void FMM_V2V(double * vi, double * vv, bool verbose);
extern "C" void Direct(int nb, double * xb, double * yb, double * zb, double * vb,
		       int nv, double * xv, double * yv, double * zv, double * vv);

void Validate(int n, double * vb, double * vd, int verbose) {
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  for (int i=0; i<n; i++) {
    diff1 += (vb[i] - vd[i]) * (vb[i] - vd[i]);
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
  const int ncrit = 16;
  const int threads = 16;
  const double eps2 = 0.0;
  double * xb = new double [Nmax];
  double * yb = new double [Nmax];
  double * zb = new double [Nmax];
  double * vb = new double [Nmax];
  double * xv = new double [Nmax];
  double * yv = new double [Nmax];
  double * zv = new double [Nmax];
  double * vv = new double [Nmax];
  double * vd = new double [Nmax];
  double * vi = new double [Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int nb = 50000 / mpisize;
  int nv = 100000 / mpisize;

  srand48(mpirank);
  for (int i=0; i<nb; i++) {
    xb[i] = drand48() - .5;
    yb[i] = drand48() - .5;
    zb[i] = drand48() - .5;
  }
  for (int i=0; i<nv; i++) {
    xv[i] = drand48() - .5;
    yv[i] = drand48() - .5;
    zv[i] = drand48() - .5;
  }

  FMM_Init(eps2, ncrit, threads, nb, xb, yb, zb, vb, nv, xv, yv, zv, vv);
  FMM_Partition(nb, xb, yb, zb, vb, nv, xv, yv, zv, vv);
  FMM_BuildTree();

  for (int i=0; i<nb; i++) {
    vb[i] = 1.0 / nb;
    vi[i] = 0;
    vd[i] = 0;
  }
  FMM_B2B(vi, vb, 1);
  Direct(100, xb, yb, zb, vd, nb, xb, yb, zb, vb);
  Validate(100, vi, vd, mpirank == 0);

  for (int i=0; i<nb; i++) {
    vb[i] = 0;
    vd[i] = 0;
  }
  for (int i=0; i<nv; i++) {
    vv[i] = 1.0 / nv;
  }
  FMM_V2B(vb, vv, true);
  Direct(100, xb, yb, zb, vd, nv, xv, yv, zv, vv);
  Validate(100, vb, vd, mpirank == 0);

  for (int i=0; i<nb; i++) {
    vb[i] = 1.0 / nb;
  }
  for (int i=0; i<nv; i++) {
    vv[i] = 0;
    vd[i] = 0;
  }
  FMM_B2V(vv, vb, 1);
  Direct(100, xv, yv, zv, vd, nb, xb, yb, zb, vb);
  Validate(100, vv, vd, mpirank == 0);

  for (int i=0; i<nv; i++) {
    vv[i] = 1.0 / nv;
    vi[i] = 0;
    vd[i] = 0;
  }
  FMM_V2V(vi, vv, 1);
  Direct(100, xv, yv, zv, vd, nv, xv, yv, zv, vv);
  Validate(100, vi, vd, mpirank == 0);

  FMM_Finalize();
  MPI_Finalize();
  delete[] xb;
  delete[] yb;
  delete[] zb;
  delete[] vb;
  delete[] xv;
  delete[] yv;
  delete[] zv;
  delete[] vv;
  delete[] vd;
  delete[] vi;
}
