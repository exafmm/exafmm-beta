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
extern "C" void FMM_Partition(int & ni, double * xi, double * yi, double * zi, double * vi,
			      int & nj, double * xj, double * yj, double * zj, double * vj);
extern "C" void FMM_BuildTree();
extern "C" void FMM_V2B(double * vi, double * vj, bool verbose);
extern "C" void Direct(int ni, double * xi, double * yi, double * zi, double * vi,
		       int nj, double * xj, double * yj, double * zj, double * vj);

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
  const int Nmax = 1000000;
  const int ncrit = 16;
  const int threads = 16;
  const double eps2 = 0.0;
  int ni = 500;
  int nj = 1000;
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
    vj[i] = 1. / nj;
  }

  FMM_Init(eps2, ncrit, threads, ni, xi, yi, zi, vi, nj, xj, yj, zj, vj);
  FMM_Partition(ni, xi, yi, zi, vi, nj, xj, yj, zj, vj);
  FMM_BuildTree();
  FMM_V2B(vi, vj, true);
  for (int i=0; i<ni; i++) {
    v2[i] = 0;
  }
  Direct(ni, xi, yi, zi, v2, nj, xj, yj, zj, vj);
  Validate(100, vi, v2, mpirank == 0);

  FMM_Finalize();
  MPI_Finalize();
  delete[] xi;
  delete[] yi;
  delete[] zi;
  delete[] vi;
  delete[] xj;
  delete[] yj;
  delete[] zj;
  delete[] vj;
  delete[] v2;
}
