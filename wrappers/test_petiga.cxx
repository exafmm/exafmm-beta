#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

extern "C" void FMM_Init();
extern "C" void FMM_Finalize();
extern "C" void FMM_Partition(int & ni, double * xi, double * yi, double * zi, double * vi,
			      int & nj, double * xj, double * yj, double * zj, double * vj);
extern "C" void FMM_Laplace(int ni, double * xi, double * yi, double * zi, double * vi,
			    int nj, double * xj, double * yj, double * zj, double * vj);
extern "C" void Direct_Laplace(int ni, double * xi, double * yi, double * zi, double * vi,
			       int nj, double * xj, double * yj, double * zj, double * vj);

int main(int argc, char ** argv) {
  const int Nmax = 1000000;
  int ni = 500;
  int nj = 1000;
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
    vj[i] = 1. / nj;
  }

  FMM_Init();
  FMM_Partition(ni, xi, yi, zi, vi, nj, xj, yj, zj, vj);
  FMM_Laplace(ni, xi, yi, zi, vi, nj, xj, yj, zj, vj);
  for (int i=0; i<ni; i++) {
    v2[i] = 0;
  }
  Direct_Laplace(ni, xi, yi, zi, v2, nj, xj, yj, zj, vj);
  double potDif = 0, potNrm = 0;
  for (int i=0; i<ni; i++) {
    potDif += (vi[i] - v2[i]) * (vi[i] - v2[i]);
    potNrm += v2[i] * v2[i];
  }
  double potDifGlob, potNrmGlob;
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
}
