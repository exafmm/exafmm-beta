#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

extern "C" void FMM_Init(int images, int threads, double theta, double cutoff, int verbose);
extern "C" void FMM_Finalize();
extern "C" void FMM_Partition(int & ni, int nimax, int * res_index, double * x, double * q, double *v, double * cycle);
extern "C" void FMM_FMM(int ni, int &nj, double * x, double * q, double * p, double * f, double * cycle);
extern "C" void FMM_Ewald(int ni, double * x, double * q, double * p, double * f,
			  int ksize, double alpha, double sigma, double cutoff, double * cycle);
extern "C" void FMM_Cutoff(int ni, double * x, double * q, double * p, double * f, double cutoff, double * cycle);

int main(int argc, char ** argv) {
  const int nimax = 1000000;
  int ni = 1000;
  int nj;
  int stringLength = 20;
  int images = 3;
  int ksize = 11;
  int threads = 16;
  int verbose = 1;
  double theta = 0.5;
  double cycle[3] = {2*M_PI, 2*M_PI, 2*M_PI};
  double alpha = 10 / cycle[0];
  double sigma = .25 / M_PI;
  double cutoff = cycle[0] / 2;
  int * res_index = new int [nimax];
  double * x = new double [3*nimax];
  double * q = new double [nimax];
  double * v = new double [3*nimax];
  double * p = new double [nimax];
  double * f = new double [3*nimax];
  double * p2 = new double [nimax];
  double * f2 = new double [3*nimax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  srand48(mpirank);
  double average = 0;
  for (int i=0; i<ni; i++) {
    x[3*i+0] = drand48() * cycle[0] - cycle[0] / 2;
    x[3*i+1] = drand48() * cycle[1] - cycle[1] / 2;
    x[3*i+2] = drand48() * cycle[2] - cycle[2] / 2;
    p[i] = f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
    res_index[i] = i + mpirank*ni;
  }
  for (int i=0; i<ni; i++) {
    q[i] = drand48() - .5;
    average += q[i];
  }
  average /= ni;
  for (int i=0; i<ni; i++) {
    q[i] -= average;
  }

  FMM_Init(images, threads, theta, cutoff, verbose);
  FMM_Partition(ni, nimax, res_index, x, q, v, cycle);
  FMM_FMM(ni, nj, x, q, p, f, cycle);
  for (int i=0; i<ni; i++) {
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
  FMM_Ewald(ni, x, q, p2, f2, ksize, alpha, sigma, cutoff, cycle);
  double potSum = 0, potSum2 = 0, accDif = 0, accNrm = 0;
  for (int i=0; i<ni; i++) {
    potSum  += p[i]  * q[i];
    potSum2 += p2[i] * q[i];
    accDif  += (f[3*i+0] - f2[3*i+0]) * (f[3*i+0] - f2[3*i+0])
      + (f[3*i+1] - f2[3*i+1]) * (f[3*i+1] - f2[3*i+1])
      + (f[3*i+2] - f2[3*i+2]) * (f[3*i+2] - f2[3*i+2]);
    accNrm  += f2[3*i+0] * f2[3*i+0] + f2[3*i+1] * f2[3*i+1] + f2[3*i+2] * f2[3*i+2];
  }
  double potSumGlob, potSumGlob2, accDifGlob, accNrmGlob;
  MPI_Reduce(&potSum,  &potSumGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potSum2, &potSumGlob2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accDif,  &accDifGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accNrm,  &accNrmGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  double potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
  double potNrmGlob = potSumGlob * potSumGlob;
  if (mpirank == 0) {
    std::cout << "--- FMM vs. Ewald  ---------------" << std::endl;
    std::cout << std::setw(stringLength) << std::left << std::scientific
  	      << "Rel. L2 Error (pot)" << " : " << std::sqrt(potDifGlob/potNrmGlob) << std::endl;
    std::cout << std::setw(stringLength) << std::left
	      << "Rel. L2 Error (acc)" << " : " << std::sqrt(accDifGlob/accNrmGlob) << std::endl;
  }
#if 0
  for (int i=0; i<ni; i++) {
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
  FMM_Cutoff(ni, x, q, p2, f2, cutoff, cycle);
  potSum = potSum2 = accDif = accNrm = 0;
  for (int i=0; i<ni; i++) {
    potSum  += p[i]  * q[i];
    potSum2 += p2[i] * q[i];
    accDif  += (f[3*i+0] - f2[3*i+0]) * (f[3*i+0] - f2[3*i+0])
      + (f[3*i+1] - f2[3*i+1]) * (f[3*i+1] - f2[3*i+1])
      + (f[3*i+2] - f2[3*i+2]) * (f[3*i+2] - f2[3*i+2]);
    accNrm  += f2[3*i+0] * f2[3*i+0] + f2[3*i+1] * f2[3*i+1] + f2[3*i+2] * f2[3*i+2];
  }
  MPI_Reduce(&potSum,  &potSumGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potSum2, &potSumGlob2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accDif,  &accDifGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&accNrm,  &accNrmGlob,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  potDifGlob = (potSumGlob - potSumGlob2) * (potSumGlob - potSumGlob2);
  potNrmGlob = potSumGlob * potSumGlob;
  if (mpirank == 0) {
    std::cout << "--- FMM vs. Cutoff  ---------------" << std::endl;
    std::cout << std::setw(stringLength) << std::left << std::scientific
              << "Rel. L2 Error (pot)" << " : " << std::sqrt(potDifGlob/potNrmGlob) << std::endl;
    std::cout << std::setw(stringLength) << std::left
              << "Rel. L2 Error (acc)" << " : " << std::sqrt(accDifGlob/accNrmGlob) << std::endl;
  }
#endif

  delete[] res_index;
  delete[] x;
  delete[] q;
  delete[] v;
  delete[] p;
  delete[] f;
  delete[] p2;
  delete[] f2;
  FMM_Finalize();
  MPI_Finalize();
}
