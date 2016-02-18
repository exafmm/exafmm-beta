#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

extern "C" void FMM_Init(int images, int threads, double theta, int verbose);
extern "C" void FMM_Finalize();
extern "C" void Partition(int & n, int * res_index, double * x, double * q, double * cycle);
extern "C" void FMM(int n, double * x, double * q, double * p, double * f, double * cycle);
extern "C" void FMM_Ewald(int n, double * x, double * q, double * p, double * f,
			  int ksize, double alpha, double sigma, double cutoff, double * cycle);
extern "C" void FMM_Cutoff(int n, double * x, double * q, double * p, double * f, double cutoff, double * cycle);

int main(int argc, char ** argv) {
  const int Nmax = 1000000;
  int Ni = 1000;
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
  int * res_index = new int [Nmax];
  double * x = new double [3*Nmax];
  double * q = new double [Nmax];
  double * p = new double [Nmax];
  double * f = new double [3*Nmax];
  double * p2 = new double [Nmax];
  double * f2 = new double [3*Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

#if 1
  srand48(mpirank);
  double average = 0;
  for (int i=0; i<Ni; i++) {
    x[3*i+0] = drand48() * cycle[0] - cycle[0] / 2;
    x[3*i+1] = drand48() * cycle[1] - cycle[1] / 2;
    x[3*i+2] = drand48() * cycle[2] - cycle[2] / 2;
    p[i] = f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
    res_index[i] = i + mpirank*Ni;
  }
  for (int i=0; i<Ni; i++) {
    q[i] = drand48() - .5;
    average += q[i];
  }
  average /= Ni;
  for (int i=0; i<Ni; i++) {
    q[i] -= average;
  }
#else
  std::stringstream name;
  name << "source" << std::setfill('0') << std::setw(4)
       << mpirank << ".dat";
  std::ifstream file(name.str().c_str(),std::ios::in);
  for (int i=0; i<Ni; i++) {
    file >> x[3*i+0];
    file >> x[3*i+1];
    file >> x[3*i+2];
    file >> q[i];
    res_index[i] = i + mpirank*Ni;
  }
  file.close();
#endif

  FMM_Init(images, threads, theta, verbose);
  Partition(Ni, res_index, x, q, cycle);
  FMM(Ni, x, q, p, f, cycle);
  for (int i=0; i<Ni; i++) {
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
  FMM_Ewald(Ni, x, q, p2, f2, ksize, alpha, sigma, cutoff, cycle);
  double potSum = 0, potSum2 = 0, accDif = 0, accNrm = 0;
  for (int i=0; i<Ni; i++) {
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
  for (int i=0; i<Ni; i++) {
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
  FMM_Cutoff(Ni, x, q, p2, f2, cutoff, cycle);
  potSum = potSum2 = accDif = accNrm = 0;
  for (int i=0; i<Ni; i++) {
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
  delete[] p;
  delete[] f;
  delete[] p2;
  delete[] f2;
  FMM_Finalize();
  MPI_Finalize();
}
