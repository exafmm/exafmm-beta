#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

extern "C" void FMM_Init(int images, int threads, int verbose);
extern "C" void FMM_Finalize();
extern "C" void FMM_Partition(int & n, int * ibody, int * icell, float * x, float * q, float cycle);
extern "C" void FMM_Coulomb(int n, int * icell, float * x, float * q, float * p, float * f, float cycle);
extern "C" void Ewald_Coulomb(int n, float * x, float * q, float * p, float * f,
			      int ksize, float alpha, float sigma, float cutoff, float cycle);
extern "C" void Direct_Coulomb(int n, float * x, float * q, float * p, float * f, float cycle);

int main(int argc, char ** argv) {
  const int Nmax = 1000000;
  int Ni = 500;
  int stringLength = 20;
  int images = 3;
  int ksize = 11;
  int threads = 16;
  int verbose = 1;
  float cycle = 2 * M_PI;
  float alpha = 10 / cycle;
  float sigma = .25 / M_PI;
  float cutoff = cycle / 2;
  int * ibody = new int [Nmax];
  int * icell = new int [Nmax];
  float * x = new float [3*Nmax];
  float * q = new float [Nmax];
  float * p = new float [Nmax];
  float * f = new float [3*Nmax];
  float * p2 = new float [Nmax];
  float * f2 = new float [3*Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

#if 1
  srand48(mpirank);
  double average = 0;
  for (int i=0; i<Ni; i++) {
    x[3*i+0] = drand48() * cycle - cycle / 2;
    x[3*i+1] = drand48() * cycle - cycle / 2;
    x[3*i+2] = drand48() * cycle - cycle / 2;
    p[i] = f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
    ibody[i] = i + mpirank*Ni;
    icell[i] = ibody[i];
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
    ibody[i] = i + mpirank*Ni;
    icell[i] = ibody[i];
  }
  file.close();
#endif

#if EXAFMM_CLUSTER // Use preassigned cells
  int ncrit = 32;
  int level = Ni >= ncrit ? 1 + int(log2(Ni / ncrit)/3) : 0;
  float diameter = cycle / (1 << level);
  for (int i=0; i<Ni; i++) {
    int iX[3] = {0, 0, 0};
    for (int d=0; d<3; d++) iX[d] = int((x[3*i+d] + cycle / 2) / diameter);
    int key = 0;
    for (int l=0; l<level; l++) {
      for (int d=0; d<3; d++) key += (iX[d] & 1) << (3 * l + d);
      for (int d=0; d<3; d++) iX[d] >>= 1;
    }
    icell[i] = key;
  }
  int numBucket = 1 << (3 * level);
  int *bucket = new int [numBucket];
  for( int i=0; i<numBucket; i++ ) bucket[i] = 0;
  for( int i=0; i<Ni; i++ ) bucket[icell[i]]++;
  for( int i=1; i<numBucket; i++ ) bucket[i] += bucket[i-1];
  for( int i=Ni-1; i>=0; --i ) {
    bucket[icell[i]]--;
    int inew = bucket[icell[i]];
    p2[inew] = q[i];
    for (int d=0; d<3; d++) f2[3*inew+d] = x[3*i+d];
  }
  for (int i=0; i<Ni; i++) {
    q[i] = p2[i];
    for (int d=0; d<3; d++) x[3*i+d] = f2[3*i+d];
  }
  delete[] bucket;
  for (int i=0; i<Ni; i++) {
    int iX[3] = {0, 0, 0};
    for (int d=0; d<3; d++) iX[d] = int((x[3*i+d] + cycle / 2) / diameter);
    int key = 0;
    for (int l=0; l<level; l++) {
      for (int d=0; d<3; d++) key += (iX[d] & 1) << (3 * l + d);
      for (int d=0; d<3; d++) iX[d] >>= 1;
    }
    icell[i] = key;
  }
#endif
  FMM_Init(images, threads, verbose);
  FMM_Partition(Ni, ibody, icell, x, q, cycle);
  FMM_Coulomb(Ni, icell, x, q, p, f, cycle);
  for (int i=0; i<Ni; i++) {
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
#if 1
  Ewald_Coulomb(Ni, x, q, p2, f2, ksize, alpha, sigma, cutoff, cycle);
#else
  Direct_Coulomb(Ni, x, q, p2, f2, cycle);
#endif
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

  delete[] ibody;
  delete[] icell;
  delete[] x;
  delete[] q;
  delete[] p;
  delete[] f;
  delete[] p2;
  delete[] f2;
  FMM_Finalize();
  MPI_Finalize();
}
