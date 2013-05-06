#include <mpi.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
using boost::math::cyl_bessel_k;
using boost::math::tgamma;

extern "C" void FMM(int n, double* x, double* q, double *p, double nu, double rho);

extern "C" void MPI_Shift(double *var, int n, int mpisize, int mpirank) {
  double *buf = new double [n];
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  MPI_Request sreq, rreq;
  MPI_Isend(var, n, MPI_DOUBLE, irecv, 1, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(buf, n, MPI_DOUBLE, isend, 1, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  for (int i=0; i<n; i++) {
    var[i] = buf[i];
  }
  delete[] buf;
}

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);
  int mpisize, mpirank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  const int N = 10000;
  const double size = 2 * M_PI;
  const double nu = 1.5;
  const double rho = 10;
  double *xi = new double [3*N];
  double *qi = new double [N];
  double *pi = new double [N];
  double *pd = new double [N];
  double *xj = new double [3*N];
  double *qj = new double [N];

  srand48(mpirank);
  for (int i=0; i<N; i++) {
    xi[3*i+0] = drand48() * size - M_PI;
    xi[3*i+1] = drand48() * size - M_PI;
    xi[3*i+2] = drand48() * size - M_PI;
    qi[i] = 1. / N;
    pi[i] = 0;
    pd[i] = 0;
    xj[3*i+0] = xi[3*i+0];
    xj[3*i+1] = xi[3*i+1];
    xj[3*i+2] = xi[3*i+2];
    qj[i] = qi[i];
  }

  FMM(N, xi, qi, pi, nu, rho);
  for (int i=0; i<N; i++) {
    xj[3*i+0] = xi[3*i+0];
    xj[3*i+1] = xi[3*i+1];
    xj[3*i+2] = xi[3*i+2];
  }
  for (int irank=0; irank<mpisize; irank++) {
    if (mpirank==0) std::cout << "Direct loop          : " << irank+1 << "/" << mpisize << std::endl;
    MPI_Shift(xj, 3*N, mpisize, mpirank);
    MPI_Shift(qj, N, mpisize, mpirank);
    double coef = 1 / (std::pow(2,nu-1) * tgamma(nu));
    for (int i=0; i<100; i++) {
      double P = 0;
      for (int j=0; j<N; j++) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R = std::sqrt(2 * nu * (dx * dx + dy * dy + dz * dz)) / rho;
        if( irank == mpisize-1 && i == j ) P += qj[j];
        else P += qj[j] * std::pow(R,nu) * cyl_bessel_k(nu,R) * coef;
      }
      pd[i] += P;
    }
  }
  double Pd = 0, Pn = 0;
  for (int i=0; i<100; i++) {
    Pd += (pi[i] - pd[i]) * (pi[i] - pd[i]);
    Pn += pd[i] * pd[i];
  }
  std::cout << std::fixed << std::setprecision(7);
  std::cout << "Potential @ rank " << mpirank << "   : " << sqrtf(Pd/Pn) << std::endl;

  delete[] xi;
  delete[] qi;
  delete[] pi;
  delete[] pd;
  delete[] xj;
  delete[] qj;

  MPI_Finalize();
}
