#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

extern "C" void FMM(int ni, double * xi, double * pi, double * fi, int nj, double * xj, double * qj, double size, int periodicflag);

extern "C" void MPI_Shift(double * var, int n, int mpisize, int mpirank) {
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
  const int N = 1000000 / mpisize;
  const int size = 2 * M_PI;
  const int stringLength = 20;
  double *xi = new double [3*N];
  double *pi = new double [N];
  double *fi = new double [3*N];
  double *pd = new double [N];
  double *fd = new double [3*N];
  double *xj = new double [3*N];
  double *qj = new double [N];

  srand48(mpirank);
  double average = 0;
  for (int i=0; i<N; i++) {
    xi[3*i+0] = drand48() * size - size / 2;
    xi[3*i+1] = drand48() * size - size / 2;
    xi[3*i+2] = drand48() * size - size / 2;
    pi[i] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[i] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
    xj[3*i+0] = drand48() * size - size / 2;
    xj[3*i+1] = drand48() * size - size / 2;
    xj[3*i+2] = drand48() * size - size / 2;
    qj[i] = drand48() - .5;
    average += qj[i];
  }
  average /= N;
  for (int i=0; i<N; i++) {
    qj[i] -= average;
  }

  FMM(N, xi, pi, fi, N, xj, qj, size, 0);
  if (mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
  for (int irank=0; irank<mpisize; irank++) {
    if (mpirank==0) std::cout << "Direct loop          : " << irank+1 << "/" << mpisize << std::endl;
    MPI_Shift(xj, 3*N, mpisize, mpirank);
    MPI_Shift(qj, N, mpisize, mpirank);
    for (int i=0; i<100; i++) {
      double P = 0, Fx = 0, Fy = 0, Fz = 0;
      for (int j=0; j<N; j++) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        double invR = 1 / std::sqrt(R2);
        if( irank == mpisize-1 && i == j ) invR = 0;
        double invR3 = qj[j] * invR * invR * invR;
        P += qj[j] * invR;
        Fx += dx * invR3;
        Fy += dy * invR3;
        Fz += dz * invR3;
      }
      pd[i] += P;
      fd[3*i+0] -= Fx;
      fd[3*i+1] -= Fy;
      fd[3*i+2] -= Fz;
    }
  }
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  for (int i=0; i<100; i++) {
    diff1 += (pi[i] - pd[i]) * (pi[i] - pd[i]);
    norm1 += pd[i] * pd[i];
    diff2 += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
          + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
          + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    norm2 += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
  }
  MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpirank == 0) {
    std::cout << "--- FMM vs. direct ---------------" << std::endl;
    std::cout << std::setw(stringLength) << std::left
	      << "Rel. L2 Error (pot)" << " : " << std::sqrt(diff3/norm3) << std::endl;
    if( std::abs(diff3) > 0 ) {
      std::cout << std::setw(stringLength) << std::left
	        << "Rel. L2 Error (acc)" << " : " << std::sqrt(diff4/norm4) << std::endl;
    }
  }    

  delete[] xi;
  delete[] pi;
  delete[] fi;
  delete[] pd;
  delete[] fd;
  delete[] xj;
  delete[] qj;

  MPI_Finalize();
}
