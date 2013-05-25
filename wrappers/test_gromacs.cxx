#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

extern "C" void fmm(int n, double * x, double * q, double * p, double * f, double cycle, int images);
extern "C" void ewald(int n, double * x, double * q, double * p, double * f, int ksize, double alpha, double cycle);

extern "C" void MPI_Shift(double * var, int n, int mpisize, int mpirank) {
  double * buf = new double [n];
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

int main() {
  const int N = 1000;
  const int stringLength = 20;
  const int images = 3;
  const int ksize = 11;
  const double cycle = 2 * M_PI;
  const double alpha = 10 / cycle;
  double *x = new double [3*N];
  double *q = new double [N];
  double *p = new double [N];
  double *f = new double [3*N];
  double *x2 = new double [3*N];
  double *q2 = new double [N];
  double *p2 = new double [N];
  double *f2 = new double [3*N];

  srand48(0);
  double average = 0;
  for (int i=0; i<N; i++) {
    x[3*i+0] = drand48() * cycle - cycle / 2;
    x[3*i+1] = drand48() * cycle - cycle / 2;
    x[3*i+2] = drand48() * cycle - cycle / 2;
    x2[3*i+0] = x[3*i+0];
    x2[3*i+1] = x[3*i+1];
    x2[3*i+2] = x[3*i+2];
    p[i] = f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
  for (int i=0; i<N; i++) {
    q[i] = drand48() - .5;
    average += q[i];
  }
  average /= N;
  for (int i=0; i<N; i++) {
    q[i] -= average;
    q2[i] = q[i];
  }

  fmm(N, x, q, p, f, cycle, images);
#if 1
  ewald(N, x2, q2, p2, f2, ksize, alpha, cycle);
#else
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  double dipole[3] = {0, 0, 0};
  for (int i=0; i<N; i++) {
    for (int d=0; d<3; d++) dipole[d] += x[3*i+d] * q[i];
  }
  double norm = 0;
  for (int d=0; d<3; d++) {
    norm += dipole[d] * dipole[d];
  }
  double coef = 4 * M_PI / (3 * cycle * cycle * cycle);
  double Xperiodic[3];
  for (int i=0; i<N; i++) {
    double P = 0, Fx = 0, Fy = 0, Fz = 0;
    for (int ix=-prange; ix<=prange; ix++) {
      for (int iy=-prange; iy<=prange; iy++) {
        for (int iz=-prange; iz<=prange; iz++) {
          Xperiodic[0] = ix * cycle;
          Xperiodic[1] = iy * cycle;
          Xperiodic[2] = iz * cycle;
	  for (int j=0; j<N; j++) {
	    double dx = x[3*i+0] - x2[3*j+0] - Xperiodic[0];
	    double dy = x[3*i+1] - x2[3*j+1] - Xperiodic[1];
	    double dz = x[3*i+2] - x2[3*j+2] - Xperiodic[2];
	    double R2 = dx * dx + dy * dy + dz * dz;
	    double invR = 1 / std::sqrt(R2);
	    if (R2 == 0) invR = 0;
	    double invR3 = q2[j] * invR * invR * invR;
	    P += q2[j] * invR;
	    Fx += dx * invR3;
	    Fy += dy * invR3;
	    Fz += dz * invR3;
	  }
	}
      }
    }
    p2[i] += P - coef * norm / N / q2[i];
    f2[3*i+0] -= Fx + coef * dipole[0];
    f2[3*i+1] -= Fy + coef * dipole[1];
    f2[3*i+2] -= Fz + coef * dipole[2];
  }
#endif
  double diff2 = 0, norm2 = 0;
  double v = 0, v2 = 0;
  for (int i=0; i<N; i++) {
    v += p[i] * q[i];
    v2 += p2[i] * q2[i];
    diff2 += (f[3*i+0] - f2[3*i+0]) * (f[3*i+0] - f2[3*i+0])
           + (f[3*i+1] - f2[3*i+1]) * (f[3*i+1] - f2[3*i+1])
           + (f[3*i+2] - f2[3*i+2]) * (f[3*i+2] - f2[3*i+2]);
    norm2 += f2[3*i+0] * f2[3*i+0] + f2[3*i+1] * f2[3*i+1] + f2[3*i+2] * f2[3*i+2];
  }
  double diff1 = (v - v2) * (v - v2);
  double norm1 = v2 * v2;
  std::cout << "--- FMM vs. Ewald  ---------------" << std::endl;
  std::cout << std::setw(stringLength) << std::left
	    << "Rel. L2 Error (pot)" << " : " << std::sqrt(diff1/norm1) << std::endl;
  if( std::abs(diff2) > 0 ) {
    std::cout << std::setw(stringLength) << std::left
	      << "Rel. L2 Error (acc)" << " : " << std::sqrt(diff2/norm2) << std::endl;
  }

  delete[] x;
  delete[] q;
  delete[] p;
  delete[] f;
  delete[] x2;
  delete[] q2;
  delete[] p2;
  delete[] f2;
}
