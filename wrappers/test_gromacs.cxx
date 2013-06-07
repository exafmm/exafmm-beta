#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

extern "C" void FMM_Init(int images);
extern "C" void FMM_Partition(int & n, double * x, double * q, double cycle);
extern "C" void FMM(int n, double * x, double * q, double * p, double * f, double cycle);
extern "C" void FMM_Ewald(int n, double * x, double * q, double * p, double * f, int ksize, double alpha, double cycle);

extern "C" void MPI_Shift(double * var, int &nold, int mpisize, int mpirank) {
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
  int nnew;
  MPI_Request sreq, rreq;

  MPI_Isend(&nold, 1, MPI_DOUBLE, irecv, 0, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(&nnew, 1, MPI_DOUBLE, isend, 0, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);

  double * buf = new double [nnew];
  MPI_Isend(var, nold, MPI_DOUBLE, irecv, 1, MPI_COMM_WORLD, &sreq);
  MPI_Irecv(buf, nnew, MPI_DOUBLE, isend, 1, MPI_COMM_WORLD, &rreq);
  MPI_Wait(&sreq, MPI_STATUS_IGNORE);
  MPI_Wait(&rreq, MPI_STATUS_IGNORE);
  for (int i=0; i<nnew; i++) {
    var[i] = buf[i];
  }
  nold = nnew;
  delete[] buf;
}

int main(int argc, char ** argv) {
  const int Nmax = 1000000;
  int Ni = 1000;
  int stringLength = 20;
  int images = 0;
  int ksize = 11;
  double cycle = 2 * M_PI;
  double alpha = 10 / cycle;
  double *x = new double [3*Nmax];
  double *q = new double [Nmax];
  double *p = new double [Nmax];
  double *f = new double [3*Nmax];
  double *x2 = new double [3*Nmax];
  double *q2 = new double [Nmax];
  double *p2 = new double [Nmax];
  double *f2 = new double [3*Nmax];

  int mpisize, mpirank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  srand48(mpirank);
  double average = 0;
  for (int i=0; i<Ni; i++) {
    x[3*i+0] = drand48() * cycle - cycle / 2;
    x[3*i+1] = drand48() * cycle - cycle / 2;
    x[3*i+2] = drand48() * cycle - cycle / 2;
    p[i] = f[3*i+0] = f[3*i+1] = f[3*i+2] = 0;
  }
  for (int i=0; i<Ni; i++) {
    q[i] = drand48() - .5;
    average += q[i];
  }
  average /= Ni;
  for (int i=0; i<Ni; i++) {
    q[i] -= average;
  }

  FMM_Init(images);
  FMM_Partition(Ni, x, q, cycle);
  FMM(Ni, x, q, p, f, cycle);
  for (int i=0; i<Ni; i++) {
    x2[3*i+0] = x[3*i+0];
    x2[3*i+1] = x[3*i+1];
    x2[3*i+2] = x[3*i+2];
    q2[i] = q[i];
    p2[i] = f2[3*i+0] = f2[3*i+1] = f2[3*i+2] = 0;
  }
#if 0
  FMM_Ewald(Ni, x2, q2, p2, f2, ksize, alpha, cycle);
#else
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  double Xperiodic[3];
  int Nj = Ni, Nj3 = 3 * Ni;
  if (mpirank == 0) std::cout << "--- MPI direct sum ---------------" << std::endl;
  for (int irank=0; irank<mpisize; irank++) {
    if (mpirank == 0) std::cout << "Direct loop          : " << irank+1 << "/" << mpisize << std::endl;
    MPI_Shift(x2, Nj3, mpisize, mpirank);
    MPI_Shift(q2, Nj,  mpisize, mpirank);
    for (int i=0; i<Ni; i++) {
      double pp = 0, fx = 0, fy = 0, fz = 0;
      for (int ix=-prange; ix<=prange; ix++) {
	for (int iy=-prange; iy<=prange; iy++) {
	  for (int iz=-prange; iz<=prange; iz++) {
	    Xperiodic[0] = ix * cycle;
	    Xperiodic[1] = iy * cycle;
	    Xperiodic[2] = iz * cycle;
	    for (int j=0; j<Nj; j++) {
	      double dx = x[3*i+0] - x2[3*j+0] - Xperiodic[0];
	      double dy = x[3*i+1] - x2[3*j+1] - Xperiodic[1];
	      double dz = x[3*i+2] - x2[3*j+2] - Xperiodic[2];
	      double R2 = dx * dx + dy * dy + dz * dz;
	      double invR = 1 / std::sqrt(R2);
	      if (irank == mpisize-1 && i == j) invR = 0;
	      double invR3 = q2[j] * invR * invR * invR;
	      pp += q2[j] * invR;
	      fx += dx * invR3;
	      fy += dy * invR3;
	      fz += dz * invR3;
	    }
	  }
	}
      }
      p2[i] += pp;
      f2[3*i+0] -= fx;
      f2[3*i+1] -= fy;
      f2[3*i+2] -= fz;
    }
  }
  double localDipole[3] = {0, 0, 0};
  for (int i=0; i<Ni; i++) {
    for (int d=0; d<3; d++) localDipole[d] += x[3*i+d] * q[i];
  }
  int N;
  MPI_Allreduce(&Ni, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  double globalDipole[3];
  MPI_Allreduce(localDipole, globalDipole, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double norm = 0;
  for (int d=0; d<3; d++) {
    norm += globalDipole[d] * globalDipole[d];
  }
  double coef = 4 * M_PI / (3 * cycle * cycle * cycle);
  for (int i=0; i<Ni; i++) {
      p2[i] -= coef * norm / N / q[i];
      f2[3*i+0] -= coef * globalDipole[0];
      f2[3*i+1] -= coef * globalDipole[1];
      f2[3*i+2] -= coef * globalDipole[2];    
  }
#endif
  double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0, diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
  double v = 0, v2 = 0;
  for (int i=0; i<Ni; i++) {
    v += p[i] * q[i];
    v2 += p2[i] * q[i];
    diff2 += (f[3*i+0] - f2[3*i+0]) * (f[3*i+0] - f2[3*i+0])
           + (f[3*i+1] - f2[3*i+1]) * (f[3*i+1] - f2[3*i+1])
           + (f[3*i+2] - f2[3*i+2]) * (f[3*i+2] - f2[3*i+2]);
    norm2 += f2[3*i+0] * f2[3*i+0] + f2[3*i+1] * f2[3*i+1] + f2[3*i+2] * f2[3*i+2];
  }
  diff1 = (v - v2) * (v - v2);
  norm1 = v2 * v2;
  MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (mpirank == 0) {
    std::cout << "--- FMM vs. Ewald  ---------------" << std::endl;
    std::cout << std::setw(stringLength) << std::left
  	      << "Rel. L2 Error (pot)" << " : " << std::sqrt(diff3/norm3) << std::endl;
    if (std::abs(diff4) > 0) {
      std::cout << std::setw(stringLength) << std::left
	        << "Rel. L2 Error (acc)" << " : " << std::sqrt(diff4/norm4) << std::endl;
    }
  }

  delete[] x;
  delete[] q;
  delete[] p;
  delete[] f;
  delete[] x2;
  delete[] q2;
  delete[] p2;
  delete[] f2;
  MPI_Finalize();
}
