#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/time.h>

extern "C" void FMM_Init(double eps2, double kreal, double kimag, int ncrit, int threads, const char * path,
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
extern "C" void FMM_Verify_Accuracy(int & t, double potRel, double accRel);
extern "C" bool FMM_Only_Accuracy();
extern "C" void FMM_Verify_Time(int & t, double totalFMM);
extern "C" void FMM_Verify_End();

double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

void Verify(int n, double * vb, double * vd, int & t, int verbose) {
  double potDif = 0, potNrm = 0;
  for (int i=0; i<n; i++) {
    potDif += (vb[i] - vd[i]) * (vb[i] - vd[i]);
    potNrm += vd[i] * vd[i];
  }
  double potDifGlob = 0, potNrmGlob = 0;
  MPI_Reduce(&potDif, &potDifGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&potNrm, &potNrmGlob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  double potRel = std::sqrt(potDifGlob/potNrmGlob);
  if (verbose) {
    std::cout << "--- FMM vs. direct ---------------" << std::endl;
    std::cout << std::setw(20) << std::left << std::scientific
              << "Rel. L2 Error" << " : " << potRel << std::endl;
  }
  FMM_Verify_Accuracy(t, potRel, potRel);
}

int main(int argc, char ** argv) {
  const int Nmax = 10000000;
  const int ncrit = 1000;
  const int threads = 16;
  const int verbose = 0;
  const double eps2 = 0.0;
  const double kreal = 1.0;
  const double kimag = 0.1;
  const char * path = "./";
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
  int nb = 5000 / mpisize;
  int nv = 10000 / mpisize;
  int nit = 10;

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

  FMM_Init(eps2, kreal, kimag, ncrit, threads, path, nb, xb, yb, zb, vb, nv, xv, yv, zv, vv);
  FMM_Partition(nb, xb, yb, zb, vb, nv, xv, yv, zv, vv);
  FMM_BuildTree();

  for (int t=0; t<10; t++) {
    for (int i=0; i<nb; i++) {
      vb[i] = 1.0 / nb;
      vi[i] = 0;
      vd[i] = 0;
    }
    FMM_B2B(vi, vb, verbose);
    Direct(100, xb, yb, zb, vd, nb, xb, yb, zb, vb);
    Verify(100, vi, vd, t, mpirank == 0);
    if (t == -1) break;
  }
  FMM_Verify_End();

  if (!FMM_Only_Accuracy()) {
    for (int t=0; t<10; t++) {
      double tic = get_time();
      for (int it=0; it<nit; it++) {
        for (int i=0; i<nb; i++) {
          vb[i] = 1.0 / nb;
          vi[i] = 0;
          vd[i] = 0;
        }
        FMM_B2B(vi, vb, verbose);
      }
      double toc = get_time();
      FMM_Verify_Time(t, (toc-tic)/nit);
      if (t == -1) break;
    }
    FMM_Verify_End();
  }

  for (int t=0; t<10; t++) {
    for (int i=0; i<nb; i++) {
      vb[i] = 0;
      vd[i] = 0;
    }
    for (int i=0; i<nv; i++) {
      vv[i] = 1.0 / nv;
    }
    FMM_V2B(vb, vv, verbose);
    Direct(100, xb, yb, zb, vd, nv, xv, yv, zv, vv);
    Verify(100, vb, vd, t, mpirank == 0);
    if (t == -1) break;
  }
  FMM_Verify_End();

  if (!FMM_Only_Accuracy()) {
    for (int t=0; t<10; t++) {
      double tic = get_time();
      for (int it=0; it<nit; it++) {
        for (int i=0; i<nb; i++) {
          vb[i] = 0;
          vd[i] = 0;
        }
        for (int i=0; i<nv; i++) {
          vv[i] = 1.0 / nv;
        }
        FMM_V2B(vb, vv, verbose);
      }
      double toc = get_time();
      FMM_Verify_Time(t, (toc-tic)/nit);
      if (t == -1) break;
    }
    FMM_Verify_End();
  }

  for (int t=0; t<10; t++) {
    for (int i=0; i<nb; i++) {
      vb[i] = 1.0 / nb;
    }
    for (int i=0; i<nv; i++) {
      vv[i] = 0;
      vd[i] = 0;
    }
    FMM_B2V(vv, vb, verbose);
    Direct(100, xv, yv, zv, vd, nb, xb, yb, zb, vb);
    Verify(100, vv, vd, t, mpirank == 0);
    if (t == -1) break;
  }
  FMM_Verify_End();

  if (!FMM_Only_Accuracy()) {
    for (int t=0; t<10; t++) {
      double tic = get_time();
      for (int it=0; it<nit; it++) {
        for (int i=0; i<nb; i++) {
          vb[i] = 1.0 / nb;
        }
        for (int i=0; i<nv; i++) {
          vv[i] = 0;
          vd[i] = 0;
        }
        FMM_B2V(vv, vb, verbose);
      }
      double toc = get_time();
      FMM_Verify_Time(t, (toc-tic)/nit);
      if (t == -1) break;
    }
    FMM_Verify_End();
  }

  for (int t=0; t<10; t++) {
    for (int i=0; i<nv; i++) {
      vv[i] = 1.0 / nv;
      vi[i] = 0;
      vd[i] = 0;
    }
    FMM_V2V(vi, vv, verbose);
    Direct(100, xv, yv, zv, vd, nv, xv, yv, zv, vv);
    Verify(100, vi, vd, t, mpirank == 0);
    if (t == -1) break;
  }
  FMM_Verify_End();

  if (!FMM_Only_Accuracy()) {
    for (int t=0; t<10; t++) {
      double tic = get_time();
      for (int it=0; it<nit; it++) {
        for (int i=0; i<nv; i++) {
          vv[i] = 1.0 / nv;
          vi[i] = 0;
          vd[i] = 0;
        }
        FMM_V2V(vi, vv, verbose);
      }
      double toc = get_time();
      FMM_Verify_Time(t, (toc-tic)/nit);
      if (t == -1) break;
    }
    FMM_Verify_End();
  }

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
