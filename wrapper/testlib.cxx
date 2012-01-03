/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "mr3.h"

const double R2MIN = 0.25;
const double R2MAX = 64;

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double rscale, int tblno, double size, int periodicflag);

extern "C" void FMMcalcvdw_ij(int ni, double* xi, int* atypei, double* fi,
  int nj, double* xj, int* atypej, int nat, double* gscale, double* rscale,
  int tblno, double size, int periodicflag);

extern "C" void MPI_ShiftI(int *var, int n, int mpisize, int mpirank);
extern "C" void MPI_Shift(double *var, int n, int mpisize, int mpirank);

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);
  int mpisize, mpirank;
  MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpirank);
  const int N = 10000;
  const int nat = 16;
  const double size = 2;
  double *xi     = new double [3*N];
  double *qi     = new double [N];
  double *pi     = new double [3*N];
  double *fi     = new double [3*N];
  double *pd     = new double [3*N];
  double *fd     = new double [3*N];
  double *xj     = new double [3*N];
  double *qj     = new double [N];
  int *atypei    = new int [N];
  int *atypej    = new int [N];
  double *rscale = new double [64*64];
  double *gscale = new double [64*64];

  srand48(mpirank);
  for( int i=0; i!=N; ++i ) {
    xi[3*i+0] = drand48() * size;
    xi[3*i+1] = drand48() * size;
    xi[3*i+2] = drand48() * size;
    qi[i] = drand48()*2.0-1.0;
    pi[3*i+0] = pi[3*i+1] = pi[3*i+2] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[3*i+0] = pd[3*i+1] = pd[3*i+2] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
    atypei[i] = drand48() * nat;
  }
  for( int i=0; i!=N; ++i ) {
    xj[3*i+0] = drand48() * size;
    xj[3*i+1] = drand48() * size;
    xj[3*i+2] = drand48() * size;
    qj[i] = drand48()*2.0-1.0;
    atypej[i] = drand48() * nat;
  }
  srand48(0);
  for( int i=0; i!=nat; ++i ) {
    gscale[i*nat+i] = drand48();
    drand48();
    rscale[i*nat+i] = drand48();
  }
  for( int i=0; i!=nat; ++i ) {
    for( int j=0; j!=nat; ++j ) {
      if( i != j ) {
        gscale[i*nat+j] = sqrt(gscale[i*nat+i]*gscale[j*nat+j]);
        rscale[i*nat+j] = (sqrt(rscale[i*nat+i]) + sqrt(rscale[j*nat+j])) * 0.5;
        rscale[i*nat+j] *= rscale[i*nat+j];
      }
    }
  }

  FMMcalccoulomb_ij(N, xi, qi, pi, N, xj, qj, 0.0, 1, size, 0);
  FMMcalccoulomb_ij(N, xi, qi, fi, N, xj, qj, 0.0, 0, size, 0);
  for( int irank=0; irank!=mpisize; ++irank ) {
    MPI_Shift(xj,3*N,mpisize,mpirank);
    MPI_Shift(qj,N,mpisize,mpirank);
#if 1
    MR3calccoulomb_ij(N, xi, qi, pd, N, xj, qj, 1.0, 1, size, 0);
    MR3calccoulomb_ij(N, xi, qi, fd, N, xj, qj, 1.0, 0, size, 0);
#else
    for( int i=0; i!=N; ++i ) {
      double P = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=0; j!=N; ++j ) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        double invR = 1 / std::sqrt(R2);
        if( R2 == 0 ) invR = 0;
        double invR3 = qj[j] * invR * invR * invR;
        P += qj[j] * invR;
        Fx += dx * invR3;
        Fy += dy * invR3;
        Fz += dz * invR3;
      }
      pd[3*i+0] += P;
      fd[3*i+0] += Fx;
      fd[3*i+1] += Fy;
      fd[3*i+2] += Fz;
    }
#endif
  }
  double Pd = 0, Pn = 0, Fd = 0, Fn = 0;
  for( int i=0; i!=N; ++i ) {
    Pd += (pi[3*i+0] - pd[3*i+0]) * (pi[3*i+0] - pd[3*i+0]);
    Pn += pd[3*i+0] * pd[3*i+0];
    Fd += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
        + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
        + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    Fn += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
    pi[3*i+0] = pi[3*i+1] = pi[3*i+2] = 0;
    fi[3*i+0] = fi[3*i+1] = fi[3*i+2] = 0;
    pd[3*i+0] = pd[3*i+1] = pd[3*i+2] = 0;
    fd[3*i+0] = fd[3*i+1] = fd[3*i+2] = 0;
  }
  std::cout << "Coulomb       potential @ rank " << mpirank << " : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Coulomb       force     @ rank " << mpirank << " : " << sqrtf(Fd/Fn) << std::endl;

  FMMcalcvdw_ij(N,xi,atypei,pi,N,xj,atypej,nat,gscale,rscale,3,size,0);
  FMMcalcvdw_ij(N,xi,atypei,fi,N,xj,atypej,nat,gscale,rscale,2,size,0);
  for( int irank=0; irank!=mpisize; ++irank ) {
    MPI_Shift(xj,3*N,mpisize,mpirank);
    MPI_ShiftI(atypej,N,mpisize,mpirank);
#if 0
    MR3calcvdw_ij(N,xi,atypei,pd,N,xj,atypej,nat,gscale,rscale,3,size,0);
    MR3calcvdw_ij(N,xi,atypei,fd,N,xj,atypej,nat,gscale,rscale,2,size,0);
#else
    for( int i=0; i!=N; ++i ) {
      double P = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=0; j!=N; ++j ) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        if( R2 != 0 ) {
          double rs = rscale[atypei[i]*nat+atypej[j]];
          double gs = gscale[atypei[i]*nat+atypej[j]];
          double R2s = R2 * rs;
          if( R2MIN <= R2s && R2s < R2MAX ) {
            double invR2 = 1.0 / R2s;
            double invR6 = invR2 * invR2 * invR2;
            double dtmp = gs * invR6 * invR2 * (2.0 * invR6 - 1.0);
            P += gs * invR6 * (invR6 - 1.0);
            Fx += dx * dtmp;
            Fy += dy * dtmp;
            Fz += dz * dtmp;
          }
        }
      }
      pd[3*i+0] += P;
      fd[3*i+0] += Fx;
      fd[3*i+1] += Fy;
      fd[3*i+2] += Fz;
    }
#endif
  }
  Pd = Pn = Fd = Fn = 0;
  for( int i=0; i!=N; ++i ) {
    Pd += (pi[3*i+0] - pd[3*i+0]) * (pi[3*i+0] - pd[3*i+0]);
    Pn += pd[3*i+0] * pd[3*i+0];
    Fd += (fi[3*i+0] - fd[3*i+0]) * (fi[3*i+0] - fd[3*i+0])
        + (fi[3*i+1] - fd[3*i+1]) * (fi[3*i+1] - fd[3*i+1])
        + (fi[3*i+2] - fd[3*i+2]) * (fi[3*i+2] - fd[3*i+2]);
    Fn += fd[3*i+0] * fd[3*i+0] + fd[3*i+1] * fd[3*i+1] + fd[3*i+2] * fd[3*i+2];
  }
  std::cout << "Van der Waals potential @ rank " << mpirank << " : " << sqrtf(Pd/Pn) << std::endl;
  std::cout << "Van der Waals force     @ rank " << mpirank << " : " << sqrtf(Fd/Fn) << std::endl;

  delete[] xi;
  delete[] qi;
  delete[] pi;
  delete[] fi;
  delete[] pd;
  delete[] fd;
  delete[] xj;
  delete[] qj;
  delete[] atypei;
  delete[] atypej;
  delete[] rscale;
  delete[] gscale;

  MPI_Finalize();
}
