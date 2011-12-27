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
#include "parallelfmm.h"

#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 100.0f

void MPI_Shift(int *var, int n) {
  int *buf = new int [n];
  const int isend = (MPIRANK + 1          ) % MPISIZE;
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;
  MPI_Request sreq, rreq;

  MPI_Isend(var,n,MPI_INT,irecv,1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(buf,n,MPI_INT,isend,1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  int i;
  for( i=0; i!=n; ++i ) {
    var[i] = buf[i];
  }
  delete[] buf;
}

void MPI_Shift(double *var, int n) {
  double *buf = new double [n];
  const int isend = (MPIRANK + 1          ) % MPISIZE;
  const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;
  MPI_Request sreq, rreq;

  MPI_Isend(var,n,MPI_DOUBLE,irecv,1,MPI_COMM_WORLD,&sreq);
  MPI_Irecv(buf,n,MPI_DOUBLE,isend,1,MPI_COMM_WORLD,&rreq);
  MPI_Wait(&sreq,MPI_STATUS_IGNORE);
  MPI_Wait(&rreq,MPI_STATUS_IGNORE);
  int i;
  for( i=0; i!=n; ++i ) {
    var[i] = buf[i];
  }
  delete[] buf;
}

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double, int tblno, double size, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  THETA = .5;
  Bodies bodies(ni),jbodies(nj);
  Cells cells,jcells;
  ParallelFMM<Laplace> FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    B->SRC  = qi[i];
    switch (tblno) {
    case 0 :
      B->TRG[1] = -fi[3*i+0];
      B->TRG[2] = -fi[3*i+1];
      B->TRG[3] = -fi[3*i+2];
      break;
    case 1 :
      B->TRG[0] = fi[3*i+0];
      break;
    }
    B->IBODY = i;
    B->IPROC = MPIRANK;
  }

  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
    int i = B-jbodies.begin();
    B->X[0] = xj[3*i+0];
    B->X[1] = xj[3*i+1];
    B->X[2] = xj[3*i+2];
    B->SRC  = qj[i];
  }

  FMM.initialize();
  FMM.setGlobDomain(bodies,size/2,size/2);
  FMM.octsection(bodies);
  FMM.octsection(jbodies);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);
  FMM.commBodies(jcells);
  FMM.commCells(jbodies,jcells);

  FMM.downward(cells,jcells);
  FMM.unpartition(bodies);
  std::sort(bodies.begin(),bodies.end());

#if 1
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    qi[i]     = B->SRC;
    switch (tblno) {
    case 0 :
      fi[3*i+0] = -B->TRG[1];
      fi[3*i+1] = -B->TRG[2];
      fi[3*i+2] = -B->TRG[3];
      break;
    case 1 :
      fi[3*i+0] = B->TRG[0];
      break;
    }
  }
#else
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Shift(xj,3*nj);
    MPI_Shift(qj,nj);
    switch (tblno) {
    case 0 :
      for( int i=0; i<ni; i++ ) {
        for( int j=0; j<nj; j++ ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          double invR = 1 / sqrtf(R2);
          if( R2 == 0 ) invR = 0;
          double invR3 = qj[j] * invR * invR * invR;
          fi[3*i+0] += dx * invR3;
          fi[3*i+1] += dy * invR3;
          fi[3*i+2] += dz * invR3;
        }
      }
      break;
    case 1:
      for( int i=0; i<ni; i++ ) {
        for( int j=0; j<nj; j++ ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          double invR = 1 / sqrtf(R2);
          if( R2 == 0 ) invR = 0;
          fi[3*i+0] += qj[j] * invR;
        }
      }
      break;
    }
  }
#endif
}

extern "C" void FMMcalcvdw_ij(int ni, double* xi, int* atypei, double* fi,
  int nj, double* xj, int* atypej, int nat, double* gscale, double* rscale,
  int tblno, double size, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  THETA = .5;
  Bodies bodies(ni),jbodies(nj);
  Cells cells,jcells;
  SerialFMM<VanDerWaals> FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    B->SRC  = atypei[i] + .5;
    switch (tblno) {
    case 2 :
      B->TRG[1] = fi[3*i+0];
      B->TRG[2] = fi[3*i+1];
      B->TRG[3] = fi[3*i+2];
      break;
    case 3 :
      B->TRG[0] = fi[3*i+0];
      break;
    }
    B->IBODY = i;
  }

  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
    int i = B-jbodies.begin();
    B->X[0] = xj[3*i+0];
    B->X[1] = xj[3*i+1];
    B->X[2] = xj[3*i+2];
    B->SRC  = atypej[i] + .5;
  }


  FMM.setDomain(bodies,size/2,size/2);
  FMM.setVanDerWaals(nat,rscale,gscale);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);
  FMM.downward(cells,jcells);
  std::sort(bodies.begin(),bodies.end());

#if 1
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    atypei[i] = B->SRC;
    switch (tblno) {
    case 2 :
      fi[3*i+0] = B->TRG[1];
      fi[3*i+1] = B->TRG[2];
      fi[3*i+2] = B->TRG[3];
      break;
    case 3 :
      fi[3*i+0] = B->TRG[0];
      break;
    }
  }
#else
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Shift(xj,3*nj);
    MPI_Shift(atypej,nj);
    switch (tblno) {
    case 2 :
      for( int i=0; i<ni; i++ ) {
        for( int j=0; j<nj; j++ ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          if( R2 != 0 ) {
            double rs = rscale[atypei[i]*nat+atypej[j]];
            double gs = gscale[atypei[i]*nat+atypej[j]];
            double rrs = R2 * rs;
            double invR = 1.0 / std::sqrt(rrs);
            double invR2 = invR * invR;
            double invR6 = invR2 * invR2 * invR2;
            double dtmp = gs * invR6 * invR * (2.0 * invR6 - 1.0);
            fi[3*i+0] += dtmp * dx;
            fi[3*i+1] += dtmp * dy;
            fi[3*i+2] += dtmp * dz;
          }
        }
      }
      break;
    case 3:
      for( int i=0; i<ni; i++ ) {
        for( int j=0; j<nj; j++ ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          if( R2 != 0 ) {
            double rs = rscale[atypei[i]*nat+atypej[j]];
            double gs = gscale[atypei[i]*nat+atypej[j]];
            double rrs = R2 * rs;
            double invR2 = 1.0 / rrs;
            double invR6 = invR2 * invR2 * invR2;
            fi[3*i+0] += gs * invR6 * (invR6 - 1.0);
          }
        }
      }
      break;
    }
  }
#endif
}
