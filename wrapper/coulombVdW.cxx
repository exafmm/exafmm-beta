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

extern "C" void MPI_ShiftI(int *var, int n, int mpisize, int mpirank) {
  int *buf = new int [n];
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
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

extern "C" void MPI_Shift(double *var, int n, int mpisize, int mpirank) {
  double *buf = new double [n];
  const int isend = (mpirank + 1          ) % mpisize;
  const int irecv = (mpirank - 1 + mpisize) % mpisize;
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

extern "C" void FMMcalccoulomb_fp_ij(int ni, double* xi, double* qi, double* fi, double* p,
  int nj, double* xj, double* qj, double size, int periodicflag) {
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
    B->TRG[1] = -fi[3*i+0];
    B->TRG[2] = -fi[3*i+1];
    B->TRG[3] = -fi[3*i+2];
    B->TRG[0] = p[i];
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
  FMM.finalize();

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    qi[i]     = B->SRC;
    fi[3*i+0] = -B->TRG[1];
    fi[3*i+1] = -B->TRG[2];
    fi[3*i+2] = -B->TRG[3];
    p[i] = B->TRG[0];
  }
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
  FMM.finalize();

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
    MPI_Shift(xj,3*nj,MPISIZE,MPIRANK);
    MPI_Shift(qj,nj,MPISIZE,MPIRANK);
    switch (tblno) {
    case 0 :
      for( int i=0; i!=ni; ++i ) {
        double Fx = 0, Fy = 0, Fz = 0;
        for( int j=0; j!=nj; ++j ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          double invR = 1 / std::sqrt(R2);
          if( R2 == 0 ) invR = 0;
          double invR3 = qj[j] * invR * invR * invR;
          Fx += dx * invR3;
          Fy += dy * invR3;
          Fz += dz * invR3;
        }
        fi[3*i+0] += Fx;
        fi[3*i+1] += Fy;
        fi[3*i+2] += Fz;
      }
      break;
    case 1:
      for( int i=0; i!=ni; ++i ) {
        double Po = 0;
        for( int j=0; j!=nj; ++j ) {
          double dx = xi[3*i+0] - xj[3*j+0];
          double dy = xi[3*i+1] - xj[3*j+1];
          double dz = xi[3*i+2] - xj[3*j+2];
          double R2 = dx * dx + dy * dy + dz * dz;
          double invR = 1 / std::sqrt(R2);
          if( R2 == 0 ) invR = 0;
          Po += qj[j] * invR;
        }
        fi[3*i+0] += Po;
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
  ParallelFMM<VanDerWaals> FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    B->SRC  = atypei[i] + .5;
    switch (tblno) {
    case 2 :
      B->TRG[1] = -fi[3*i+0];
      B->TRG[2] = -fi[3*i+1];
      B->TRG[3] = -fi[3*i+2];
      break;
    case 3 :
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
    B->SRC  = atypej[i] + .5;
  }

  FMM.initialize();
  FMM.setGlobDomain(bodies,size/2,size/2);
  FMM.setVanDerWaals(nat,rscale,gscale);
  FMM.octsection(bodies);
  FMM.octsection(jbodies);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);
  FMM.commBodies(jcells);
  FMM.commCells(jbodies,jcells);
  FMM.downward(cells,jcells);
  FMM.unpartition(bodies);
  std::sort(bodies.begin(),bodies.end());
  FMM.finalize();

#if 1
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    atypei[i] = B->SRC;
    switch (tblno) {
    case 2 :
      fi[3*i+0] = -B->TRG[1];
      fi[3*i+1] = -B->TRG[2];
      fi[3*i+2] = -B->TRG[3];
      break;
    case 3 :
      fi[3*i+0] = B->TRG[0];
      break;
    }
  }
#else
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Shift(xj,3*nj,MPISIZE,MPIRANK);
    MPI_ShiftI(atypej,nj,MPISIZE,MPIRANK);
    switch (tblno) {
    case 2 :
      for( int i=0; i!=ni; ++i ) {
        double Fx = 0, Fy = 0, Fz = 0;
        for( int j=0; j!=nj; ++j ) {
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
              Fx += dx * dtmp;
              Fy += dy * dtmp;
              Fz += dz * dtmp;
            }
          }
        }
        fi[3*i+0] += Fx;
        fi[3*i+1] += Fy;
        fi[3*i+2] += Fz;
      }
      break;
    case 3:
      for( int i=0; i!=ni; ++i ) {
        double Po = 0;
        for( int j=0; j!=nj; ++j ) {
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
              Po += gs * invR6 * (invR6 - 1.0);
            }
          }
        }
        fi[3*i+0] += Po;
      }
      break;
    }
  }
#endif
}
