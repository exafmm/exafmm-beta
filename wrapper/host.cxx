#include "parallelfmm.h"

#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 100.0f

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

extern "C" void FMMcalccoulomb_ij_host(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double rscale, int tblno, double size, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  THETA = 1/sqrtf(3);
  vect shift = size/2;
  Bodies bodies(ni),jbodies(nj);
  Cells cells,jcells;
  ParallelFMM FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0]   = xi[3*i+0];
    B->X[1]   = xi[3*i+1];
    B->X[2]   = xi[3*i+2];
    B->SRC[0] = qi[i];
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
    B->X[0]   = xj[3*i+0];
    B->X[1]   = xj[3*i+1];
    B->X[2]   = xj[3*i+2];
    B->SRC[0] = qj[i];
  }

  FMM.initialize();
  FMM.setGlobDomain(bodies,shift,size/2);
  FMM.octsection(bodies);
  FMM.octsection(jbodies);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);
  FMM.commBodies(jcells);
  FMM.commCells(jbodies,jcells);

  FMM.downward(cells,jcells,1);
  FMM.unpartition(bodies);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    qi[i]     = B->SRC[0];
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

/*
  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Shift(xj,3*nj);
    MPI_Shift(qj,nj);
    for( int i=0; i<ni; i++ ) {
      for( int j=0; j<nj; j++ ) {
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double r2 = dx * dx + dy * dy + dz * dz;
        if( r2 != 0 ) {
          fi[3*i+0] += qj[j] / std::sqrt(r2);
        }
      }
    }
  }
*/
}

extern "C" void FMMcalcvdw_ij_host(int ni, double* xi, int* atypei, double* fi,
  int nj, double* xj, int* atypej, int nat, double* gscale, double* rscale,
  int tblno, double size, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
  THETA = 1/sqrtf(3);
  vect shift = size/2;
  Bodies bodies(ni),jbodies(nj);
  Cells cells,jcells;
  SerialFMM FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0]   = xi[3*i+0];
    B->X[1]   = xi[3*i+1];
    B->X[2]   = xi[3*i+2];
    B->SRC[1] = atypei[i] + .5;
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
    B->X[0]   = xj[3*i+0];
    B->X[1]   = xj[3*i+1];
    B->X[2]   = xj[3*i+2];
    B->SRC[1] = atypej[i] + .5;
  }


  FMM.setKernel("CoulombVdW");
  FMM.setDomain(bodies,shift,size/2);
  FMM.setVanDerWaals(nat,rscale,gscale);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);
  FMM.downward(cells,jcells,1);
  std::sort(bodies.begin(),bodies.end());

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    xi[3*i+0] = B->X[0];
    xi[3*i+1] = B->X[1];
    xi[3*i+2] = B->X[2];
    atypei[i] = B->SRC[1];
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

/*
  for( int i=0; i<ni; i++ ) {
    for( int j=0; j<nj; j++ ) {
      double dx = xi[3*i+0] - xj[3*j+0];
      double dy = xi[3*i+1] - xj[3*j+1];
      double dz = xi[3*i+2] - xj[3*j+2];
      double r2 = dx * dx + dy * dy + dz * dz;
      if( r2 != 0 ) {
        double rs = rscale[atypei[i]*nat+atypej[j]];
        double gs = gscale[atypei[i]*nat+atypej[j]];
        double rrs = r2 * rs;
        double r1 = 1.0 / rrs;
        double r6 = r1 * r1 * r1;
        double dtmp = gs * r6 * r1 * (2.0 * r6 - 1.0);
        fi[3*i+0] += dtmp * dx;
        fi[3*i+1] += dtmp * dy;
        fi[3*i+2] += dtmp * dz;
      }
    }
  }
*/
}
