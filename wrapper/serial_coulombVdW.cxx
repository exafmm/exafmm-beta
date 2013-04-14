#include "serialfmm.h"
extern "C" {
#include "md.h"
#include "vtgrapeproto.h"
}

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double, int tblno, double size, int periodicflag) {
  std::cout << "tblno: " << tblno << std::endl;
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 5;
  THETA = .5;
  Bodies bodies(ni),jbodies(nj);
  Cells cells,jcells;
  SerialFMM<Laplace> FMM;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    B->X[0] = xi[3*i+0];
    B->X[1] = xi[3*i+1];
    B->X[2] = xi[3*i+2];
    if( B->X[0] < -size/2 ) B->X[0] += size;
    if( B->X[1] < -size/2 ) B->X[1] += size;
    if( B->X[2] < -size/2 ) B->X[2] += size;
    if( B->X[0] >  size/2 ) B->X[0] -= size;
    if( B->X[1] >  size/2 ) B->X[1] -= size;
    if( B->X[2] >  size/2 ) B->X[2] -= size;
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
    if( B->X[0] < -size/2 ) B->X[0] += size;
    if( B->X[1] < -size/2 ) B->X[1] += size;
    if( B->X[2] < -size/2 ) B->X[2] += size;
    if( B->X[0] >  size/2 ) B->X[0] -= size;
    if( B->X[1] >  size/2 ) B->X[1] -= size;
    if( B->X[2] >  size/2 ) B->X[2] -= size;
    B->SRC  = qj[i];
  }

  FMM.initialize();
  FMM.setDomain(bodies,0,size/2);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);

  FMM.downward(cells,jcells);
  std::sort(bodies.begin(),bodies.end());
  FMM.writeTime();
  FMM.finalize();

#if 1
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
//    xi[3*i+0] = B->X[0];
//    xi[3*i+1] = B->X[1];
//    xi[3*i+2] = B->X[2];
//    qi[i]     = B->SRC;
    switch (tblno) {
    case 0 :
      fi[3*i+0] = -B->SRC * B->TRG[1];
      fi[3*i+1] = -B->SRC * B->TRG[2];
      fi[3*i+2] = -B->SRC * B->TRG[3];
      break;
    case 1 :
      fi[3*i+0] = 0.5 * B->SRC * B->TRG[0];
      break;
    }
  }
#else
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
      fi[3*i+0] += qi[i] * Fx;
      fi[3*i+1] += qi[i] * Fy;
      fi[3*i+2] += qi[i] * Fz;
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
      fi[3*i+0] += 0.5 * qi[i] * Po;
    }
    break;
  }
#endif
  // This is the correction factor from FMM to MD Ewald.
  double fc[3];
  for( int d=0; d!=3; ++d ) fc[d]=0;
  for( int i=0; i!=ni; ++i ) { 
    for( int d=0; d!=3; ++d ) { 
      fc[d] += qi[i] * xi[3*i+d];
    }   
  }
  if( tblno == 0 ) { 
    for( int i=0; i!=ni; ++i ) { 
      for( int d=0; d!=3; ++d ) {
        fi[3*i+d] -= 4.0 * M_PI * qi[i] * fc[d] / (3.0 * size * size * size);
      }
    }   
  } else {
    for( int i=0; i!=ni; ++i ) { 
      fi[3*i+0] += M_PI / (3.0 * size * size * size)
                * (fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2]) / ni;
    }   
  }
}

extern "C" void FMMcalcvdw_ij(int ni, double* xi, int* atypei, double* fi,
  int nj, double* xj, int* atypej, int nat, double* gscale, double* rscale,
  int tblno, double size, int periodicflag) {
  std::cout << "tblno: " << tblno << std::endl;
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
    if( B->X[0] < -size/2 ) B->X[0] += size;
    if( B->X[1] < -size/2 ) B->X[1] += size;
    if( B->X[2] < -size/2 ) B->X[2] += size;
    if( B->X[0] >  size/2 ) B->X[0] -= size;
    if( B->X[1] >  size/2 ) B->X[1] -= size;
    if( B->X[2] >  size/2 ) B->X[2] -= size;
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
    if( B->X[0] < -size/2 ) B->X[0] += size;
    if( B->X[1] < -size/2 ) B->X[1] += size;
    if( B->X[2] < -size/2 ) B->X[2] += size;
    if( B->X[0] >  size/2 ) B->X[0] -= size;
    if( B->X[1] >  size/2 ) B->X[1] -= size;
    if( B->X[2] >  size/2 ) B->X[2] -= size;
    B->SRC  = atypej[i] + .5;
  }

  FMM.initialize();
  FMM.setDomain(bodies,0,size/2);
  FMM.setVanDerWaals(nat,rscale,gscale);
  FMM.bottomup(bodies,cells);
  FMM.bottomup(jbodies,jcells);

  FMM.downward(cells,jcells);
  std::sort(bodies.begin(),bodies.end());
  FMM.writeTime();
  FMM.finalize();

#if 1
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
//    xi[3*i+0] = B->X[0];
//    xi[3*i+1] = B->X[1];
//    xi[3*i+2] = B->X[2];
//    atypei[i] = B->SRC;
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
#endif
}

extern "C" void fmmcalccoulomb_ij_(int *ni, double* xi, double* qi, double* fi,
  int *nj, double* xj, double* qj, double *rscale, int *tblno, double *size, int *periodicflag) {
  std::cout << "Starting FMM" << std::endl;
  FMMcalccoulomb_ij(*ni,xi,qi,fi,*nj,xj,qj,*rscale,*tblno-6,*size,*periodicflag);
}

extern "C" void fmmcalccoulomb_ij_exlist_(int *ni, double* xi, double* qi, double* fi,
  int *nj, double* xj, double* qj, double *rscale, int *tblno, double *size, int *periodicflag,
  int *numex, int* natex) {
  std::cout << "Starting FMM" << std::endl;
  FMMcalccoulomb_ij(*ni,xi,qi,fi,*nj,xj,qj,*rscale,*tblno-6,*size,*periodicflag);
  switch (*tblno-6) {
  case 0 :
#if 0
    for( int i=0,ic=0; i<*ni; i++ ) {
      for( int j=0; j<numex[i]; j++,ic++ ) natex[ic]--;
    }
//    MR3calccoulomb_nlist_ij_host(*ni,xi,qi,fi,*nj,xj,qj,
//                                *rscale,*tblno-6,*size,*periodicflag&3,numex,natex,-1.0);
//    MR3calccoulomb_nlist_ij_emu(*ni,xi,qi,fi,*nj,xj,qj,
//			        *rscale,*tblno-6,*size,*periodicflag&3,numex,natex,-1.0);
    for( int i=0,ic=0; i<*ni; i++ ) {
      for( int j=0; j<numex[i]; j++,ic++ ) natex[ic]++;
    }
#else
    for( int i=0,ic=0; i!=*ni; ++i ) {
      double Fx = 0, Fy = 0, Fz = 0;
      for( int in=0; in!=numex[i]; in++,ic++ ) {
        int j = natex[ic]-1;
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
        fi[3*j+0] += qi[i] * Fx;
        fi[3*j+1] += qi[i] * Fy;
        fi[3*j+2] += qi[i] * Fz;
      }
      fi[3*i+0] -= qi[i] * Fx;
      fi[3*i+1] -= qi[i] * Fy;
      fi[3*i+2] -= qi[i] * Fz;
    }
#endif
    break;
  case 1:
    double Total = 0;
    int ic = 0;
    for( int i=0; i!=*ni; ++i ) {
      double Po = 0;
      for( int in=0; in!=numex[i]; in++,ic++ ) {
        int j = natex[ic]-1;
        double dx = xi[3*i+0] - xj[3*j+0];
        double dy = xi[3*i+1] - xj[3*j+1];
        double dz = xi[3*i+2] - xj[3*j+2];
        double R2 = dx * dx + dy * dy + dz * dz;
        double invR = 1 / std::sqrt(R2);
        assert( i != j );
        if( R2 == 0 ) invR = 0;
        Po += qj[j] * invR;
      }
      fi[3*i+0] -= qi[i] * Po;
      Total -= qi[i] * Po;
    }
    std::cout << "Total exclusion: " << Total << std::endl;
    std::cout << "Total # of exclusion: " << ic << std::endl;
    break;
  }
}


extern "C" void fmmcalcvdw_ij_(int *ni, double* xi, int* atypei, double* fi,
  int *nj, double* xj, int* atypej, int *nat, double* gscale, double* rscale,
  int *tblno, double *size, int *periodicflag) {
  std::cout << "Starting FMM" << std::endl;
  for( int i=0; i<*ni; i++ ) atypei[i]--;
  for( int i=0; i<*nj; i++ ) atypej[i]--;
  FMMcalcvdw_ij(*ni,xi,atypei,fi,*nj,xj,atypej,*nat,gscale,rscale,*tblno,*size,*periodicflag);
  for( int i=0; i<*ni; i++ ) atypei[i]++;
  for( int i=0; i<*nj; i++ ) atypej[i]++;
}

extern "C" void fmmcalcvdw_ij_exlist_(int *ni, double* xi, int* atypei, double* fi,
  int *nj, double* xj, int* atypej, int *nat, double* gscale, double* rscale,
  int *tblno, double *size, int *periodicflag, int* numex, int* natex) {
  std::cout << "Starting FMM" << std::endl;
  for( int i=0; i<*ni; i++ ) atypei[i]--;
  for( int i=0; i<*nj; i++ ) atypej[i]--;
  for( int i=0,ic=0; i<*ni; i++ ) {
    for( int j=0; j<numex[i]; j++,ic++ ) natex[ic]--;
  }
#if 0
  MR3calcvdw_ij(*ni,xi,atypei,fi,*nj,xj,atypej,*nat,gscale,rscale,*tblno,*size,*periodicflag);
  MR3calcvdw_nlist_ij_emu(*ni,xi,atypei,fi,*nj,xj,atypej,*nat,gscale,rscale,*tblno,
                          *size,(*periodicflag & 3),numex,natex,-1.0);
#else
  FMMcalcvdw_ij(*ni,xi,atypei,fi,*nj,xj,atypej,*nat,gscale,rscale,*tblno,*size,*periodicflag);
//  MR3calcvdw_nlist_ij_host(*ni,xi,atypei,fi,*nj,xj,atypej,*nat,gscale,rscale,*tblno,
//                          *size,(*periodicflag & 3),numex,natex,-1.0);
#endif
  for( int i=0; i<*ni; i++ ) atypei[i]++;
  for( int i=0; i<*nj; i++ ) atypej[i]++;
  for( int i=0,ic=0; i<*ni; i++ ) {
    for( int j=0; j<numex[i]; j++,ic++ ) natex[ic]++;
  }
}

