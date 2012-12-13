#include "serialfmm.h"

extern "C" void FMMcalccoulomb_ij(int ni, double* xi, double* qi, double* fi,
  int nj, double* xj, double* qj, double, int tblno, double size, int periodicflag) {
  std::cout << "tblno: " << tblno << std::endl;
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 3;
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
//    printf("%f %f %f\n", xi[0], xi[1], xi[2]);
//    printf("%f %f %f\n", xi[96], xi[97], xi[98]);
//    printf("%f %f %f\n", xi[240], xi[241], xi[242]);

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
  FMM.finalize();

#if 1
   std::cout<< "check 1" <<std::endl;
  //for( B_iter B=bodies.begin(); B!=bodies.end()-ni+300; ++B ) {
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();

 //   std::cout << xi[3*i+0] << " " << B->X[0] << std::endl;
    //assert(fabs(xi[3*i+0] - B->X[0])< 1e-6);
    //xi[3*i+0] = B->X[0];
   // xi[3*i+1] = B->X[1];
    //xi[3*i+2] = B->X[2];

    qi[i]     = B->SRC;
    switch (tblno) {
    case 0 :
      fi[3*i+0] = -B->SRC * B->TRG[1];
      fi[3*i+1] = -B->SRC * B->TRG[2];
      fi[3*i+2] = -B->SRC * B->TRG[3];
      break;
    case 1 :
      fi[3*i+0] = B->SRC * B->TRG[0];
      break;
    }
  }

//    printf("%f %f %f\n", xi[0], xi[1], xi[2]);
//    printf("%f %f %f\n", xi[96], xi[97], xi[98]);
//   printf("%f %f %f\n", xi[240], xi[241], xi[242]);

  double fc[3];
  for( int d=0; d!=3; ++d ) fc[d]=0;
  for( int i=0; i!=ni; ++i ) { 
    for( int d=0; d!=3; ++d ) { 
      fc[d] += qi[i] * xi[3*i+d];
    }   
  }
  if( (tblno % 2) == 1 ) { 
    for( int i=0; i!=ni; ++i ) { 
      fi[3*i+0] += 2.0 * M_PI / (3.0 * size * size * size)
                * (fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2]) / ni;
    }   
  } else {
    for( int i=0; i!=ni; ++i ) { 
      for( int d=0; d!=3; ++d ) {
        fi[3*i+d] -= 4.0 * M_PI * qi[i] * fc[d] / (3.0 * size * size * size);
      }
    }   
  }
#else
   std::cout<< "check not 1" <<std::endl;
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
   //   std::cout << i << " " << fi[3*i+0] << std::endl;
    }
    break;
  } 
#endif 
}

extern "C" void getaccfmm_(int *ni, double* xi, double* qi, double* fi,
  int *nj, double* xj, double* qj, double *rscale, int *tblno, double *size, int *periodicflag) {
  std::cout << "Starting FMM" << std::endl;
  // std::cout<< "check" <<std::endl;
  FMMcalccoulomb_ij(*ni,xi,qi,fi,*nj,xj,qj,*rscale,*tblno,*size,*periodicflag);
}
