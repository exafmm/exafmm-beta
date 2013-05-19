#include "serialfmm.h"

/* to change real between float and double go to include/types.h */
typedef real            rvec[3];

extern "C" void FMMcalccoulomb_ij(
                   int ni, rvec* xi, real* qi, 
                   rvec* fi, real* energy, real convfctr,
                   real size, int periodicflag) {
  IMAGES = ((periodicflag & 0x1) == 0) ? 0 : 5;
  THETA = .5;
  Bodies bodies(ni);
  Cells cells,jcells;
  SerialFMM<Laplace> FMM;
  int i;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {

    i = B-bodies.begin();

    B->X[0] = xi[i][0];
    B->X[1] = xi[i][1];
    B->X[2] = xi[i][2];

    if( B->X[0] < -size/2 ) B->X[0] += size;
    if( B->X[1] < -size/2 ) B->X[1] += size;
    if( B->X[2] < -size/2 ) B->X[2] += size;
    if( B->X[0] >  size/2 ) B->X[0] -= size;
    if( B->X[1] >  size/2 ) B->X[1] -= size;
    if( B->X[2] >  size/2 ) B->X[2] -= size;

    B->SRC    =  qi[i];
    B->TRG[1] = 0.0; // -fi[i][0]; // gromacs input for force is unwanted MISC
    B->TRG[2] = 0.0; // -fi[i][1]; // gromacs input for force is unwanted MISC
    B->TRG[3] = 0.0; // -fi[i][2]; // gromacs input for force is unwanted MISC
    B->TRG[0] = 0.0; //  pi[i];    // gromacs only wants one output
    B->IBODY  = i;
    B->IPROC  = MPIRANK;

  }

  FMM.initialize();
  FMM.setDomain(bodies,0,size/2);
  FMM.bottomup(bodies,cells);
  jcells = cells;
  FMM.downward(cells,jcells);
  std::sort(bodies.begin(),bodies.end());
  FMM.writeTime();
  FMM.finalize();

  // dipole correction
  real dipole[3]={0.0,0.0,0.0};
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    for( int d=0; d!=3; ++d ) {
      dipole[d] += (B->X[d] - size/2) * B->SRC;
    }
  }
  real coef     = 4*M_PI / (3 * size  * size  * size );
  real coefTP   = coef * ( dipole[0] * dipole[0] + 
                              dipole[1] * dipole[1] + 
                              dipole[2] * dipole[2] );
  //real coefP    = coefTP / bodies.size();
  real coefF[3] = {coef*dipole[0], coef*dipole[1], coef*dipole[2]};

  // sending the data back
  real potenr=0.0;
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    int i = B-bodies.begin();
    fi[i][0] = -B->SRC * (B->TRG[1] - coefF[0]) * convfctr;
    fi[i][1] = -B->SRC * (B->TRG[2] - coefF[1]) * convfctr;
    fi[i][2] = -B->SRC * (B->TRG[3] - coefF[2]) * convfctr;
    // pi[i]    = 0.5 * (B->SRC * B->TRG[0] - coefP) * convfctr;
    potenr  += B->SRC * B->TRG[0];
  }
  *energy = (potenr - coefTP) * 0.5 * convfctr;
}
