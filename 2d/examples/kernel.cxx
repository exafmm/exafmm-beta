#include <fstream>
#include "kernel.h"
#include <iostream>
#include <vector>

int main() {
  Bodies bodies(1), bodies2(1), jbodies(1);
  Cells cells(4);
  Kernel kernel;
  const real_t theta = 0.5;
  const real_t R = 2 / theta;

  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    B->X[0] = 2 * drand48();
    B->X[1] = 2 * drand48();
    B->SRC  = drand48();
  }
  C_iter Cj = cells.begin();
  Cj->X = 1;
  Cj->BODY = jbodies.begin();
  Cj->NCBODY = jbodies.size();
  Cj->M = 0;
  kernel.P2M(Cj);

#if 1
  C_iter CJ = cells.begin()+1;
  CJ->CHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 0;
  CJ->M = 0;
  kernel.M2M(CJ,cells.begin());

  C_iter CI = cells.begin()+2;
  CI->X = R + 4;
  CI->M = 1;
  CI->L = 0;
  kernel.M2L(CI,CJ,false);

  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->PARENT = 2;
  Ci->M = 1;
  Ci->L = 0;
  kernel.L2L(Ci,cells.begin());
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->M = 1;
  Ci->L = 0;
  kernel.M2L(Ci,Cj,false);
#endif

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] = R + 2 + 2 * drand48();
    B->X[1] = R + 2 + 2 * drand48();
    B->SRC  = drand48();
    B->TRG  = 0;
  }
  Ci->BODY = bodies.begin();
  Ci->NCBODY = bodies.size();
  kernel.L2P(Ci);

  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    *B = bodies[B-bodies2.begin()];
    B->TRG = 0;
  }
  Cj->NDBODY = jbodies.size();
  Ci->NDBODY = bodies2.size();
  Ci->BODY = bodies2.begin();
  kernel.P2P(Ci,Cj,false);
  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    B->TRG /= B->SRC;
  }

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double diff = 0, norm = 0;
  for (B_iter B=bodies.begin(),B2=bodies2.begin(); B!=bodies.end(); B++,B2++) {
    diff += (B->TRG - B2->TRG) * (B->TRG - B2->TRG);
    norm += B2->TRG * B2->TRG;
  }
  double err = std::sqrt(diff/norm);
  std::cout << P << " " << err << std::endl;
  file << P << " " << err << std::endl;
  file.close();
  return 0;
}
