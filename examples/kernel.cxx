#include <fstream>
#include "kernel.h"
#include <vector>
#include "verify.h"
using namespace exafmm;

int main() {
  Bodies bodies(1), bodies2(1), jbodies(1);
  const real_t theta = 0.5;
  const real_t R = 2 / theta;
  kernel::eps2 = 0.0;
  kernel::Xperiodic = 0;

  Cells cells(4);
  Verify verify;
  jbodies[0].X = 0;
#if EXAFMM_BIOTSAVART
  jbodies[0].SRC[0] = drand48();
  jbodies[0].SRC[1] = drand48();
  jbodies[0].SRC[2] = drand48();
  jbodies[0].SRC[3] = 0.1;
#else
  jbodies[0].SRC = 1;
#endif
  C_iter Cj = cells.begin();
  Cj->X = 3. / 16;
  Cj->BODY = jbodies.begin();
  Cj->NBODY = jbodies.size();
  Cj->M = 0;
  kernel::P2M(Cj);

#if 1
  C_iter CJ = cells.begin()+1;
  CJ->ICHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 3. / 8;
  CJ->M = 0;
  kernel::M2M(CJ, cells.begin());

  C_iter CI = cells.begin()+2;
  CI->X = 21. / 8;
  CI->M = 1;
  CI->L = 0;
#if EXAFMM_MASS
  for (int i=1; i<NTERM; i++) CJ->M[i] /= CJ->M[0];
#endif
  kernel::M2L(CI, CJ, false);

  C_iter Ci = cells.begin()+3;
  Ci->X = 45. / 16;
  Ci->IPARENT = 2;
  Ci->M = 1;
  Ci->L = 0;
  kernel::L2L(Ci, cells.begin());
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = 45. / 16;
  Ci->M = 1;
  Ci->L = 0;
#if EXAFMM_MASS
  for (int i=1; i<NTERM; i++) Cj->M[i] /= Cj->M[0];
#endif
  kernel::M2L(Ci, Cj, false);
#endif

  bodies[0].X = 3;
  bodies[0].SRC = 1;
  bodies[0].TRG = 0;
  Ci->BODY = bodies.begin();
  Ci->NBODY = bodies.size();
  kernel::L2P(Ci);

  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    *B = bodies[B-bodies2.begin()];
    B->TRG = 0;
  }
  Cj->NBODY = jbodies.size();
  Ci->NBODY = bodies2.size();
  Ci->BODY = bodies2.begin();
  kernel::P2P(Ci, Cj, false);
  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    B->TRG /= B->SRC;
  }
  std::cout << bodies[0].TRG[0] << " " << bodies2[0].TRG[0] << std::endl;

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double potDif = verify.getDifScalar(bodies, bodies2);
  double potNrm = verify.getNrmScalar(bodies);
  double accDif = verify.getDifVector(bodies, bodies2);
  double accNrm = verify.getNrmVector(bodies);
  std::cout << P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
  file << P << " " << std::sqrt(potDif/potNrm) << "  " << std::sqrt(accDif/accNrm) << std::endl;
  file.close();
  return 0;
}
