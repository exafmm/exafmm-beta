#include "args.h"
#include "bound_box.h"
#include "logger.h"
#include "verify.h"

#include "arrays.h"
#include "constants.h"
#include "tree.h"
#include "translate.h"
#include "kernel.h"
#include "fmm.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  Verify verify;
  const int numBodies=args.numBodies;
  wavek = complex_t(10.,1.);
  vec3 * Xj = new vec3 [numBodies];
  complex_t * qj = new complex_t [numBodies];
  complex_t * pi = new complex_t [numBodies];
  cvec3 * Fi = new cvec3 [numBodies];
  complex_t * pi2 = new complex_t [numBodies];
  cvec3 * Fi2 = new cvec3 [numBodies];
  logger::verbose = args.verbose;
  bodies.resize(numBodies);
  B_iter B = bodies.begin();
  for (int i=0; i<numBodies; i++,B++) {
    Xj[i][0] = drand48();
    Xj[i][1] = drand48();
    Xj[i][2] = drand48();
    qj[i] = Xj[i][0] + I * Xj[i][1];
    B->X = Xj[i];
    B->SRC  = qj[i];
  }
  logger::startTimer("FMM");
  fmm(numBodies,Xj);
  logger::stopTimer("FMM");
  const int numTarget = 100;
  Bodies bodies2(numTarget);
  for (int i=0; i<numTarget; i++) {
    pi2[i] = 0.0;
    Fi2[i] = 0.0;
    bodies2[i] = bodies[i];
    bodies2[i].TRG = 0;
  }
  cells.resize(2);
  C_iter Ci = cells.begin();
  C_iter Cj = cells.begin() + 1;
  Ci->BODY = bodies2.begin();
  Ci->NBODY = bodies2.size();
  Cj->BODY = bodies.begin();
  Cj->NBODY = bodies.size();
  logger::startTimer("Direct");
  P2P(Ci, Cj);
  logger::stopTimer("Direct");
  std::complex<double> potDif = 0, potNrm = 0, accDif = 0, accNrm = 0;
  for (int i=0; i<numTarget; i++) {
    B_iter B = bodies.begin() + i;
    B_iter B2 = bodies2.begin() + i;
    potDif += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);
    potNrm += B2->TRG[0] * B2->TRG[0];
    accDif += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1])
      + (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2])
      + (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);
    accNrm += B2->TRG[1] * B2->TRG[1]
      + B2->TRG[2] * B2->TRG[2]
      + B2->TRG[3] * B2->TRG[3];
  }
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
  delete[] Xj;
  delete[] qj;
  delete[] pi;
  delete[] Fi;
  delete[] pi2;
  delete[] Fi2;
}
