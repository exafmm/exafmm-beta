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
  vec3 * Xj = new vec3 [numBodies];
  real_t * qj = new real_t [numBodies];
  real_t * pi = new real_t [numBodies];
  vec3 * Fi = new vec3 [numBodies];
  real_t * pi2 = new real_t [numBodies];
  vec3 * Fi2 = new vec3 [numBodies];
  logger::verbose = args.verbose;
  for (int i=0; i<numBodies; i++) {
    Xj[i][0] = drand48();
    Xj[i][1] = drand48();
    Xj[i][2] = drand48();
    qj[i] = drand48();
  }
  logger::startTimer("FMM");
  fmm(numBodies, Xj, qj, pi, Fi);
  logger::stopTimer("FMM");
  const int numTarget = 100;
  for (int i=0; i<numTarget; i++) {
    pi2[i] = 0.0;
    Fi2[i] = 0.0;
  }
  int icell[10], jcell[10];
  icell[7] = 0;
  icell[8] = numTarget;
  jcell[7] = 0;
  jcell[8] = numBodies;
  logger::startTimer("Direct");
  P2P(icell, pi2, Fi2, jcell, Xj, qj);
  logger::stopTimer("Direct");
  real_t potDif = 0, potNrm = 0, accDif = 0, accNrm = 0;
  for (int i=0; i<numTarget; i++) {
    potDif += (pi[i] - pi2[i]) * (pi[i] - pi2[i]);
    potNrm += pi2[i] * pi2[i];
    accDif += (Fi[i][0] - Fi2[i][0]) * (Fi[i][0] - Fi2[i][0])
      + (Fi[i][1] - Fi2[i][1]) * (Fi[i][1] - Fi2[i][1])
      + (Fi[i][2] - Fi2[i][2]) * (Fi[i][2] - Fi2[i][2]);
    accNrm += Fi2[i][0] * Fi2[i][0]
      + Fi2[i][1] * Fi2[i][1]
      + Fi2[i][2] * Fi2[i][2];
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
