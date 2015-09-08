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
  logger::verbose = args.verbose;
  bodies.resize(numBodies);
  B_iter B = bodies.begin();
  for (int i=0; i<numBodies; i++,B++) {
    Xj[i][0] = drand48();
    Xj[i][1] = drand48();
    Xj[i][2] = drand48();
    B->SRC = Xj[i][0] + I * Xj[i][1];
    B->X = Xj[i];
  }
  logger::startTimer("Total FMM");
  fmm(numBodies,Xj);
  logger::stopTimer("Total FMM");
  const int numTarget = 100;
  Bodies bodies2(numTarget);
  for (int i=0; i<numTarget; i++) {
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
  logger::startTimer("Total Direct");
  real_t eps2 = 0;
  vec3 Xperiodic = 0;
  bool mutual = false;
  kernel::P2P(Ci, Cj, eps2, Xperiodic, mutual);
  logger::stopTimer("Total Direct");
  std::complex<double> potDif = verify.getDifScalar(bodies2, bodies);
  std::complex<double> potNrm = verify.getNrmScalar(bodies2);
  std::complex<double> accDif = verify.getDifVector(bodies2, bodies);
  std::complex<double> accNrm = verify.getNrmVector(bodies2);
  logger::printTitle("FMM vs. direct");
  verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
  verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
  delete[] Xj;
  delete[] qj;
}
