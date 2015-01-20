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
  const complex_t wavek(10.,1.),imag(0.,1.);
  vec3 * Xj = new vec3 [numBodies];
  complex_t * qj = new complex_t [numBodies];
  complex_t * pi = new complex_t [numBodies];
  cvec3 * Fi = new cvec3 [numBodies];
  complex_t * pi2 = new complex_t [numBodies];
  cvec3 * Fi2 = new cvec3 [numBodies];
  logger::verbose = args.verbose;
  FILE * fid = fopen("data","r");
  for (int i=0; i<numBodies; i++) {
    if (fscanf(fid,"%lf %lf %lf",&Xj[i][0],&Xj[i][1],&Xj[i][2])) {
      qj[i] = Xj[i][0] + imag * Xj[i][1];
    }
  }
  fclose(fid);
  logger::startTimer("FMM");
  fmm(wavek,numBodies,Xj,qj,pi,Fi);
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
  P2P(icell, pi2, Fi2, jcell, Xj, qj, wavek);
  logger::stopTimer("Direct");
  real_t potDif = 0, potNrm = 0, accDif = 0, accNrm = 0;
  for (int i=0; i<numTarget; i++) {
    potDif += abs(pi[i] - pi2[i]) * abs(pi[i] - pi2[i]);
    potNrm += abs(pi2[i]) * abs(pi2[i]);
    accDif += abs(Fi[i][0] - Fi2[i][0]) * abs(Fi[i][0] - Fi2[i][0])
      + abs(Fi[i][1] - Fi2[i][1]) * abs(Fi[i][1] - Fi2[i][1])
      + abs(Fi[i][2] - Fi2[i][2]) * abs(Fi[i][2] - Fi2[i][2]);
    accNrm += abs(Fi2[i][0]) * abs(Fi2[i][0])
      + abs(Fi2[i][1]) * abs(Fi2[i][1])
      + abs(Fi2[i][2]) * abs(Fi2[i][2]);
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
