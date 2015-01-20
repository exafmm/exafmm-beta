#include "arrays.h"
#include "args.h"
#include "bound_box.h"
#include "constants.h"
#include "tree.h"
#include "translate.h"
#include "kernel.h"
#include "fmm.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  const int numBodies=args.numBodies;
  const complex_t wavek(10.,1.),imag(0.,1.);
  vec3 * Xj = new vec3 [numBodies];
  complex_t * qj = new complex_t [numBodies];
  complex_t * pi = new complex_t [numBodies];
  cvec3 * Fi = new cvec3 [numBodies];
  complex_t * pi2 = new complex_t [numBodies];
  cvec3 * Fi2 = new cvec3 [numBodies];
  FILE * fid = fopen("data","r");
  for (int i=0; i<numBodies; i++) {
    fscanf(fid,"%lf %lf %lf",&Xj[i][0],&Xj[i][1],&Xj[i][2]);
    qj[i] = Xj[i][0] + imag * Xj[i][1];
  }
  fclose(fid);
  fmm(wavek,numBodies,Xj,qj,pi,Fi);
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
  P2P(icell, pi2, Fi2, jcell, Xj, qj, wavek);
  real_t pdiff = 0, pnorm = 0, fdiff = 0, fnorm = 0;
  for (int i=0; i<numTarget; i++) {
    pdiff += abs(pi[i] - pi2[i]) * abs(pi[i] - pi2[i]);
    pnorm += abs(pi2[i]) * abs(pi2[i]);
    fdiff += abs(Fi[i][0] - Fi2[i][0]) * abs(Fi[i][0] - Fi2[i][0])
      + abs(Fi[i][1] - Fi2[i][1]) * abs(Fi[i][1] - Fi2[i][1])
      + abs(Fi[i][2] - Fi2[i][2]) * abs(Fi[i][2] - Fi2[i][2]);
    fnorm += abs(Fi2[i][0]) * abs(Fi2[i][0])
      + abs(Fi2[i][1]) * abs(Fi2[i][1])
      + abs(Fi2[i][2]) * abs(Fi2[i][2]);
  }
  std::cout << "Err pot= " << sqrt(pdiff/pnorm) << std::endl; 
  std::cout << "Err acc= " << sqrt(fdiff/fnorm) << std::endl; 
  delete[] Xj;
  delete[] qj;
  delete[] pi;
  delete[] Fi;
  delete[] pi2;
  delete[] Fi2;
}
