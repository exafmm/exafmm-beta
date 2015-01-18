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
  delete[] Xj;
  delete[] qj;
  delete[] pi;
  delete[] Fi;
  delete[] pi2;
  delete[] Fi2;
}
