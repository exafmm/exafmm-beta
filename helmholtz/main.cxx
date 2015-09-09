#include "args.h"
#include "bound_box.h"
#include "dataset.h"
#include "logger.h"
#include "verify.h"

#include "arrays.h"
#include "constants.h"
#include "tree.h"
#include "kernel.h"
#include "fmm.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  Bodies bodies, bodies2;
  Dataset data;
  Verify verify;
  const int numBodies=args.numBodies;
  wavek = complex_t(10.,1.) / (2 * M_PI);
  vec3 * Xj = new vec3 [numBodies];
  logger::verbose = args.verbose;
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  bodies.resize(numBodies);
  B_iter B = bodies.begin();
  for (int i=0; i<numBodies; i++,B++) {
    Xj[i] = B->X;
  }
  logger::startTimer("Total FMM");
  int * permutation = new int [numBodies];
  levelOffset = new int [maxLevel];
  vec3 X0;
  real_t R0;
  logger::startTimer("Tree");
  getBounds(Xj, numBodies, X0, R0);
  int numCells, numLevels;
  Cells cells = buildTree(Xj, numBodies, numCells, permutation, numLevels, X0, R0);
  Bodies buffer(numBodies);
  for (int i=0; i<numBodies; i++) {
    buffer[i] = bodies[permutation[i]];
  }
  B = buffer.begin();
  for (C_iter C=cells.begin(); C!=cells.end(); C++) {
    C->BODY = B + C->IBODY;
  }
  logger::stopTimer("Tree");
  evaluate(cells);
  for (int i=0; i<numBodies; i++) {
    bodies[permutation[i]].TRG = buffer[i].TRG;
  }
  delete[] listOffset;
  delete[] lists;
  delete[] levelOffset;
  delete[] permutation;
  logger::stopTimer("Total FMM");
  const int numTarget = 100;
  bodies2.resize(numTarget);
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
}
