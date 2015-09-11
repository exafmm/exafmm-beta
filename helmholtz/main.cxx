#include "args.h"
#include "bound_box.h"
#include "dataset.h"
#include "logger.h"
#include "verify.h"

#include "build_tree.h"
#include "kernel.h"
#include "up_down_pass.h"
#include "traversal.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  Bodies bodies, bodies2, jbodies, buffer;
  BoundBox boundBox(args.nspawn);
  Bounds bounds;
  Dataset data;
  Verify verify;

  kernel::wavek = complex_t(10.,1.) / (2 * M_PI);
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  buffer.reserve(bodies.size());
#if IneJ
  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] += M_PI;
    B->X[0] *= 0.5;
  }
  jbodies = data.initBodies(args.numBodies, args.distribution, 1);
  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    B->X[0] -= M_PI;
    B->X[0] *= 0.5;
  }
#endif
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    logger::startDAG();
    bounds = boundBox.getBounds(bodies);
#if IneJ
    bounds = boundBox.getBounds(jbodies,bounds);
#endif
    Cells cells = buildTree(bodies, buffer, bounds);
    upwardPass(cells);
    evaluate(cells);
    downwardPass(cells);
    logger::printTitle("Total runtime");
    logger::stopDAG();
    logger::stopPAPI();
    logger::stopTimer("Total FMM");
    jbodies = bodies;
    const int numTargets = 100;
    data.sampleBodies(bodies, numTargets);
    bodies2 = bodies;
    data.initTarget(bodies);
    cells.resize(2);
    C_iter Ci = cells.begin();
    C_iter Cj = cells.begin() + 1;
    Ci->BODY = bodies.begin();
    Ci->NBODY = bodies.size();
    Cj->BODY = jbodies.begin();
    Cj->NBODY = jbodies.size();
    logger::startTimer("Total Direct");
    real_t eps2 = 0;
    vec3 Xperiodic = 0;
    bool mutual = false;
    kernel::P2P(Ci, Cj, eps2, Xperiodic, mutual);
    logger::stopTimer("Total Direct");
    std::complex<double> potDif = verify.getDifScalar(bodies, bodies2);
    std::complex<double> potNrm = verify.getNrmScalar(bodies);
    std::complex<double> accDif = verify.getDifVector(bodies, bodies2);
    std::complex<double> accNrm = verify.getNrmVector(bodies);
    logger::printTitle("FMM vs. direct");
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
  }
}
