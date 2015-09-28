#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal2.h"
#include "up_down_pass.h"
#include "verify.h"

int main(int argc, char ** argv) {
  Args args(argc,argv);
  Bodies bodies, bodies2, jbodies, buffer;
  BoundBox boundBox(args.nspawn);
  Bounds bounds;
  BuildTree buildTree(args.ncrit, args.nspawn);
  Cells cells, jcells;
  Dataset data;
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;

  kernel::eps2 = 0.0;
  kernel::wavek = complex_t(10.,1.) / real_t(2 * M_PI);
  kernel::setup();
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, P);
  bodies = data.initBodies(args.numBodies, args.distribution, 0);
  buffer.reserve(bodies.size());
  if (args.IneJ) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->X[0] += M_PI;
      B->X[0] *= 0.5;
    }
    jbodies = data.initBodies(args.numBodies, args.distribution, 1);
    for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
      B->X[0] -= M_PI;
      B->X[0] *= 0.5;
    }
  }
  for (int t=0; t<args.repeat; t++) {
    logger::printTitle("FMM Profiling");
    logger::startTimer("Total FMM");
    logger::startPAPI();
    logger::startDAG();
    bounds = boundBox.getBounds(bodies);
    if (args.IneJ) {
      bounds = boundBox.getBounds(jbodies,bounds);
    }
    cells = buildTree.buildTree(bodies, buffer, bounds);
    upDownPass.upwardPass(cells);
    listBasedTraversal(cells);
    upDownPass.downwardPass(cells);
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
    kernel::Xperiodic = 0;
    bool mutual = false;
    kernel::P2P(Ci, Cj, mutual);
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
