#include "base_mpi.h"
#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "ewald.h"
#include "traversal.h"
#include "tree_mpi.h"
#include "up_down_pass.h"
#include "verify.h"
#include "serialfmm.h"

int main(int argc, char ** argv) {
  Args args(argc, argv);
  BaseMPI baseMPI;
  BoundBox boundBox(args.nspawn);
  Dataset data;
  BuildTree buildTree(args.ncrit, args.nspawn);
  Traversal traversal(args.nspawn, args.images);
  UpDownPass upDownPass(args.theta, args.useRmax, args.useRopt);
  Verify verify;

  const int numBodies = args.numBodies;
  const int ksize = 11;
  const real cycle = 2 * M_PI;
  const real_t alpha = 10 / cycle;
  const real_t sigma = .25 / M_PI;
  const real_t cutoff = cycle * alpha / 3;
  Ewald ewald(ksize, alpha, sigma, cutoff, cycle);
  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength, PP);

  logger::printTitle("FMM Profiling");
  logger::startTimer("Total FMM");
  for( int it=0; it<1; it++ ) {
    srand48(0);
    Bodies bodies = data.initBodies(args.numBodies, args.distribution, baseMPI.mpirank, baseMPI.mpisize);
    
    Bodies jbodies = bodies;

    logger::startTimer("Total Direct");
    const int numTargets = numBodies;
    data.sampleBodies(bodies, numTargets);
    data.initTarget(bodies);
    traversal.direct(bodies, jbodies, cycle);
    traversal.normalize(bodies);
    vec3 localDipole = upDownPass.getDipole(bodies, 0);
    vec3 globalDipole = baseMPI.allreduceVec3(localDipole);
    int numBodies = baseMPI.allreduceInt(bodies.size());
    upDownPass.dipoleCorrection(bodies, globalDipole, numBodies, cycle);
    logger::stopTimer("Total Direct");

    logger::startTimer("Total Ewald");
    Bounds bounds = boundBox.getBounds(bodies);
    Cells cells = buildTree.buildTree(bodies, bounds);
    Bodies bodies2 = bodies;
    data.initTarget(bodies);
    bounds = boundBox.getBounds(jbodies);
    Cells jcells = buildTree.buildTree(jbodies, bounds);
    ewald.wavePart(bodies, jbodies);
    ewald.realPart(cells, jcells);
    double potSum = verify.getSumScalar(bodies);
    double potSum2 = verify.getSumScalar(bodies2);
    double accDif = verify.getDifVector(bodies, bodies2);
    double accNrm = verify.getNrmVector(bodies);
    logger::printTitle("FMM vs. direct");
    double potDif = (potSum - potSum2) * (potSum - potSum2);
    double potNrm = potSum * potSum;
    verify.print("Rel. L2 Error (pot)",std::sqrt(potDif/potNrm));
    verify.print("Rel. L2 Error (acc)",std::sqrt(accDif/accNrm));
  }
}
