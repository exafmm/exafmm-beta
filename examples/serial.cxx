#include "args.h"
#include "bounds.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args ARGS(argc, argv);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset DAT;
  Logger LOG;
  Bounds BND(ARGS.NSPAWN,ARGS.IMAGES);
  BuildTree BLD(ARGS.NCRIT,ARGS.NSPAWN);
  UpDownPass UDP(ARGS.IMAGES,ARGS.THETA);
  Traversal TRV(ARGS.NSPAWN,ARGS.IMAGES);
  LOG.printNow = true;
  BND.printNow = true;
  BLD.printNow = true;
  UDP.printNow = true;
  TRV.printNow = true;
#if AUTO
  TRV.timeKernels();
#endif
#if _OPENMP
#pragma omp parallel
#pragma omp master
#endif
#ifdef MANY
  for( int it=0; it<25; it++ ) {
    int numBodies = int(pow(10,(it+24)/8.0));
#else
  {
    int numBodies = ARGS.numBodies;
#endif // MANY
    if(LOG.printNow) std::cout << std::endl
      << "Num bodies           : " << numBodies << std::endl;
    std::cout << "--- Profiling --------------------" << std::endl;
    bodies.resize(numBodies);
    DAT.initBodies(bodies, ARGS.distribution);
    BND.setLocal(bodies);
    BLD.buildTree(bodies, cells, BND.localBox);                 //TODO : make it work without this
    BLD.resetTimer();
    LOG.startTimer("Total FMM");
    BLD.buildTree(bodies, cells, BND.localBox);
    UDP.upwardPass(cells);
    LOG.startPAPI();
    TRV.dualTreeTraversal(cells, cells, BND.CYCLE, ARGS.mutual);
    LOG.stopPAPI();
    UDP.downwardPass(cells);
    std::cout << "--- Total runtime ----------------" << std::endl;
    LOG.stopTimer("Total FMM", LOG.printNow);
    BND.writeTime();
    BLD.writeTime();
    UDP.writeTime();
    TRV.writeTime();
    BND.resetTimer();
    BLD.resetTimer();
    UDP.resetTimer();
    TRV.resetTimer();
    jbodies = bodies;
    if (int(bodies.size()) > ARGS.numTarget) DAT.sampleBodies(bodies, ARGS.numTarget);
    Bodies bodies2 = bodies;
    DAT.initTarget(bodies2);
    LOG.startTimer("Total Direct");
    UDP.direct(bodies2, jbodies, BND.CYCLE);
    UDP.normalize(bodies2);
    LOG.stopTimer("Total Direct", LOG.printNow);
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    DAT.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    if(LOG.printNow) {
      DAT.printError(diff1, norm1, diff2, norm2);
      BLD.printTreeData(cells);
      TRV.printTraversalData();
    }
  }
#ifdef VTK
  for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) B->ICELL = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.ICELL = 1;
    body.X     = C->X;
    body.SRC   = 0;
    jbodies.push_back(body);
  }
  int Ncell = 0;
  vtkPlot vtk;
  vtk.setDomain(M_PI,0);
  vtk.setGroupOfPoints(jbodies, Ncell);
  vtk.plot(Ncell);
#endif
  return 0;
}
