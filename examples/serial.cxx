#include "args.h"
#include "dataset.h"
#include "serialfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args ARGS(argc, argv);
  Bodies bodies, jbodies;
  Cells cells, jcells;
  Dataset DATA;
  SerialFMM FMM;
  FMM.NCRIT = ARGS.NCRIT;
  FMM.NSPAWN = ARGS.NSPAWN;
  FMM.IMAGES = ARGS.IMAGES;
  FMM.THETA = ARGS.THETA;
  FMM.printNow = true;
#if AUTO
  FMM.timeKernels();
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
    if(FMM.printNow) std::cout << std::endl
      << "Num bodies           : " << numBodies << std::endl;
    bodies.resize(numBodies);
    DATA.initBodies(bodies, ARGS.distribution);
    Box box = FMM.setBounds(bodies);
    FMM.buildTree(bodies, cells, box);                          //TODO : make it work without this
    FMM.resetTimer();
    FMM.startTimer("Total FMM");
    FMM.buildTree(bodies, cells, box);
    FMM.upwardPass(cells);
    FMM.startPAPI();
    FMM.dualTreeTraversal(cells, cells, FMM.periodicCycle, ARGS.mutual);
    FMM.stopPAPI();
    FMM.downwardPass(cells);
    std::cout << "----------------------------------" << std::endl;
    FMM.stopTimer("Total FMM", FMM.printNow);
    FMM.eraseTimer("Total FMM");
    FMM.writeTime();
    FMM.resetTimer();
    jbodies = bodies;
    if (int(bodies.size()) > ARGS.numTarget) DATA.sampleBodies(bodies, ARGS.numTarget);
    Bodies bodies2 = bodies;
    DATA.initTarget(bodies2);
    FMM.startTimer("Total Direct");
    FMM.direct(bodies2, jbodies);
    FMM.normalize(bodies2);
    std::cout << "----------------------------------" << std::endl;
    FMM.stopTimer("Total Direct", FMM.printNow);
    FMM.eraseTimer("Total Direct");
    std::cout << "----------------------------------" << std::endl;
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    DATA.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    if(FMM.printNow) {
      DATA.printError(diff1, norm1, diff2, norm2);
      FMM.printTreeData(cells);
      FMM.printTraversalData();
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
