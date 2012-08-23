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
#ifdef MANY
  for( int it=0; it<25; it++ ) {
    int numBodies = int(pow(10,(it+24)/8.0));
#else
  {
    int numBodies = ARGS.numBodies;
#endif // MANY
    if(FMM.printNow) std::cout << std::endl
      << "N                    : " << numBodies << std::endl;
    bodies.resize(numBodies);
    DATA.initBodies(bodies, ARGS.distribution);
    FMM.setBounds(bodies);
    FMM.buildTree(bodies, cells);                               // TODO : make it work without this
    FMM.resetTimer();
    FMM.startTimer("FMM");
    FMM.buildTree(bodies, cells);
    FMM.upwardPass(cells);
    FMM.startPAPI();
    FMM.evaluate(cells,cells,ARGS.mutual);
    FMM.stopPAPI();
    FMM.downwardPass(cells);
    FMM.stopTimer("FMM", FMM.printNow);
    FMM.eraseTimer("FMM");
    FMM.writeTime();
    FMM.resetTimer();
    jbodies = bodies;
    if (int(bodies.size()) > ARGS.numTarget) bodies.resize(ARGS.numTarget);
    Bodies bodies2 = bodies;
    DATA.initTarget(bodies2);
    FMM.startTimer("Direct sum");
    FMM.direct(bodies2, jbodies);
    FMM.normalize(bodies2);
    FMM.stopTimer("Direct sum",FMM.printNow);
    FMM.eraseTimer("Direct sum");
    double diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    DATA.evalError(bodies, bodies2, diff1, norm1, diff2, norm2);
    if(FMM.printNow) DATA.printError(diff1, norm1, diff2, norm2);
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
}
