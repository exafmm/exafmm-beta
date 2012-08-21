#include "dataset.h"
#include "options.h"
#include "serialfmm.h"
#ifdef VTK
#include "vtk.h"
#endif

int main(int argc, char ** argv) {
  Args args[1];
  parse_cmdline_args(argc, argv, args);
  showArgs(args);
  NCRIT = args->ncrit;
  NSPAWN = args->nspawn;
  IMAGES = args->images;
  THETA = args->theta;

#if _OPENMP
#pragma omp parallel
#pragma omp single
#endif
  {
    Bodies bodies, jbodies;
#if PARALLEL_EVERYTHING
    Bodies t_bodies;
#endif
    Cells cells, jcells;
    Dataset DATA;
    SerialFMM FMM;
    FMM.printNow = true;
#if AUTO
    FMM.timeKernels();
#endif
    if(FMM.printNow) std::cout << "N                    : " << args->numBodies << std::endl;
    bodies.resize(args->numBodies);
    gendata(DATA, bodies, args->distribution);
#if PARALLEL_EVERYTHING
    t_bodies.resize(bodies.size());
#endif
    for( int i=0; i<2; i++ ) {
      FMM.startTimer("FMM");
#if PARALLEL_EVERYTHING
      DATA.initTargetRec(bodies);
      FMM.setBoundsRec(bodies);
      FMM.buildTreeRec(bodies,t_bodies,cells);
      FMM.upwardPassRec(cells);
#else
      FMM.setBounds(bodies);
      FMM.buildTree(bodies,cells);
      FMM.upwardPass(cells);
#endif

      if (args->buildOnly == 0) {
        FMM.startPAPI();
        FMM.evaluate(cells,cells,args->mutual);
        FMM.stopPAPI();
#if 1
        FMM.downwardPassRec(cells);
#else
        FMM.downwardPass(cells);
#endif
      }
      FMM.stopTimer("FMM",FMM.printNow);
      FMM.eraseTimer("FMM");
      FMM.writeTime();
      FMM.resetTimer();
    }
    
    if (!args->buildOnly) {
      jbodies = bodies;
      if (int(bodies.size()) > args->numTarget) bodies.resize(args->numTarget);
      Bodies bodies2 = bodies;
      DATA.initTarget(bodies2);
      FMM.startTimer("Direct sum");
      FMM.direct(bodies2,jbodies);
      FMM.normalize(bodies2);
      FMM.stopTimer("Direct sum",FMM.printNow);
      FMM.eraseTimer("Direct sum");
      real_t diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
      DATA.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
      if(FMM.printNow) DATA.printError(diff1,norm1,diff2,norm2);
    }
#ifdef VTK
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) B->ICELL = 0;
    int Ncell = 0;
    vtkPlot vtk;
    vtk.setDomain(M_PI,0);
    vtk.setGroupOfPoints(jbodies,Ncell);
    vtk.plot(Ncell);
#endif
  }
}
