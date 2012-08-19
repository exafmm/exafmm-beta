#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <errno.h>

#include "dataset.h"
#include "serialfmm.h"

static int gendata(Dataset& data, Bodies& bodies, const char * distribution) {
  switch (distribution[0]) {
  case 'l':
    data.lattice(bodies);
    return 1;
  case 'c':
    data.cube(bodies);
    return 1;
  case 's':
    data.sphere(bodies);
    return 1;
#if PLUMMER_DISTRIBUTION
  case 'p':
    data.plummer(bodies);
    return 1;			// OK
#endif
  default:
    fprintf(stderr, "unknown data distribution %s\n", distribution);
    return 0;			// NG
  }
}

int main(int argc, char ** argv) {
  exafmm_config o[1];
  if (parse_cmdline_args(argc, argv, o) == NULL) {
    return EXIT_FAILURE;
  }
  show_exafmm_config(o);
  IMAGES = o->images;
  THETA = o->theta;
  NCRIT = o->ncrit;
#if SIMDIZATION
  SIMDIZE = o->simdize;
#endif
#if IMPL_MUTUAL
  splitBothThreshold = o->splitBothThreshold;
#endif

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
#if HYBRID
    FMM.timeKernels();
#endif
    if(FMM.printNow) std::cout << "N                    : " << o->numBodies << std::endl;
    bodies.resize(o->numBodies);
    gendata(DATA, bodies, o->distribution);
#if PARALLEL_EVERYTHING
    if (o->parallelEverything) t_bodies.resize(bodies.size());
#endif
    for (int i = 0; i < o->steps; i++) {
      FMM.startTimer("FMM");
#if PARALLEL_EVERYTHING
      if (o->parallelEverything) {
	DATA.initTargetRec(bodies);
	FMM.setBoundsRec(bodies);
	FMM.buildTreeRec(bodies,t_bodies,cells);
	FMM.upwardPassRec(cells);
      } else {
	DATA.initTarget(bodies);
	FMM.setBounds(bodies);
	FMM.buildTree(bodies,cells);
	FMM.upwardPass(cells);
      }
#else
      FMM.setBounds(bodies);
      FMM.buildTree(bodies,cells);
      FMM.upwardPass(cells);
#endif

      if (o->buildOnly == 0) {
	FMM.startPAPI();
#if IneJ
#if IMPL_MUTUAL
	FMM.evaluate(cells,cells,o->mutual);
#else
	FMM.evaluate(cells,cells);
#endif
#else
	FMM.evaluate(cells);
#endif
	FMM.stopPAPI();
#if PARALLEL_EVERYTHING
	if (o->parallelEverything)
	  FMM.downwardPassRec(cells);
	else
	  FMM.downwardPass(cells);
#else
	FMM.downwardPass(cells);
#endif
      }
      FMM.stopTimer("FMM",FMM.printNow);
      FMM.eraseTimer("FMM");
      FMM.writeTime();
      FMM.resetTimer();
    }
    
    if (o->buildOnly == 0 && o->evalError) {
      jbodies = bodies;
      if (bodies.size() > o->evalError) bodies.resize(o->evalError);
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
  }
}
