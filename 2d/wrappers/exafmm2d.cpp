#include "mex.h"

#include "args.h"
#include "boundbox.h"
#include "buildtree.h"
#include "dataset.h"
#include "logger.h"
#include "traversal.h"
#include "updownpass.h"

void mexFunction (int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[]) {
  mxArray *xj_in_m, *yj_in_m, *sj_in_m;                         // Source info
  mxArray *xi_in_m, *yi_in_m, *si_in_m;                         // Target info
  mxArray *ti_out_m;                                            // Target result
  double *xj, *yj, *sj;
  double *xi, *yi, *si;
  double *ti;
  const mwSize *dims;                                           // Array dimensions
  int dimx, dimy;
  int nsrc, ntrg;                                               // Number of bodies

  xj_in_m = mxDuplicateArray(prhs[0]);                          // Get source info
  yj_in_m = mxDuplicateArray(prhs[1]);
  sj_in_m = mxDuplicateArray(prhs[2]);
  xi_in_m = mxDuplicateArray(prhs[3]);                          // Get target info
  yi_in_m = mxDuplicateArray(prhs[4]);
  si_in_m = mxDuplicateArray(prhs[5]);
  
  dims = mxGetDimensions(prhs[0]);                              // Get number of sources
  nsrc = (int)dims[0]*(int)dims[1];

  dims = mxGetDimensions(prhs[3]);                              // Get number of targets
  dimy = (int)dims[0]; dimx = (int)dims[1];
  ntrg = (int)dims[0]*(int)dims[1];

  ti_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);  // Create output array
  
  xj = mxGetPr(xj_in_m);                                        // Associate pointers
  yj = mxGetPr(yj_in_m);
  sj = mxGetPr(sj_in_m);
  xi = mxGetPr(xi_in_m);
  yi = mxGetPr(yi_in_m);
  si = mxGetPr(si_in_m);
  ti = mxGetPr(ti_out_m);
  
  int nspawn = 1000;                                            // ExaFMM setup
  int ncrit = 64;
  real_t theta = 0.4;
  int images = 0;
  bool verbose = true;
  Logger logger;
  const real_t eps2 = 0.0;
  const real_t cycle = 2 * M_PI;
  BoundBox boundbox(nspawn);
  BuildTree tree(ncrit,nspawn);
  UpDownPass pass(theta,eps2);
  Traversal traversal(nspawn,images,eps2);
  if (verbose) {
    logger.verbose = true;
    boundbox.verbose = true;
    tree.verbose = true;
    pass.verbose = true;
    traversal.verbose = true;
  }
  logger.printTitle("FMM Profiling");
  logger.startTimer("Total FMM");
  logger.startPAPI();

  Bodies bodies(ntrg);                                          // Define targets
  for(int i=0;i<ntrg;i++) {
    bodies[i].X[0] = xi[i];
    bodies[i].X[1] = yi[i];
    bodies[i].SRC = si[i];
    bodies[i].TRG = 0;
    bodies[i].IBODY = i;
  }
  Bodies jbodies(nsrc);                                         // Define sources
  for(int i=0;i<nsrc;i++) {
    jbodies[i].X[0] = xj[i];
    jbodies[i].X[1] = yj[i];
    jbodies[i].SRC = sj[i];
    jbodies[i].TRG = 0;
    jbodies[i].IBODY = i;
  }

  Bounds bounds = boundbox.getBounds(bodies);                   // Apply ExaFMM
  bounds = boundbox.getBounds(jbodies,bounds);
  Cells cells = tree.buildTree(bodies, bounds);
  pass.upwardPass(cells);
  Cells jcells = tree.buildTree(jbodies, bounds);
  pass.upwardPass(jcells);
  traversal.dualTreeTraversal(cells, jcells, cycle);
  pass.downwardPass(cells);

  logger.printTitle("Total runtime");
  logger.stopPAPI();
  logger.stopTimer("Total FMM");
  boundbox.writeTime();
  tree.writeTime();
  pass.writeTime();
  traversal.writeTime();
  tree.printTreeData(cells);
  traversal.printTraversalData();
  logger.printPAPI();

  for(int i=0; i<ntrg; i++) {                                   // Loop over targets
    ti[bodies[i].IBODY] = bodies[i].TRG;                        //  Recover results
  }                                                             // End loop
  return;
}
