#include "args.h"
#include "bound_box.h"
#include "build_tree.h"
#include "dataset.h"
#include "logger.h"
using namespace exafmm;
#include "laplace_cartesian_cpu.h"

int main(int argc, char ** argv) {
  Args args(argc, argv);
  typedef LaplaceCartesianCPU<4,0> Kernel;
  typedef typename Kernel::Bodies Bodies;                       //!< Vector of bodies
  typedef typename Kernel::Cells Cells;                         //!< Vector of cells
  typedef typename Kernel::B_iter B_iter;                       //!< Iterator of body vector
  typedef typename Kernel::C_iter C_iter;                       //!< Iterator of cell vector

  Bodies bodies, bodies2, jbodies, buffer;
  BoundBox<Kernel> boundBox(args.nspawn);
  Bounds bounds;
  BuildTree<Kernel> buildTree(args.ncrit, args.nspawn);
  Cells cells, jcells;
  Dataset<Kernel> data;
  num_threads(args.threads);

  logger::verbose = args.verbose;
  logger::printTitle("FMM Parameters");
  args.print(logger::stringLength);
  buffer.reserve(args.numBodies);
  double * grow1 = new double [args.repeat+1];
  double * link1 = new double [args.repeat+1];
  double * grow2 = new double [args.repeat+1];
  double * link2 = new double [args.repeat+1];
  for (int t=0; t<args.repeat+1; t++) {
    std::cout << t << std::endl;
    bodies = data.initBodies(args.numBodies, args.distribution, 0);
    bounds = boundBox.getBounds(bodies);
    cells = buildTree.buildTree(bodies, buffer, bounds);
    grow1[t] = logger::timer["Grow tree"];
    link1[t] = logger::timer["Link tree"];
    logger::resetTimer();
    cells = buildTree.buildTree(bodies, buffer, bounds);
    grow2[t] = logger::timer["Grow tree"];
    link2[t] = logger::timer["Link tree"];
    logger::resetTimer();
  }
  double grow1ave = 0, link1ave = 0, grow2ave = 0, link2ave = 0;
  for (int t=0; t<args.repeat; t++) {
    grow1ave += grow1[t+1];
    link1ave += link1[t+1];
    grow2ave += grow2[t+1];
    link2ave += link2[t+1];
  }
  grow1ave /= args.repeat;
  link1ave /= args.repeat;
  grow2ave /= args.repeat;
  link2ave /= args.repeat;
  double grow1std = 0, link1std = 0, grow2std = 0, link2std = 0;
  for (int t=0; t<args.repeat; t++) {
    grow1std += (grow1[t+1] - grow1ave) * (grow1[t+1] - grow1ave);
    link1std += (link1[t+1] - link1ave) * (link1[t+1] - link1ave);
    grow2std += (grow2[t+1] - grow2ave) * (grow2[t+1] - grow2ave);
    link2std += (link2[t+1] - link2ave) * (link2[t+1] - link2ave);
  }
  grow1std /= args.repeat;
  link1std /= args.repeat;
  grow2std /= args.repeat;
  link2std /= args.repeat;
  std::cout << "Grow1: " << grow1ave << "+-" << std::sqrt(grow1std)
	    << " Link1: " << link1ave << "+-" << link1std << std::endl;
  std::cout << "Grow2: " << grow2ave << "+-" << std::sqrt(grow2std)
	    << " Link2: " << link2ave << "+-" << link2std << std::endl;
  std::ofstream fid("time.dat", std::ios::app);
  fid << args.numBodies << " " << args.threads << " " << grow1ave << " " << grow2ave << std::endl;
  delete[] grow1;
  delete[] link1;
  delete[] grow2;
  delete[] link2;
  fid.close();
  return 0;
}
