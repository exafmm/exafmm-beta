#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

class BuildTreeFromCluster {
public:
  typedef std::vector<vec3> vec3s;

  BuildTreeFromCluster() {}

  int getNumCells(Bodies & bodies) {
    int icell = -1;
    int numCells = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->ICELL;
      if (index != icell) {
        numCells++;
        icell = index;
      }
    }
    return numCells;
  }

  vec3s getXmin(Bodies & bodies, int numCells) {
    vec3s Xmin(numCells);
    int i = 0;
    int icell = bodies.begin()->ICELL;
    Xmin[0] = bodies.begin()->X;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->ICELL;
      if (index != icell) {
        i++;
        icell = index;
	Xmin[i] = B->X;
      }
      Xmin[i] = min(Xmin[i], B->X);
    }
    return Xmin;
  }

  Bodies setClusterCenter(Bodies & bodies, real_t cycle) {
    int numCells = getNumCells(bodies);
    vec3s Xmin = getXmin(bodies, numCells);
    Bodies cluster(numCells);
    int icell = bodies.begin()->ICELL;
    int numBodies = 0;
    B_iter C=cluster.begin();
    int c = 0;
    C->X = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->ICELL;
      if (index != icell) {
	assert(index > icell);
	C->X /= numBodies;
	C->ICELL = icell;
	C++;
	c++;
	C->X = 0;
	numBodies = 0;
	icell = index;
      }
      vec3 Xperiodic = 0;
      for (int d=0; d<3; d++) {
	Xperiodic[d] = cycle * (B->X[d] - Xmin[c][d] > cycle / 2);
      }
      C->X += B->X - Xperiodic;
      numBodies++;
    }
    C->X /= numBodies;
    C->ICELL = icell;
    return cluster;
  }

  void upwardPass(C_iter C, C_iter C0) {
    for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {
      upwardPass(CC, C0);
      C->NBODY += CC->NBODY;
      C->R = std::max(C->R, 2 * CC->R);
    }
  }

  void attachClusterBodies(Bodies & bodies, Cells & cells) {
    B_iter B0 = bodies.begin();
    B_iter B = B0;
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {
      if (C->NCHILD == 0) {
	C->BODY = B;
	C->IBODY = B - B0;
	C->ICELL = B - B0;
	C->NBODY = 0;
	C->R = 0;
	while (B->ICELL == C->BODY->ICELL) {
	  C->R = std::max(C->R, sqrtf(norm(B->X - C->X)));
	  B++;
	  C->NBODY++;
	}
      } else {
	C->NBODY = 0;
	C->R = 0;
      }
    }
    upwardPass(cells.begin(), cells.begin());
  }
};
#endif
