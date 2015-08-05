#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

class BuildTreeFromCluster {
public:
  BuildTreeFromCluster() {}

  Bodies setClusterCenter(Bodies & bodies) {
    int icell = -1;
    int numCells = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->ICELL;
      if (index != icell) {
	numCells++;
	icell = index;
      }
    }
    Bodies cluster(numCells);
    icell = bodies.begin()->ICELL;
    int numBodies = 0;
    B_iter C=cluster.begin();
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->ICELL;
      if (index != icell) {
	C->X /= numBodies;
	C->ICELL = icell;
	C++;
	numBodies = 0;
	icell = index;
      }
      C->X += B->X;
      C->ICELL = icell;
      numBodies++;
    }
    C->X /= numBodies;
    return cluster;
  }

  void attachClusterBodies(Bodies & bodies, Cells & cells) {
    B_iter B0 = bodies.begin();
    B_iter B = B0;
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {
      if (C->NCHILD == 0) {
	C->BODY = B;
	C->ICELL = B - B0;
	C->NBODY = 0;
	while (B->ICELL == C->BODY->ICELL) {
	  B++;
	  C->NBODY++;
	}
      }
    }
  }
};
#endif
