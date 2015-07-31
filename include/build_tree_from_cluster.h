#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

class BuildTreeFromCluster {
private:
  const int mask;
  const int nspawn;
  
public:
  BuildTreeFromCluster(int _mask, int _nspawn) : mask(_mask), nspawn(_nspawn) {}

  Bodies setClusterCenter(Bodies & bodies) {
    int ibody = -1;
    int numCells = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->IBODY & mask;
      if (index != ibody) {
	numCells++;
	ibody = index;
      }
    }
    Bodies cluster(numCells);
    ibody = bodies.begin()->IBODY & mask;
    int numBodies = 0;
    B_iter C=cluster.begin();
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int index = B->IBODY & mask;
      if (index != ibody) {
	C->X /= numBodies;
	C->IBODY = ibody;
	C++;
	numBodies = 0;
	ibody = index;
      }
      C->X += B->X;
      C->IBODY = ibody;
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
	C->IBODY = B - B0;
	C->NBODY = 0;
	while (B->IBODY == C->BODY->IBODY) {
	  B++;
	  C->NBODY++;
	}
      }
    }
  }
};
#endif
