#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

class BuildTreeFromCluster {
private:
  const int mask;
  const int nspawn;
  
public:
  BuildTreeFromCluster(int _mask, int _nspawn) : mask(_mask), nspawn(_nspawn) {}

  void setClusterCenter(Bodies bodies) {
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
	numBodies = 0;
	C++;
      }
      C->X += B->X;
      numBodies++;
    }
    C->X /= numBodies;
  }
};
#endif
