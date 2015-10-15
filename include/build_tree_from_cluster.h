#ifndef build_tree_from_cluster_h
#define build_tree_from_cluster_h
#include "build_tree.h"

namespace exafmm {
  class BuildTreeFromCluster {
  public:
    typedef std::vector<int> ints;
    typedef std::vector<vec3> vec3s;
    ints iwrap;

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
      iwrap.resize(bodies.size());
      int numCells = getNumCells(bodies);
      vec3s Xmin = getXmin(bodies, numCells);
      Bodies cluster(numCells);
      int icell = bodies.begin()->ICELL;
      int numBodies = 0;
      B_iter C=cluster.begin();
      int b = 0, c = 0;
      C->X = 0;
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
	int index = B->ICELL;
	if (index != icell) {
	  assert(index > icell);
	  C->X /= numBodies;
	  C->IBODY = B - bodies.begin() - numBodies;
	  C->ICELL = icell;
	  C++;
	  c++;
	  C->X = 0;
	  numBodies = 0;
	  icell = index;
	}
	iwrap[b] = 0;
	for (int d=0; d<3; d++) {
	  int flag = B->X[d] - Xmin[c][d] > cycle / 2;
	  B->X[d] -= cycle * flag;
	  iwrap[b] |= flag << d;
	}
	C->X += B->X;
	numBodies++;
      }
      C->X /= numBodies;
      C->IBODY = bodies.size() - numBodies;
      C->ICELL = icell;
      return cluster;
    }

    void upwardPass(C_iter C, C_iter C0) {
      for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) {
	upwardPass(CC, C0);
	C->NBODY += CC->NBODY;
	C->R = std::max(C->R, 2 * CC->R);
	C->X += CC->X;
      }
      if (C->NCHILD != 0) C->X /= C->NCHILD;
    }

    void attachClusterBodies(Bodies & bodies, Cells & cells, real_t cycle) {
      B_iter B0 = bodies.begin();
      int numLeafs = 0;
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {
	if (C->NCHILD == 0) numLeafs++;
      }
      int dimLeafs = cbrt(numLeafs);
      assert(dimLeafs*dimLeafs*dimLeafs == numLeafs);
      int numLevels = log2(dimLeafs) + 1;
      for (C_iter C=cells.begin(); C!=cells.end(); C++) {
	if (C->NCHILD == 0) {
	  uint64_t icell = C->BODY->ICELL;
	  B_iter B = B0 + C->BODY->IBODY;
	  int ix = icell / dimLeafs / dimLeafs;
	  int iy = icell / dimLeafs % dimLeafs;
	  int iz = icell % dimLeafs;
	  int key = 0;
	  for (int l=0; l<numLevels; l++) {
	    key += (ix & 1) << (3 * l);
	    key += (iy & 1) << (3 * l + 1);
	    key += (iz & 1) << (3 * l + 2);
	    ix >>= 1;
	    iy >>= 1;
	    iz >>= 1;
	  }
	  C->X = C->BODY->X;
	  C->BODY = B;
	  C->IBODY = B - B0;
	  C->NBODY = 0;
	  C->ICELL = key;
	  C->R = cycle / dimLeafs / 2;
	  while (B->ICELL == icell) {
	    C->NBODY++;
	    if (B==bodies.end()-1) break;
	    B++;
	  }
	} else {
	  C->NBODY = 0;
	  C->ICELL = 0;
	  C->R = 0;
	  C->X = 0;
	}
      }
      upwardPass(cells.begin(), cells.begin());
    }

    void shiftBackBodies(Bodies & bodies, real_t cycle) {
      int b = 0;
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
	unwrap(B->X,cycle,iwrap[b]);
      }
    }
  };
}
#endif
