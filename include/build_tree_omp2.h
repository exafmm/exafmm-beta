#ifndef build_tree_omp2_h
#define build_tree_omp2_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class BuildTree {
private:
  const int ncrit;
  int numLevels;

private:
  void reorder(Box box, int level, int * iX, vec3 * Xj,
	       int * permutation, int n, int * iwork, int * nbody) {
    int offset[9];
    vec3 X;
    real_t R = box.R / (1 << level);
    for (int d=0; d<3; d++) {
      X[d] = box.X[d] - box.R + iX[d] * R * 2 + R;
    }
    for (int i=0; i<8; i++) nbody[i] = 0;
    for (int i=0; i<n; i++) {
      int j = permutation[i];
      int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);
      nbody[octant]++;
    }
    offset[0] = 0;
    for (int i=0; i<8; i++) {
      offset[i+1] = offset[i] + nbody[i];
      nbody[i] = 0;
    }
    for (int i=0; i<n; i++) {
      int j = permutation[i];
      int octant = (Xj[j][2] > X[2]) * 4 + (Xj[j][1] > X[1]) * 2 + (Xj[j][0] > X[0]);
      iwork[offset[octant]+nbody[octant]] = permutation[i];
      nbody[octant]++;
    }
    for (int i=0; i<n; i++) {
      permutation[i] = iwork[i];
    }
  }

  uint64_t getKey(ivec3 iX, int level) {
    uint64_t index = ((1 << 3 * level) - 1) / 7;
    for (int l=0; l<level; l++) {
      for (int d=0; d<3; d++) {
	index += (iX[d] & 1) << (3 * l + d);
	iX[d] >>= 1;
      }
    }
    return index;
  }

  int getLevel(uint64_t key) {
    int level = -1;
    while( int(key) >= 0 ) {
      level++;
      key -= 1 << 3*level;
    }
    return level;
  }

  Box bounds2box(Bounds bounds) {
    vec3 Xmin = bounds.Xmin;
    vec3 Xmax = bounds.Xmax;
    Box box;
    for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2;
    box.R = 0;
    for (int d=0; d<3; d++) {
      box.R = std::max(box.X[d] - Xmin[d], box.R);
      box.R = std::max(Xmax[d] - box.X[d], box.R);
    }
    box.R *= 1.00001;
    return box;
  }

  void growTree(Bodies & bodies, int (* nodes)[10], int & numCells,
		int * permutation, int & numLevels, Box box) {
    logger::startTimer("Grow tree");
    const int maxLevel = 30;
    const int numBodies = bodies.size();
    int nbody8[8];
    int * iwork = new int [numBodies];
    int * levelOffset = new int [maxLevel];
    vec3 * Xj = new vec3 [numBodies];
    nodes[0][0] = 0;
    nodes[0][1] = 0;
    nodes[0][2] = 0;
    nodes[0][3] = 0;
    nodes[0][4] = 0;
    nodes[0][5] = 0;
    nodes[0][6] = 0;
    nodes[0][7] = 0;
    nodes[0][8] = numBodies;
    levelOffset[0] = 0;
    levelOffset[1] = 1;
    for (int i=0; i<numBodies; i++) {
      permutation[i] = i;
      Xj[i] = bodies[i].X;
    }
    numCells = 1;
    numLevels = 0;
    for (int level=0; level<maxLevel; level++) {
      for (int iparent=levelOffset[level]; iparent<levelOffset[level+1]; iparent++) {
	int nbody = nodes[iparent][8];
	if (nbody > ncrit) {
	  int ibody = nodes[iparent][7];
	  reorder(box, level, &nodes[iparent][1], Xj, &permutation[ibody], nbody, iwork, nbody8);
	  int nchild = 0;
	  int offset = ibody;
	  nodes[iparent][5] = numCells;
	  for (int i=0; i<8; i++) {
	    nodes[numCells][0] = level + 1;
	    nodes[numCells][1] = nodes[iparent][1] * 2 + i % 2;
	    nodes[numCells][2] = nodes[iparent][2] * 2 + (i / 2) % 2;
	    nodes[numCells][3] = nodes[iparent][3] * 2 + i / 4;
	    nodes[numCells][4] = iparent;
	    nodes[numCells][5] = 0;
	    nodes[numCells][6] = 0;
	    nodes[numCells][7] = offset;
	    nodes[numCells][8] = nbody8[i];
	    nchild++;
	    offset += nbody8[i];
	    numCells++;
	    numLevels=level+1;
	  }
	  nodes[iparent][6] = nchild;
	}
      }
      levelOffset[level+2] = numCells;
      if (levelOffset[level+1] == levelOffset[level+2]) break;
    }
    delete[] Xj;
    delete[] levelOffset;
    delete[] iwork;
    logger::stopTimer("Grow tree");
  }

  Cells linkTree(Bodies & bodies, Bodies & buffer, int (* nodes)[10], int numCells,
		 int * permutation, Box box) {
    logger::startTimer("Link tree");
    int numBodies = bodies.size();
    Cells cells(numCells);
    C_iter C = cells.begin();
    ivec3 iX;
    for (int i=0; i<numCells; i++,C++) {
      int level = nodes[i][0];
      iX[0]      = nodes[i][1];
      iX[1]      = nodes[i][2];
      iX[2]      = nodes[i][3];
      C->ICELL   = getKey(iX, level);
      C->IPARENT = nodes[i][4];
      C->ICHILD  = nodes[i][5];
      C->NCHILD  = nodes[i][6];
      C->IBODY   = nodes[i][7];
      C->NBODY   = nodes[i][8];
      real_t R = box.R / (1 << level);
      C->R = R;
      for (int d=0; d<3; d++) {
	C->X[d] = box.X[d] - box.R + iX[d] * R * 2 + R;
      }
    }
    buffer.resize(numBodies);
    for (int i=0; i<numBodies; i++) {
      buffer[i] = bodies[permutation[i]];
    }
    bodies = buffer;
    B_iter B = bodies.begin();
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {
      C->BODY = B + C->IBODY;
    }
    logger::stopTimer("Link tree");
    return cells;
  }

public:
  BuildTree(int _ncrit, int ) : ncrit(_ncrit) {}

  Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
    int numCells;
    int numBodies = bodies.size();
    int (* nodes)[10] = new int [numBodies][10]();
    int * permutation = new int [numBodies];
    Box box = bounds2box(bounds);
    growTree(bodies, nodes, numCells, permutation, numLevels, box);
    Cells cells = linkTree(bodies, buffer, nodes, numCells, permutation, box);
    delete[] permutation;
    delete[] nodes;
    return cells;
  }

  void printTreeData(Cells & cells) {
    if (logger::verbose && !cells.empty()) {                    // If verbose flag is true
      logger::printTitle("Tree stats");                         //  Print title
      std::cout  << std::setw(logger::stringLength) << std::left//  Set format
		 << "Bodies"     << " : " << cells.front().NBODY << std::endl// Print number of bodies
		 << std::setw(logger::stringLength) << std::left//  Set format
		 << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
		 << std::setw(logger::stringLength) << std::left//  Set format
		 << "Tree depth" << " : " << numLevels << std::endl;//  Print number of levels
    }                                                           // End if for verbose flag
  }
};

#endif
