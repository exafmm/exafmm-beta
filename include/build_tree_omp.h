#ifndef build_tree_omp_h
#define build_tree_omp_h
#include <algorithm>
#include "logger.h"
#include "thread.h"
#include "types.h"

class BuildTree {
private:
  const int ncrit;
  int numLevels;

private:
  //! Transform Xmin & Xmax to X (center) & R (radius)
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

  //! Calculate the Morton key
  inline void getKey(Bodies &bodies, uint64_t * key, Bounds bounds, int level) {
    Box box = bounds2box(bounds);
    float d = 2 * box.R / (1 << level);
#pragma omp parallel for
    for (int b=0; b<int(bodies.size()); b++) {
      B_iter B=bodies.begin()+b;
      int ix = (B->X[0] - bounds.Xmin[0]) / d;
      int iy = (B->X[1] - bounds.Xmin[1]) / d;
      int iz = (B->X[2] - bounds.Xmin[2]) / d;
      int id = 0;
      for( int l=0; l!=level; ++l ) {
	id += (ix & 1) << (3 * l);
	id += (iy & 1) << (3 * l + 1);
	id += (iz & 1) << (3 * l + 2);
	ix >>= 1;
	iy >>= 1;
	iz >>= 1;
      }
      key[b] = id;
      B->ICELL = id;
    }
  }

  void radixSort(uint64_t * key, int * value, uint64_t * buffer, int * permutation, int size) {
    const int bitStride = 8;
    const int stride = 1 << bitStride;
    const int mask = stride - 1;
    int bucket[stride];
    uint64_t maxKey = 0;
    for (int i=0; i<size; i++)
      if (key[i] > maxKey)
	maxKey = key[i];
    while (maxKey > 0) {
      for (int i=0; i<stride; i++)
	bucket[i] = 0;
      for (int i=0; i<size; i++)
	bucket[key[i] & mask]++;
      for (int i=1; i<stride; i++)
	bucket[i] += bucket[i-1];
      for (int i=size-1; i>=0; i--)
	permutation[i] = --bucket[key[i] & mask];
      for (int i=0; i<size; i++)
	buffer[permutation[i]] = value[i];
      for (int i=0; i<size; i++)
	value[i] = buffer[i];
      for (int i=0; i<size; i++)
	buffer[permutation[i]] = key[i];
      for (int i=0; i<size; i++)
	key[i] = buffer[i] >> bitStride;
      maxKey >>= bitStride;
    }
  }

  void permute(Bodies & bodies, Bodies & buffer, int * index) {
    const int n = bodies.size();
    for (int b=0; b<n; b++)
      bodies[b] = buffer[index[b]];
  }

  void bodies2leafs(Bodies & bodies, Cells & cells, Bounds bounds, int level) {
    int I = -1;
    C_iter C;
    cells.reserve(1 << (3 * level));
    Box box = bounds2box(bounds);
    float d = 2 * box.R / (1 << level);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      int IC = B->ICELL;
      int ix = (B->X[0] - bounds.Xmin[0]) / d;
      int iy = (B->X[1] - bounds.Xmin[1]) / d;
      int iz = (B->X[2] - bounds.Xmin[2]) / d;
      if( IC != I ) {
	Cell cell;
	cell.NCHILD = 0;
	cell.NBODY  = 0;
	cell.ICHILD = 0;
	cell.BODY   = B;
	cell.ICELL  = IC;
	cell.X[0]   = d * (ix + .5) + bounds.Xmin[0];
	cell.X[1]   = d * (iy + .5) + bounds.Xmin[1];
	cell.X[2]   = d * (iz + .5) + bounds.Xmin[2];
	cell.R      = d * .5;
	cells.push_back(cell);
	C = cells.end()-1;
	I = IC;
      }
      C->NBODY++;
    }
  }

  void leafs2cells(Cells & cells, Bounds bounds, int level) {
    int begin = 0, end = cells.size();
    Box box = bounds2box(bounds);
    float d = 2 * box.R / (1 << level);
    for (int l=1; l<=level; l++) {
      int div = (1 << (3 * l));
      d *= 2;
      int I = -1;
      int p = end - 1;
      for (int c=begin; c!=end; c++) {
	B_iter B = cells[c].BODY;
	int IC = B->ICELL / div;
	int ix = (B->X[0] - bounds.Xmin[0]) / d;
	int iy = (B->X[1] - bounds.Xmin[1]) / d;
	int iz = (B->X[2] - bounds.Xmin[2]) / d;
	if (IC != I) {
	  Cell cell;
	  cell.NCHILD = 0;
	  cell.NBODY  = 0;
	  cell.ICHILD = c;
	  cell.BODY   = cells[c].BODY;
	  cell.ICELL  = IC;
	  cell.X[0]   = d * (ix + .5) + bounds.Xmin[0];
	  cell.X[1]   = d * (iy + .5) + bounds.Xmin[1];
	  cell.X[2]   = d * (iz + .5) + bounds.Xmin[2];
	  cell.R      = d * .5;
	  cells.push_back(cell);
	  p++;
	  I = IC;
	}
	cells[p].NCHILD++;
	cells[p].NBODY += cells[c].NBODY;
	cells[c].IPARENT = p;
      }
      begin = end;
      end = cells.size();
    }
    cells.back().IPARENT = 0;
  }

  void reverseOrder(Cells & cells, int * permutation) {
    const int numCells = cells.size();
    int ic = numCells - 1;
    for (int c=0; c<numCells; c++,ic--) {
      permutation[c] = ic;
    }
    for (C_iter C=cells.begin(); C!=cells.end(); C++) {
      C->ICHILD = permutation[C->ICHILD] - C->NCHILD + 1;
      C->IPARENT = permutation[C->IPARENT];
    }
    std::reverse(cells.begin(), cells.end());
  }

public:
  BuildTree(int _ncrit, int) : ncrit(_ncrit), numLevels(0) {}

  Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
    const int numBodies = bodies.size();
    const int level = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
    numLevels = level;
    uint64_t * key = new uint64_t [numBodies];
    uint64_t * key_buffer = new uint64_t [numBodies];
    int * index = new int [numBodies];
    int * permutation = new int [numBodies];
    Cells cells;
    for (int b=0; b<int(bodies.size()); b++) {
      index[b] = b;
    }

    logger::startTimer("Morton key");
    getKey(bodies, key, bounds, level);
    logger::stopTimer("Morton key");

    logger::startTimer("Radix sort");
    radixSort(key, index, key_buffer, permutation, numBodies);
    logger::stopTimer("Radix sort");

    logger::startTimer("Copy buffer");
    buffer = bodies;
    logger::stopTimer("Copy buffer");

    logger::startTimer("Permutation");
    permute(bodies, buffer, index);
    logger::stopTimer("Permutation");

    logger::startTimer("Bodies to leafs");
    bodies2leafs(bodies, cells, bounds, level);
    logger::stopTimer("Bodies to leafs");

    logger::startTimer("Leafs to cells");
    leafs2cells(cells, bounds, level);
    logger::stopTimer("Leafs to cells");

    logger::startTimer("Reverse order");
    reverseOrder(cells, permutation);
    logger::stopTimer("Reverse order");

    delete[] key;
    delete[] key_buffer;
    delete[] index;
    delete[] permutation;
    return cells;
  }

  //! Print tree structure statistics
  void printTreeData(Cells & cells) {
    if (logger::verbose && !cells.empty()) {
      logger::printTitle("Tree stats");
      std::cout  << std::setw(logger::stringLength) << std::left
		 << "Bodies"     << " : " << cells.front().NBODY << std::endl
		 << std::setw(logger::stringLength) << std::left
		 << "Cells"      << " : " << cells.size() << std::endl
		 << std::setw(logger::stringLength) << std::left
		 << "Tree depth" << " : " << numLevels << std::endl;
    }
  }
};
#endif
