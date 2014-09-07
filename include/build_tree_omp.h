#ifndef build_tree_h
#define build_tree_h
#include <algorithm>
#include "logger.h"
#include "thread.h"
#include "types.h"
#ifndef _OPENMP
int omp_get_num_threads() {return 1;}
int omp_get_thread_num() {return 0;}
#else
#include <omp.h>
#endif

class BuildTree {
private:
  int maxlevel;

private:
  //! Transform Xmin & Xmax to X (center) & R (radius)
  Box bounds2box(Bounds bounds) {
    vec3 Xmin = bounds.Xmin;                                    // Set local Xmin
    vec3 Xmax = bounds.Xmax;                                    // Set local Xmax
    Box box;                                                    // Bounding box
    for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2; // Calculate center of domain
    box.R = 0;                                                  // Initialize localRadius
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      box.R = std::max(box.X[d] - Xmin[d], box.R);              //  Calculate min distance from center
      box.R = std::max(Xmax[d] - box.X[d], box.R);              //  Calculate max distance from center
    }                                                           // End loop over dimensions
    box.R *= 1.00001;                                           // Add some leeway to radius
    return box;                                                 // Return box.X and box.R
  }

  //! Calculate the Morton key
  inline void getKey(Bodies &bodies, uint64_t * key, Bounds bounds, int level) {
    Box box = bounds2box(bounds);
    float d = 2 * box.R / (1 << level);                         // Cell size at current level
#pragma omp parallel for
    for (int b=0; b<int(bodies.size()); b++) {                  // Loop over bodies
      B_iter B=bodies.begin()+b;                                //  Body iterator
      int ix = (B->X[0] - bounds.Xmin[0]) / d;                  //  Index in x dimension
      int iy = (B->X[1] - bounds.Xmin[1]) / d;                  //  Index in y dimension
      int iz = (B->X[2] - bounds.Xmin[2]) / d;                  //  Index in z dimension
      int id = 0;                                               //  Initialize Morton key
      for( int l=0; l!=level; ++l ) {                           //  Loop over levels
	id += (ix & 1) << (3 * l);                              //   Interleave x bit
	id += (iy & 1) << (3 * l + 1);                          //   Interleave y bit
	id += (iz & 1) << (3 * l + 2);                          //   Interleave z bit
	ix >>= 1;                                               //   Shift x index
	iy >>= 1;                                               //   Shift y index
	iz >>= 1;                                               //   Shift z index
      }                                                         //  End loop over levels
      key[b] = id;                                              //  Store Morton key in array
      B->ICELL = id;                                            //  Store Morton key in body struct
    }                                                           // End loop over bodies
  }

  void radixSort(uint64_t * key, int * value, uint64_t * buffer, int * permutation, int size) {
    const int bitStride = 8;
    const int stride = 1 << bitStride;
    const int mask = stride - 1;
    int numThreads;
    int (*bucketPerThread)[stride];
    uint64_t maxKey = 0;
    uint64_t * maxKeyPerThread;
#pragma omp parallel
    {
      numThreads = omp_get_num_threads();
#pragma omp single
      {
	bucketPerThread = new int [numThreads][stride]();
	maxKeyPerThread = new uint64_t [numThreads];
	for (int i=0; i<numThreads; i++)
	  maxKeyPerThread[i] = 0;
      }
#pragma omp for
      for (int i=0; i<size; i++)
	if (key[i] > maxKeyPerThread[omp_get_thread_num()])
	  maxKeyPerThread[omp_get_thread_num()] = key[i];
#pragma omp single
      for (int i=0; i<numThreads; i++)
	if (maxKeyPerThread[i] > maxKey) maxKey = maxKeyPerThread[i];
      while (maxKey > 0) {
	int bucket[stride] = {0};
#pragma omp single
	for (int t=0; t<numThreads; t++)
	  for (int i=0; i<stride; i++)
	    bucketPerThread[t][i] = 0;
#pragma omp for
	for (int i=0; i<size; i++)
	  bucketPerThread[omp_get_thread_num()][key[i] & mask]++;
#pragma omp single
	{
	  for (int t=0; t<numThreads; t++)
	    for (int i=0; i<stride; i++)
	      bucket[i] += bucketPerThread[t][i];
	  for (int i=1; i<stride; i++)
	    bucket[i] += bucket[i-1];
	  for (int i=size-1; i>=0; i--)
	    permutation[i] = --bucket[key[i] & mask];
	}
#pragma omp for
	for (int i=0; i<size; i++)
	  buffer[permutation[i]] = value[i];
#pragma omp for
	for (int i=0; i<size; i++)
	  value[i] = buffer[i];
#pragma omp for
	for (int i=0; i<size; i++)
	  buffer[permutation[i]] = key[i];
#pragma omp for
	for (int i=0; i<size; i++)
	  key[i] = buffer[i] >> bitStride;
#pragma omp single
	maxKey >>= bitStride;
      }
    }
    delete[] bucketPerThread;
    delete[] maxKeyPerThread;
  }

  void permute(Bodies & bodies, int * index) {
    const int n = bodies.size();
    Bodies buffer = bodies;
#pragma omp parallel for
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
  BuildTree(int, int) : maxlevel(0) {}

  Cells buildTree(Bodies & bodies, Bounds bounds) {
    const int numBodies = bodies.size();
    const int level = 6;
    maxlevel = level;
    uint64_t * key = new uint64_t [numBodies];
    uint64_t * buffer = new uint64_t [numBodies];
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
    radixSort(key, index, buffer, permutation, numBodies);
    logger::stopTimer("Radix sort");

    logger::startTimer("Permutation");
    permute(bodies, index);
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
    delete[] buffer;
    delete[] index;
    delete[] permutation;
    return cells;
  }

  //! Print tree structure statistics
  void printTreeData(Cells & cells) {
    if (logger::verbose && !cells.empty()) {                    // If verbose flag is true
      logger::printTitle("Tree stats");                         //  Print title
      std::cout  << std::setw(logger::stringLength) << std::left//  Set format
		 << "Bodies"     << " : " << cells.front().NBODY << std::endl// Print number of bodies
		 << std::setw(logger::stringLength) << std::left//  Set format
		 << "Cells"      << " : " << cells.size() << std::endl// Print number of cells
		 << std::setw(logger::stringLength) << std::left//  Set format
		 << "Tree depth" << " : " << maxlevel << std::endl;//  Print number of levels
    }                                                           // End if for verbose flag
  }
};
#endif
