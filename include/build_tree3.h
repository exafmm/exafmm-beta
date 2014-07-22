#ifndef build_tree_h
#define build_tree_h
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h> // TODO: Is this needed?
#include "logger.h"
#include "thread.h"
#include "types.h"
#include "utils.h"
#ifndef _OPENMP
int omp_get_num_threads() {return 1;}
int omp_get_thread_num() {return 0;}
#else
#include <omp.h>
#endif

#define NP 128

class BuildTree {
private:
  int maxlevel;

private:
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
    bounds.Xmin = box.X - box.R;
    bounds.Xmax = box.X + box.R;
    return box;
  }

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
    }                                                           // End loop over bodies
  }

  void permuteBlock(Body * bodies, Bodies & buffer, uint * index, int N){
    for(int i=0; i<N; i++){
      bodies[i] = buffer[index[i]];
    }
  }

  void permute(Bodies & bodies, Bodies & buffer, uint * index, int N){
    int M = N / NP;
    for(int i=0; i<NP-1; i++){
      cilk_spawn permuteBlock(&bodies[i*M], buffer, &index[i*M], M);
    }
    permuteBlock(&bodies[(NP-1)*M], buffer, &index[(NP-1)*M], N-(NP-1)*M);
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
    const int N = numBodies; // TODO: Change to numBodies
    const int level = 5;
    maxlevel = level;

    float *X = (float*)malloc(3*N*sizeof(float));
    uint* codes = (uint*)malloc(3*N*sizeof(uint)); // TODO: Use new
    unsigned long int* mcodes = (unsigned long int*)malloc(N*sizeof(unsigned long int));
    uint long* scodes = (uint long*)malloc(N*sizeof(uint long)); // TODO: Use standard types
    uint* pointIds = (uint*)malloc(N*sizeof(uint));
    uint* index = (uint*)malloc(N*sizeof(uint));
    uint long* bins = (uint long*)malloc(N*sizeof(uint long*));
    int* levels = (int*)malloc(N*sizeof(int));
    int b = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
      X[3*b+0] = B->X[0];
      X[3*b+1] = B->X[1];
      X[3*b+2] = B->X[2];
      index[b] = b;
    }
    int maxlev = 6;
    int maxheight = 6;
    int nbins = (1 << maxlev);
    compute_quantization_codes_T(codes, X, N, nbins);
    morton_encoding_T(mcodes, codes, N);
    uint64_t * key = new uint64_t [numBodies];
    uint64_t * buffer = new uint64_t [numBodies];
    uint * index2 = new uint [numBodies];
    int * permutation = new int [numBodies];
    Cells cells;
    for (int b=0; b<int(bodies.size()); b++) {
      index2[b] = b;
    }

    logger::startTimer("Morton key");
    getKey(bodies, key, bounds, level);
    logger::stopTimer("Morton key");

    logger::startTimer("Radix sort");
    for (int b=0; b<int(bodies.size()); b++) {
      mcodes[b] = key[b];
    }
    bin_sort_radix6(mcodes, scodes, index2, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));
    logger::stopTimer("Radix sort");

    Bodies bodies2 = bodies;
    logger::startTimer("Permutation");
    permute(bodies, bodies2, index2, N);
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
    delete[] index2;
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
