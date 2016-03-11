#ifndef build_tree_cilk_h
#define build_tree_cilk_h
#include <algorithm>
#include "logger.h"
#include "thread.h"
#include "types.h"

static const int blockSize 512
static const int ncrit 16
static const int nbins 64

Box bounds2box(Bounds & bounds) {
  vec3 Xmin = bounds.Xmin;
  vec3 Xmax = bounds.Xmax;
  Box box;
  for (int d=0; d<3; d++){
    box.X[d] = (Xmax[d] + Xmin[d]) / 2;
  }
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

void getKey(int numBodies, float * X, float * Xmin, float * Xmax, uint32_t * keys, int numLevels) {
  const int nbins = 1 << numLevels;
  Bounds bounds;
  for (int d=0; d<3; d++) {
    bounds.Xmin[d] = Xmin[d];
    bounds.Xmax[d] = Xmax[d];
  }
  Box box = bounds2box(bounds);
  for (int d=0; d<3; d++) {
    Xmin[d] = bounds.Xmin[d];
    Xmax[d] = bounds.Xmax[d];
  }
  float D = 2 * box.R / nbins;

  cilk_for (int b=0; b<numBodies; b++) {
    int ix = floor((X[12*b+0] - Xmin[0]) / D);
    int iy = floor((X[12*b+1] - Xmin[1]) / D);
    int iz = floor((X[12*b+2] - Xmin[2]) / D);
    keys[b] = morton::interleave(ix,iy,iz);
  }
}

void relocate(uint32_t * keys, uint32_t * buffer, uint32_t * index, uint32_t * permutation,
	      int * counter, int offset, int numBlock, int numBodies, int bitShift){
#pragma ivdep
  for (int i=0; i<numBlock; i++) {
    if (offset+i<numBodies) {
      int b = (keys[i] >> bitShift) & 0x3F;
      int c = counter[b];
      buffer[c] = keys[i];
      permutation[c] = index[i];
      counter[b]++;
    }
  }
}

void recursion(uint32_t * keys, uint32_t * buffer, uint32_t * permutation,
	       uint32_t * index, int numBodies, int bitShift) {
  //if (numBodies<=ncrit || bitShift<0) {
  if (bitShift<0) {
    permutation[0:numBodies] = index[0:numBodies];
    return;
  }

  int counter[nbins] = {0};
#pragma ivdep
  for (int i=0; i<numBodies; i++) {
    int b = (keys[i] >> bitShift) & 0x3F;
    counter[b]++;
  }

  int offset = 0;
#pragma ivdep
  for (int b=0; b<nbins; b++) {
    int size = counter[b];
    counter[b] = offset;
    offset += size;
  }

  relocate(keys, buffer, index, permutation, counter, 0, numBodies, numBodies, bitShift);
  std::swap(index, permutation);
  std::swap(keys, buffer);

  offset = 0;
  for (int b=0; b<nbins; b++) {
    int size = counter[b] - offset;
    recursion(&keys[offset], &buffer[offset], &permutation[offset], &index[offset],
	      size, bitShift-6);
    offset += size;
  }
}

void radixSort(int numBodies, uint32_t * keys, uint32_t * buffer,
	       uint32_t * permutation, uint32_t * index, int numLevels) {
  const int bitShift = 3 * (numLevels - 2);
  //if (numBodies<=ncrit || bitShift<0) {
  if (bitShift<0) {
    permutation[0:numBodies] = index[0:numBodies];
    return;
  }
  
  int numBlock = (numBodies - 1) / blockSize + 1;
  int counter[blockSize*nbins] = {0};
  cilk_for (int i=0; i<blockSize; i++) {
#pragma ivdep
    for (int j=0; j<numBlock; j++) {
      if (i*numBlock+j < numBodies) {
	int b = (keys[i*numBlock+j] >> bitShift) & 0x3F;
	counter[i*nbins+b]++;
      }
    }
  }
  
  int offset = 0;
  for (int b=0; b<nbins; b++) {
#pragma ivdep
    for (int i=0; i<blockSize; i++) {
      int size = counter[i*nbins+b];
      counter[i*nbins+b] = offset;
      offset += size;
    }
  }

  for (int i=0; i<blockSize; i++) {
    offset = i * numBlock;
    cilk_spawn relocate(&keys[offset], buffer, &index[offset], permutation,
			&counter[i*nbins], offset, numBlock, numBodies, bitShift);
  }
  cilk_sync;
  
  std::swap(index, permutation);
  std::swap(keys, buffer);
  
  offset = 0;
  for (int b=0; b<nbins; b++) {
    int size = counter[(blockSize-1)*nbins+b] - offset;
    cilk_spawn recursion(&keys[offset], &buffer[offset],
			 &permutation[offset], &index[offset], size, bitShift-6);
    offset += size;
  }
  cilk_sync;
}

void permuteBlock(float * buffer, float * bodies, uint32_t * index, int numBlock) {
#pragma ivdep
  for (int i=0; i<numBlock; i++) {
    for (int j=0; j<12; j++) {
      buffer[12*i+j] = bodies[12*index[i]+j];
    }
  }
}

void permute(int numBodies, float * bodies, float * buffer, uint32_t * index) {
  int numBlock = numBodies / blockSize;
  int offset = 0;
  for (int i=0; i<blockSize-1; i++) {
    cilk_spawn permuteBlock(&buffer[12*offset], bodies, &index[offset], numBlock);
    offset += numBlock;
  }
  permuteBlock(&buffer[12*offset], bodies, &index[offset], numBodies-offset);
}

class BuildTree {
private:
  const int ncrit;
  int numLevels;

private:
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
      for (int c=begin; c<end; c++) {
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

  void reverseOrder(Cells & cells, uint32_t * permutation) {
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
  BuildTree(int _ncrit, int) : numLevels(0), ncrit(_ncrit) {}

  Cells buildTree(Bodies & bodies, Bodies & buffer, Bounds bounds) {
    const int numBodies = bodies.size();
    const int level = numBodies >= ncrit ? 1 + int(log2(numBodies / ncrit)/3) : 0;
    numLevels = level;


    uint32_t * keys = new uint32_t [numBodies];
    uint32_t * keys_buffer = new uint32_t [numBodies];
    uint32_t * index = new uint32_t [numBodies];
    uint32_t * index2 = new uint32_t [numBodies];
    uint32_t * permutation = new uint32_t [numBodies];
    uint32_t * permutation2 = new uint32_t [numBodies];
    float * X = new float [12*numBodies];
    float * Y = new float [12*numBodies];

    int b=0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
      X[12*b+0] = B->X[0];
      X[12*b+1] = B->X[1];
      X[12*b+2] = B->X[2];
      X[12*b+3] = B->SRC;
      Y[12*b+0] = B->X[0];
      Y[12*b+1] = B->X[1];
      Y[12*b+2] = B->X[2];
      Y[12*b+3] = B->SRC;
      index[b] = b;
      index2[b] = b;
    }
    vec3 Xmin = bounds.Xmin;
    vec3 Xmax = bounds.Xmax;

    float * X2 = new float [12*numBodies];
    float Xmin2[3] = {0};
    Xmin2[0] = (float)Xmin[0]; Xmin2[1] = (float)Xmin[1]; Xmin2[2] = (float)Xmin[2];
    float Xmax2[3] = {0};
    Xmax2[0] = (float)Xmax[0]; Xmax2[1] = (float)Xmax[1]; Xmax2[2] = (float)Xmax[2];


    uint32_t * keys2 = new uint32_t [numBodies];
    for (int i=0; i<numBodies; i++) {
      for (int j=0; j<12; j++) {
	X2[12*i+j] = X[12*i+j];
      }
      keys2[i] = keys[i];
    }
    logger::startTimer("Grow tree");
    logger::startTimer("Morton key");
    getKey(numBodies, X, Xmin2, Xmax2, keys, numLevels);
    logger::stopTimer("Morton key");
    printf("Key comparison:\n");

    
    delete[] X2;
    delete[] keys2;

    logger::startTimer("Radix sort");
    radixSort(numBodies, keys, keys_buffer, permutation, index, numLevels);
    logger::stopTimer("Radix sort");
    logger::stopTimer("Grow tree",0);

    logger::startTimer("Grow tree");
    logger::startTimer("Permutation");
    permute(numBodies, Y, X, permutation);
    logger::stopTimer("Permutation");
    logger::stopTimer("Grow tree",0);

    b=0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
      B->X[0] = X[12*b+0];
      B->X[1] = X[12*b+1];
      B->X[2] = X[12*b+2];
      B->SRC  = X[12*b+3];
      B->IBODY = 0;
      B->IRANK = 0;
      B->ICELL = keys_buffer[b];
      B->WEIGHT = 0;
      B->TRG = 0;
    }

    Cells cells;
#if 1
    logger::startTimer("Link tree");
    logger::startTimer("Bodies to leafs");
    bodies2leafs(bodies, cells, bounds, level);
    logger::stopTimer("Bodies to leafs");

    logger::startTimer("Leafs to cells");
    leafs2cells(cells, bounds, level);
    logger::stopTimer("Leafs to cells");

    logger::startTimer("Reverse order");
    reverseOrder(cells, permutation);
    logger::stopTimer("Reverse order");
    logger::stopTimer("Link tree",0);
#endif

    delete[] keys;
    delete[] keys_buffer;
    delete[] index;
    delete[] permutation;
    delete[] X;
    delete[] Y;
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
		 << "Tree depth" << " : " << numLevels << std::endl;//  Print number of levels
    }                                                           // End if for verbose flag
  }
};
#endif
