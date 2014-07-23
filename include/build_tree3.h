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

#define NP 512
#define NP2 128
#define NCRIT 16
#define MAXBINS64 64

class BuildTree {
private:
  int maxlevel;

private:
  Box bounds2box(Bounds & bounds) {
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

  void compute_quantization_codes_T(uint* codes, float *X, int N, int nbins) {
    float Xmin[3] = {0};
    float Xmax[3] = {0};
    float X0[3];
    for (int b=0; b<N; b++) {
      for (int d=0; d<3; d++) {
	Xmin[d] = fmin(X[3*b+d],Xmin[d]);
	Xmax[d] = fmax(X[3*b+d],Xmax[d]);
      }
    }
    for (int d=0; d<3; d++) X0[d] = (Xmax[d] + Xmin[d]) / 2;
    float range = 0;
    for(int d=0; d<3; d++) {
      range = fmax(X0[d] - Xmin[d], range);
      range = fmax(Xmax[d] - X0[d], range);
    }
    range *= 1.00001;
    for(int d=0; d<3; d++) {
      Xmin[d] = X0[d] - range;
      Xmax[d] = X0[d] + range;
    }
    float d = range / nbins;
    cilk_for(int i=0; i<N; i++){
      codes[i*3:3] = floor((X[i*3:3] - Xmin[0:3]) / d);
    }
  }

  void relocate_data_radix6(uint* pointIds, uint* index, uint long* zcodes, uint long* codes, int* str, int P, int M, int N, int sft) {
#pragma ivdep
    for(int j=0; j<M; j++){
      if(P+j<N){
	uint ii = (zcodes[j]>>sft) & 0x3F;
	int jj = str[ii];
	codes[jj] = zcodes[j];
	pointIds[jj] = index[j];
	str[ii]=jj+1;
      }
    }
  }

  void bin_sort_serial_radix6(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv) {

    int BinSizes[MAXBINS64];
    int str[MAXBINS64];
    uint acm_sizes[MAXBINS64];
    uint* tmp_ptr;
    uint long* tmp_code;

    if(N<=NCRIT || sft<0){
      pointIds[0:N] = index[0:N];
      bins[0:N] = tid;                                  
      level[0:N] = lv-1;
      return;
    }

    BinSizes[:] = 0;
    str[:] = 0;
    acm_sizes[:] = 0;

#pragma ivdep
    for(int j=0; j<N; j++){
      uint ii = (zcodes[j]>>sft) & 0x3F;
      BinSizes[ii]++;
    }

    str[0] = 0;
    acm_sizes[0] = 0;
#pragma ivdep
    for(int i=1; i<MAXBINS64; i++){
      uint tmp = str[i-1] + BinSizes[i-1];
      str[i] = tmp;
      acm_sizes[i] = tmp;
    }

#pragma ivdep
    for(int j=0; j<N; j++){
      uint ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      pointIds[jj] = index[j];
      codes[jj] = zcodes[j];
      str[ii] = jj+1;
    }

    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    if (lv<2) {
      for(int i=0; i<MAXBINS64; i++){
	cilk_spawn bin_sort_serial_radix6(&zcodes[acm_sizes[i]], &codes[acm_sizes[i]], &pointIds[acm_sizes[i]], &index[acm_sizes[i]], &bins[acm_sizes[i]], &level[acm_sizes[i]], BinSizes[i], sft-6, 64*tid + i, lv+1);
      }
      cilk_sync;
    } else {
      for(int i=0; i<MAXBINS64; i++){
	bin_sort_serial_radix6(&zcodes[acm_sizes[i]], &codes[acm_sizes[i]], &pointIds[acm_sizes[i]], &index[acm_sizes[i]], &bins[acm_sizes[i]], &level[acm_sizes[i]], BinSizes[i], sft-6, 64*tid + i, lv+1);
      }
    }
  }

  void bin_sort_radix6(uint long *zcodes, uint long* codes, uint *pointIds, uint* index, uint long* bins, int *level, int N, int sft, int tid, int lv) {

    int BinSizes[NP*MAXBINS64];
    int str[NP*MAXBINS64];
    uint Sizes[MAXBINS64];
    uint acm_sizes[MAXBINS64];
    uint* tmp_ptr;
    uint long* tmp_code;  

    BinSizes[:] = 0;
    str[:] = 0;
    Sizes[:] = 0;
    acm_sizes[:] = 0;

    if(N<=NCRIT || sft<0){
      pointIds[0:N] = index[0:N];
      level[0] = lv-1;
      bins[0] = tid;
      return;
    }

    int M = (int)ceil((float)N / (float)NP);

    cilk_for(int i=0; i<NP; i++){
#pragma ivdep
      for(int j=0; j<M; j++){
	if(i*M+j<N){
	  uint ii = (zcodes[i*M + j]>>sft) & 0x3F;
	  BinSizes[i*MAXBINS64 + ii]++;
	}
      }
    }

    int dd = 0;
    for(int i=0; i<MAXBINS64; i++){
      str[i] = dd;
      acm_sizes[i] = dd;
#pragma ivdep
      for(int j=1; j<NP; j++){
	str[j*MAXBINS64+i] = str[(j-1)*MAXBINS64+i] + BinSizes[(j-1)*MAXBINS64+i];
	Sizes[i] += BinSizes[(j-1)*MAXBINS64+i];
      }
      dd = str[(NP-1)*MAXBINS64+i] + BinSizes[(NP-1)*MAXBINS64 + i];
      Sizes[i] += BinSizes[(NP-1)*MAXBINS64 + i];
    }
    
    for(int i=0; i<NP; i++){
      cilk_spawn relocate_data_radix6(pointIds, &index[i*M], &zcodes[i*M], codes, &str[i*MAXBINS64], i*M, M, N, sft);
    }
    cilk_sync;

    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;

    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;

    if (lv<2) {    
      for(int i=0; i<MAXBINS64; i++) {
	cilk_spawn bin_sort_serial_radix6(&zcodes[acm_sizes[i]], &codes[acm_sizes[i]], &pointIds[acm_sizes[i]], &index[acm_sizes[i]], &bins[acm_sizes[i]], &level[acm_sizes[i]], Sizes[i], sft-6, 64*tid + i, lv+1);
      }
    } else {
      for(int i=0; i<MAXBINS64; i++) {
	bin_sort_serial_radix6(&zcodes[acm_sizes[i]], &codes[acm_sizes[i]], &pointIds[acm_sizes[i]], &index[acm_sizes[i]], &bins[acm_sizes[i]], &level[acm_sizes[i]], Sizes[i], sft-6, 64*tid + i, lv+1);
      }
    }
    cilk_sync;    
  }

  void permuteBlock(Body * bodies, Bodies & buffer, uint * index, int N){
    for(int i=0; i<N; i++){
      bodies[i] = buffer[index[i]];
    }
  }

  void permute(Bodies & bodies, Bodies & buffer, uint * index, int N){
    int M = N / NP2;
    for(int i=0; i<NP2-1; i++){
      cilk_spawn permuteBlock(&bodies[i*M], buffer, &index[i*M], M);
    }
    permuteBlock(&bodies[(NP2-1)*M], buffer, &index[(NP2-1)*M], N-(NP2-1)*M);
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
    const int level = 6;
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
    int nbins = (1 << maxlev);
    uint64_t * key = new uint64_t [numBodies];
    uint64_t * buffer = new uint64_t [numBodies];
    uint * index2 = new uint [numBodies];
    int * permutation = new int [numBodies];
    Cells cells;
    for (int b=0; b<int(bodies.size()); b++) {
      index2[b] = b;
    }

    logger::startTimer("Morton key");
    compute_quantization_codes_T(codes, X, N, nbins);
    morton_encoding_T(mcodes, codes, N);
    b = 0;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, b++) {
      mcodes[b] /= 8;
      B->ICELL = mcodes[b];
    }
    logger::stopTimer("Morton key");

    logger::startTimer("Radix sort");
    bin_sort_radix6(mcodes, scodes, index2, index, bins, levels, N, 3*(maxlev-2), 0, 0);
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
