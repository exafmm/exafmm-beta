#include <cilk/cilk.h>
#include <math.h>
#include <string.h>
#include "utils.h"

#define MAXBINS 8
#define MAXBINS64 64
#define NP 512
#define NP3 64
#define th 16

void relocate_data_radix6_noindex(uint32_t (*restrict pointIds), 
				  uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
				  int* str, int P, int M, int N, int sft){
#pragma ivdep
  for(int j=0; j<M; j++){
    if(P+j<N){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      int jj = str[ii];
      codes[jj] = zcodes[j];
      pointIds[jj] = P + j;
      jj++;
      str[ii]=jj;
    }
  }
}

void parallel_cpy(uint64_t (*restrict codes), uint64_t (*restrict zcodes), 
		  uint32_t (*restrict pointIds), uint32_t (*restrict index), int N){
  codes[0:N] = zcodes[0:N];
  pointIds[0:N] = index[0:N];
}

void find_bin_sizes(uint64_t (*restrict morton_codes), 
		    int (*restrict BinSizes), int sft, int M){
#pragma ivdep
  for(int i=0; i<M; i++){
    uint32_t ii = (morton_codes[i] >> sft) & 0x3F;
    BinSizes[ii]++;
  }
}

void bin_sort_serial_radix6_bitmap(uint64_t (*restrict zcodes), 
				   uint64_t (*restrict codes), 
				   uint32_t (*restrict pointIds), 
				   uint32_t (*restrict index), 
				   uint32_t (*restrict bit_map),
				   int N, 
				   int sft, int lv, int stop, 
				   int population_threshold){
  int BinSizes[MAXBINS64] = {0};
  int BinCursor[MAXBINS64] = {0};
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;
  bit_map[0] |= (1 << (2*lv-1));
  if(N<=population_threshold || sft < stop){
    pointIds[0:N] = index[0:N];
    codes[0:N] = zcodes[0:N]; 
    return;
  }else{
#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      BinSizes[ii]++;
    }
    int offset = 0;
#pragma ivdep
    for(int i=0; i<MAXBINS64; i++){
      int ss = BinSizes[i];
      BinCursor[i] = offset;
      offset += ss;
      BinSizes[i] = offset;
    }
#pragma ivdep
    for(int j=0; j<N; j++){
      uint32_t ii = (zcodes[j]>>sft) & 0x3F;
      pointIds[BinCursor[ii]] = index[j];
      codes[BinCursor[ii]] = zcodes[j];
      BinCursor[ii]++;
    }
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;
    int cursor = 0;
    int super_box_start = 0;
    offset = 0; 
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(i+1)*MAXBINS-1];
      if(super_box_size>0){
	bit_map[cursor] |= (1 << (2*lv));
      }
      cursor += super_box_size;
      if(super_box_size > th){
	for(int j=0; j<MAXBINS; j++){
	  int size = BinSizes[i*MAXBINS+j] - offset;
	  if( size>0 ){
	    cilk_spawn bin_sort_serial_radix6_bitmap(&zcodes[offset], 
						     &codes[offset], 
						     &pointIds[offset], 
						     &index[offset], 
						     &bit_map[offset], 
						     size, sft-6, 
						     lv+1, stop,
						     population_threshold);
	  }
	  offset += size;
	}
      }else{
	cilk_spawn parallel_cpy(&codes[offset], 
				&zcodes[offset], 
				&pointIds[offset], 
				&index[offset], super_box_size);	  
	offset += super_box_size;
      }
    }
    cilk_sync;    
  }
}

void bin_sort_radix6_bitmap(uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
			    uint32_t (*restrict pointIds), uint32_t (*restrict index),
			    uint32_t (*restrict bit_map),
			    int N, int sft, int lv, int stop, 
			    int population_threshold){

  int BinSizes[NP*MAXBINS64];
  int BinCursor[NP*MAXBINS64];
  uint32_t *tmp_ptr;
  uint64_t *tmp_code;  
  BinSizes[:] = 0;
  if(lv>0){
    bit_map[0] |= (1 << (2*lv-1));
  }
  if(N<=population_threshold || sft<stop){
    pointIds[0:N] = index[0:N];
    return;
  }else{
    int M = (int)ceil((float)N / (float)NP);
    for(int i=0; i<NP; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      cilk_spawn find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size); 
    }
    cilk_sync;
    int offset = 0;
    for(int i=0; i<MAXBINS64; i++){
#pragma ivdep
      for(int j=0; j<NP; j++){
	int size = BinSizes[j*MAXBINS64 + i];
	BinCursor[j*MAXBINS64 + i] = offset;
	offset += size; 
	BinSizes[j*MAXBINS64 + i] = offset;
      }
    }
    for(int i=0; i<NP; i++){
      cilk_spawn relocate_data_radix6_noindex(pointIds, 
					      &zcodes[i*M], codes, 
					      &BinCursor[i*MAXBINS64], i*M, 
					      M, N, sft);
    }
    cilk_sync;
    tmp_ptr = index;
    index = pointIds;
    pointIds = tmp_ptr;
    tmp_code = zcodes;
    zcodes = codes;
    codes = tmp_code;
    int cursor = 0;
    int super_box_start = 0;
    for(int i=0; i<MAXBINS; i++){
      int super_box_size = BinSizes[(NP-1)*MAXBINS64 + (i+1)*MAXBINS-1] - super_box_start;
      super_box_start = BinSizes[(NP-1)*MAXBINS64 + (i+1)*MAXBINS-1];
      for(int j=0; j<MAXBINS; j++){
	int start = (i==0 && j==0) ? 0 : BinSizes[(NP-1)*MAXBINS64 + i*MAXBINS + j-1];
	int size = BinSizes[(NP-1)*MAXBINS64 + i*MAXBINS + j] - start; 
	if( size > 0 ){
	  cilk_spawn bin_sort_serial_radix6_bitmap(&zcodes[start], 
						   &codes[start], 
						   &pointIds[start],
						   &index[start], 
						   &bit_map[start], 
						   size, sft-6, 
						   lv+1, stop,
						   population_threshold);
	}	
      }
      bit_map[cursor] |= 1;
      cursor += super_box_size;
    }
  }
  cilk_sync;
}

void build_tree(uint64_t (*restrict mcodes), 
		uint64_t (*restrict scodes), 
		uint32_t (*restrict permutation_vector), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold){
  bin_sort_radix6_bitmap(mcodes, scodes, permutation_vector, 
			 index, bit_map, 
			 N, 3*(maxlev-2), 
			 0, 3*(maxlev-maxheight), population_threshold);
}

void decomposeSpace(int N, uint64_t **mcodes,
		    uint32_t *permutation_vector, uint32_t *bit_map,
		    int maxlev, int population_threshold) {
  uint64_t *scodes = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Code buffer array");
  uint32_t *index = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), "Index vector");
  start_timer();
  build_tree(*mcodes, scodes, permutation_vector,
	     index, bit_map, N, maxlev, maxlev,
	     population_threshold);
  stop_timer("Tree building");
  uint64_t *tcodes = *mcodes;
  *mcodes = scodes;
  free(tcodes);
  free(index);
}

void relocateTL(float (*restrict Y), float (*restrict X), uint32_t (*restrict index), uint32_t N){
  for(uint i=0; i<N; i++){
    Y[i*LDIM:LDIM] = X[index[i]*LDIM:LDIM];
  }
}

void rearrange_dataTL(float (*restrict Y), float (*restrict X), uint32_t (*restrict Ids), int N){
  if(N>100000){
    uint64_t M = (uint32_t)ceil((float)N / (float)NP);
    for(uint64_t i=0; i<NP; i++){
      uint64_t size = ((i+1)*M < N) ? M : N - i*M;
      cilk_spawn relocateTL(&Y[i*M*LDIM], X, &Ids[i*M], size);
    }
    cilk_sync;
  }
  else{
    cilk_for(uint64_t i=0; i<N; i++){
      Y[i*LDIM:LDIM] = X[Ids[i]*LDIM:LDIM];
    }
  }
}

void rearrangeDirect(float (*restrict Y),
		     float (*restrict X),
		     uint32_t (*restrict Ids), int N){
  cilk_for(uint i=0; i<N; i++){
    memcpy(&Y[i*LDIM], &X[Ids[i]*LDIM], LDIM*sizeof(float));
  }
}

void relocateParticles(int N, float **X, uint32_t *permutation_vector){
  float *Y  = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle buffer");
  float *tmp;
  start_timer();
  rearrange_dataTL(Y, *X, permutation_vector, N);
  stop_timer("Relocate particles");
  tmp = Y;
  Y = *X;
  *X = tmp;
  free(Y);
}
