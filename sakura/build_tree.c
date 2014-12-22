#include <cilk/cilk.h>
#include <math.h> 
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
  if(N<population_threshold || sft < stop){
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
    int super_box_size = 0;
    offset = 0; 
    int i = 0;
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
  if(N<population_threshold || sft<stop){
    pointIds[0:N] = index[0:N];
    return;
  }else{
    int M = (int)ceil((float)N / (float)NP);
    for(int i=0; i<NP; i++){
      int size = ((i+1)*M < N) ? M : (N - i*M);
      cilk_spawn find_bin_sizes(&zcodes[i*M], &BinSizes[i*MAXBINS64], sft, size); 
    }
    cilk_sync;
    int dd = 0;
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

void count_bins_per_level_bitmap(int (*restrict nodes_per_level), 
				 uint32_t (*restrict bit_map), uint32_t (*restrict masks), 
				 int N, int L){
  nodes_per_level[0:L] = 0;
  for(int i=0; i<N; i++){
    for(int j=0; j<L; j++){
      if((bit_map[i] & masks[j]) > 0){
	nodes_per_level[j]++;
      }
    }
  }  
}

int count_bins_bitmap_wrapper(int (*restrict nodes_per_level), int (*restrict nodes_block_first),
			      uint32_t (*restrict bit_map), int N, int L){
  uint32_t *masks = (uint32_t *)malloc(L*sizeof(uint32_t));
  cilk_for(int i=0; i<L; i++){
    masks[i] = 1 << i;
  }
  int M = (int)ceil((float)N / (float)NP3); 
  int *nodes_per_level_buff = (int*)malloc(NP3*L*sizeof(int));
  for(int i=0; i<NP3; i++){
    int size = ((i+1)*M<N) ? M : N - i*M;
    cilk_spawn count_bins_per_level_bitmap(&nodes_per_level_buff[i*L], 
					   &bit_map[i*M], masks, 
					   size, L);
  }
  cilk_sync;
  nodes_per_level[0:L] = 0;
  nodes_block_first[0:L] = 0;
  for(int i=1; i<NP3; i++){
    nodes_block_first[i*L:L] = nodes_block_first[(i-1)*L:L] + 
      nodes_per_level_buff[(i-1)*L:L];
  }
  nodes_per_level[0:L] = nodes_block_first[(NP3-1)*L:L] + 
    nodes_per_level_buff[(NP3-1)*L:L];
  int height = L;
  for(int i=0; i<L; i++){
    if(nodes_per_level[i] == 0){
      height = i;
      break;
    }
  } 
  
  free(masks);
  free(nodes_per_level_buff);
  
  return height;
}

void parent_children_connection_singlepass(int (**restrict node_pointers),
					   int (**restrict num_children),
					   int (**restrict node_codes), 
					   int (*restrict child_counter),
					   int (*restrict remaining_child),
					   int (*restrict block_first), 
					   uint32_t (*restrict masks),
					   uint32_t (*restrict bit_map), 
					   uint64_t (*restrict leaf_morton_codes),
					   int N,
					   int L, int Base, int maxL, int maxlev){
  int cursor[20];
  cursor[0:L] = block_first[0:L];
  child_counter[0:L] = 0;
  remaining_child[0:L] = 0;  
  for(int i=0; i<N; i++){
    for(int j=0; j<L-1; j++){
      if( (bit_map[i] & masks[j]) > 0 ){
	node_pointers[j][cursor[j]] = Base + i;
	decode_morton_code(&node_codes[j][3*cursor[j]],
			   &node_codes[j][3*cursor[j]+1],
			   &node_codes[j][3*cursor[j]+2],
			   leaf_morton_codes[i] >> (3*(maxlev - j -1)) );
	if(cursor[j]>block_first[j]){
	  num_children[j][cursor[j]-1] = child_counter[j+1];
	}else{
	  remaining_child[j] = child_counter[j+1];
	}
	child_counter[j+1] = 0;
	cursor[j]++;
	child_counter[j]++;
      }
    }
    if( (bit_map[i] & masks[L-1]) > 0 ){       
      node_pointers[L-1][cursor[L-1]] = Base + i;
      decode_morton_code(&node_codes[L-1][3*cursor[L-1]],
			 &node_codes[L-1][3*cursor[L-1]+1],
			 &node_codes[L-1][3*cursor[L-1]+2],
			 leaf_morton_codes[i] >> (3*(maxlev - L)));
      cursor[L-1]++;
      child_counter[L-1]++;
    }
  }
  for(int j=0; j<L-1; j++){
    if(cursor[j]>block_first[j]){
      num_children[j][cursor[j]-1] = child_counter[j+1];
    }else{
      remaining_child[j] = child_counter[j+1];
    }
  }
}

void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					int (**restrict node_codes), 
					int (*restrict nodes_block_first),
					uint32_t (*restrict bit_map), 
					uint64_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev){
  int *child_counter = (int *)malloc(L*NP3*sizeof(int));
  int *remaining_child = (int *)malloc(L*NP3*sizeof(int));
  uint32_t *masks = (uint32_t *)malloc(L*sizeof(uint32_t));  
  int M = (int)ceil((float)N / (float)NP3);
  cilk_for(int i=0; i<L; i++){
    masks[i] = 1 << i;
  }
  for(int i=0; i<NP3; i++){  
    int size = ((i+1)*M < N) ? M : N - i*M;
    cilk_spawn parent_children_connection_singlepass(node_pointers,
						     num_children,
						     node_codes, 
						     &child_counter[i*L],
						     &remaining_child[i*L],
						     &nodes_block_first[i*maxL],
						     masks,
						     &bit_map[i*M], 
						     &leaf_morton_codes[i*M],
						     size,
						     L, i*M, maxL, maxlev);
  }
  cilk_sync;
  for(int i=0; i<NP3; i++){
    for(int j=0; j<L; j++){ 
      if(nodes_block_first[i*maxL+j]-1 >= 0){
	num_children[j][nodes_block_first[i*maxL+j]-1] += remaining_child[i*L + j];
      }
    }
  }
  free(child_counter);
  free(remaining_child);
  free(masks);
}

void first_child_position(int *children_first, int *num_children, int N){
  children_first[0] = num_children[0];
  for(int i=1; i<N; i++){ 
    children_first[i] = children_first[i-1] + num_children[i];
  }
}

void first_child_position_wrapper(int **children_first, 
				  int **num_children, 
				  int *nodes_per_level, 
				  int L){
  cilk_for(int i=0; i<L; i++){
    if(nodes_per_level[i]>0){
      first_child_position(children_first[i], num_children[i], nodes_per_level[i]);
    }
  }
}

void build_tree(float *Y, float *X, uint64_t (*restrict mcodes), 
		uint64_t (*restrict scodes), 
		uint32_t (*restrict permutation_vector), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist){
  bin_sort_radix6_bitmap(mcodes, scodes, permutation_vector, 
			 index, bit_map, 
			 N, 3*(maxlev-2), 
			 0, 3*(maxlev-maxheight), population_threshold);
}
