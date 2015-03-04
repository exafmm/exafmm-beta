#include "cilk/cilk.h"
#include <math.h>
#include "utils.h"

inline uint64_t Compact1By2(uint64_t x){
  x &= 0x1249249249249249;
  x = (x ^ (x >>  2)) & 0x10c30c30c30c30c3;
  x = (x ^ (x >>  4)) & 0x100f00f00f00f00f;
  x = (x ^ (x >>  8)) & 0x1f0000ff0000ff;
  x = (x ^ (x >> 16)) & 0x1f00000000ffff;
  x = (x ^ (x >> 32)) & 0x1fffff;
  return x;
}

void decode_morton_code(int *x, int *y, int *z, uint64_t mcode){
  x[0] = Compact1By2(mcode >> 0);
  y[0] = Compact1By2(mcode >> 1);
  z[0] = Compact1By2(mcode >> 2);
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
  int nthreads = (N > 10000) ? NP3 : NP3/4;
  int PP3 = (N > NP3) ? nthreads : 1;
  int M = (int)ceil((float)N / (float)PP3);
  int *nodes_per_level_buff = (int*)malloc(PP3*L*sizeof(int));
  for(int i=0; i<PP3; i++){
    int size = ((i+1)*M<N) ? M : N - i*M;
    cilk_spawn count_bins_per_level_bitmap(&nodes_per_level_buff[i*L],
					   &bit_map[i*M], masks,
					   size, L);
  }
  cilk_sync;
  nodes_per_level[0:L] = 0;
  nodes_block_first[0:L] = 0;
  for(int i=1; i<PP3; i++){
    nodes_block_first[i*L:L] = nodes_block_first[(i-1)*L:L] +
      nodes_per_level_buff[(i-1)*L:L];
  }
  nodes_per_level[0:L] = nodes_block_first[(PP3-1)*L:L] +
    nodes_per_level_buff[(PP3-1)*L:L];
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
					   int L, int Base, int maxlev){
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
			   leaf_morton_codes[i] >> (3*(maxlev - j -1)));
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
  int nthreads = (N > 10000) ? NP3 : NP3/4;
  int PP3 = (N > NP3) ? nthreads : 1;
  int M = (int)ceil((float)N / (float)PP3);
  cilk_for(int i=0; i<L; i++){
    masks[i] = 1 << i;
  }
  for(int i=0; i<PP3; i++){  
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
						     L, i*M, maxlev);
  }
  cilk_sync;
  for(int i=0; i<PP3; i++){
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

int tree_formation(uint32_t *bit_map, uint64_t *scodes,
		   int *nodes_per_level, int **node_pointers,
		   int **num_children, int **children_first,
		   int **node_codes, int maxlevel, int N){
  int *nodes_block_first = (int*)sakura_malloc(NP3*maxlevel, sizeof(int), "Node block");
  start_timer();
  int height = count_bins_bitmap_wrapper(nodes_per_level, nodes_block_first, bit_map, N, maxlevel);
  stop_timer("Count nodes");
  for(int i=0; i<height; i++){
    node_pointers[i] = (int *)sakura_malloc(nodes_per_level[i], sizeof(int), "Array of pointers to data");
    num_children[i] = (int *)sakura_malloc(nodes_per_level[i],sizeof(int), "Number of children array");
    children_first[i] = (int *)sakura_calloc(nodes_per_level[i], sizeof(int), "Children first array");
    node_codes[i] = (int *)sakura_malloc(3*nodes_per_level[i], sizeof(int), "Node hash codes");
  }
  start_timer();
  parent_children_connection_wrapper(node_pointers,
				     num_children,
				     node_codes,
				     nodes_block_first,
				     bit_map, scodes,
				     N, height, maxlevel, maxlevel);
  num_children[height-1][0:nodes_per_level[height-1]] = 0;
  stop_timer("Link children");
  start_timer();
  first_child_position_wrapper(children_first,
			       num_children,
			       nodes_per_level,
			       height);
  stop_timer("Find first child");
  free(nodes_block_first);
  return(height);
}

uint64_t find_leaf_populations(int *leaf_populations, uint32_t *bit_map2, int N){
  int pointer = 0;
  uint64_t numleaves = 0;
  for (int i=1; i<N; i++){
    if(bit_map2[i]>0){
      int pop = i - pointer;
      leaf_populations[pointer] = pop;
      pointer = i;
      numleaves++;
    }
  }
  leaf_populations[pointer] = N - pointer;
  return(numleaves);
}
