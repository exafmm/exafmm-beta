#include "utils.h"

uint64_t find_leaf_populations(int *populations, uint32_t *bit_map, int N){
  int pointer = 0;
  uint64_t numleaves = 0;  
  for (int i=1; i<N; i++){
    if(bit_map[i]>0){
      int pop = i - pointer;
      populations[pointer] = pop;
      pointer = i;
      numleaves++;
    }
  }
  populations[pointer] = N - pointer;
  return(numleaves);
}

int verify_tree(int **expansions, int** c_count, 
		int** node_pointers, int *leaf_populations, 
		int node_id, int level){
  int charge = 0;
  int children_stop = c_count[level][node_id];
  int children_start = (node_id==0) ? 0 : c_count[level][node_id-1];
  int isnotleaf = children_stop - children_start;
  if(isnotleaf==0){
    charge = leaf_populations[node_pointers[level][node_id]];
  }else{
    for(int i=children_start; i<children_stop; i++){
      charge += verify_tree(expansions, c_count, node_pointers, 
			    leaf_populations, i, level+1); 
    }
  }
  expansions[level][node_id] = charge;
  return charge;
}

int verify_tree_wrapper(int **expansions, int **c_count, 
			int **node_pointers, int *leaf_populations, 
			int nnodes, int N){
  int charge = 0;
  int pass = 0;
  for(int i=0; i<nnodes; i++){
    charge += verify_tree(expansions, c_count, node_pointers, leaf_populations, i, 0); 
  }
  if(charge == N){
    pass = 1;
  }else{
    printf("charge: %d\n", charge);
  }
  return pass;
}

int verify_interactions_wrapper_iterative_singlearray(int **expansions,
						      int **c_count,
						      uint32_t **n_count, int **n_list,
						      uint32_t **f_count, int **f_list,
						      uint32_t **s_count, int  **s_list,
						      int *nodes_per_level, int N,
						      int height){
  int *nodes_sum = (int *)malloc(height*sizeof(int));
  nodes_sum[0] = nodes_per_level[0];
  for(int i=1; i<height; i++){
    nodes_sum[i] = nodes_sum[i-1] + nodes_per_level[i];
  }
  int *interactions = (int *)calloc(nodes_sum[height-1],sizeof(int));
  int pass = 1;
  int level = 0;
  for(int glb_node_id=0; glb_node_id<nodes_sum[height-1]; glb_node_id++){
    if(glb_node_id>=nodes_sum[level]){
      level++;
    }
    int offset = (level==0) ? 0 : nodes_sum[level-1];
    int node_id = glb_node_id - offset;
    int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
    int c_end = c_count[level][node_id];
    int n_begin = (node_id==0)? 0 : n_count[level][node_id-1];
    int n_end = n_count[level][node_id];
    int f_begin = (node_id==0) ? 0 : f_count[level][node_id-1];
    int f_end = f_count[level][node_id];
    int s_begin = (node_id==0) ? 0 : s_count[level][node_id-1];
    int s_end = s_count[level][node_id];
    for(int i=n_begin; i<n_end; i++){
      interactions[glb_node_id] += expansions[level][n_list[level][i]];
    }
    for(int i=f_begin; i<f_end; i++){
      interactions[glb_node_id] += expansions[level][f_list[level][i]];
    }
    for(int i=s_begin; i<s_end; i++){
      interactions[glb_node_id] += expansions[level+1][s_list[level][i]];
    }
    if(level<height){
      int offset = nodes_sum[level];
      for(int i=c_begin; i<c_end; i++){
	interactions[i+offset] += interactions[glb_node_id];
      }
    }
    else{
      if(interactions[glb_node_id] == N){
	pass &= 1;
      }
      else{
	pass &= 0;
      }
    }
  }
  free(nodes_sum);
  free(interactions);
  return pass;
}
