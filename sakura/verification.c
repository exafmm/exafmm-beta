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

int verify_interaction_iterative_singlearray(int **expansion, int *interactions,
					     int **c_count, uint32_t **n_count,
					     int **n_list, uint32_t **f_count,
					     int **f_list, uint32_t **s_count,
					     int **s_list,
					     int *nodes_per_level, int N, int height){
  int pass = 1;
  int level = 0;
  for(int glb_node_id=0; glb_node_id<nodes_per_level[height-1]; glb_node_id++){
    if(glb_node_id>=nodes_per_level[level]){
      level++;
    }
    int offset = (level==0) ? 0 : nodes_per_level[level-1];
    int node_id = glb_node_id - offset;
    int children_stop = c_count[level][node_id];
    int children_start = (node_id==0) ? 0 : c_count[level][node_id-1];
    int isnotleaf = children_stop - children_start;
    int nn_start = (node_id==0)? 0 : n_count[level][node_id-1];
    int nn_stop = n_count[level][node_id];
    int numnn = nn_stop - nn_start;
    uint32_t fn_start = (node_id==0) ? 0 : f_count[level][node_id-1];
    uint32_t fn_stop = f_count[level][node_id];
    uint32_t numfn = fn_stop - fn_start;
    uint32_t common_start = (node_id==0) ? 0 : s_count[level][node_id-1];
    uint32_t common_stop = s_count[level][node_id];
    uint32_t numcomm = common_stop - common_start;
    if(numnn>0){
      for(int i=nn_start; i<nn_stop; i++){
	interactions[glb_node_id] += expansion[level][n_list[level][i]];
      }
    }
    if(numfn>0){
      for(uint32_t i=fn_start; i<fn_stop; i++){
	interactions[glb_node_id] += expansion[level][(uint32_t)f_list[level][i]];
      }
    }
    if(numcomm>0){
      for(uint32_t i=common_start; i<common_stop; i++){
	interactions[glb_node_id] += expansion[level+1][(uint)s_list[level][i]];
      }
    }
    if(isnotleaf>0 && level<height){
      int nl_offset = nodes_per_level[level];
      for(int i=children_start; i<children_stop; i++){
	interactions[i+nl_offset] += interactions[glb_node_id];
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
  return(pass);
}


int verify_interactions_wrapper_iterative_singlearray(int **expansion, int **interactions,
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
  int *tmp_interactions = (int *)calloc(nodes_sum[height-1],sizeof(int));
  int pass = verify_interaction_iterative_singlearray(expansion, tmp_interactions,
						      c_count,
						      n_count, n_list,
						      f_count, f_list,
						      s_count, s_list,
						      nodes_sum, N, height);
  free(nodes_sum);
  free(tmp_interactions);
  return pass;
}
