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

int verify_tree(int **expansions, int** edges, 
		int** node_pointers, int *leaf_populations, 
		int node_id, int level){
  int charge = 0;
  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;
  if(isnotleaf==0){
    charge = leaf_populations[node_pointers[level][node_id]];
  }else{
    for(int i=children_start; i<children_stop; i++){
      charge += verify_tree(expansions, edges, node_pointers, 
			    leaf_populations, i, level+1); 
    }
  }
  expansions[level][node_id] = charge;
  return charge;
}

int verify_tree_wrapper(int **expansions, int **edges, 
			int **node_pointers, int *leaf_populations, 
			int nnodes, int N){
  int charge = 0;
  int pass = 0;
  for(int i=0; i<nnodes; i++){
    charge += verify_tree(expansions, edges, node_pointers, leaf_populations, i, 0); 
  }
  if(charge == N){
    pass = 1;
  }else{
    printf("charge: %d\n", charge);
  }
  return pass;
}

int verify_interactions_compressed(int **expansion, int **edges, uint32_t **nn_first, 
				   int **nn_list, uint32_t **fn_first, int **fn_list, 
				   uint32_t **common_first, int **common_list,
				   int node_id, int level, int parent_charge, 
				   int N, int tree_height){
  int pass = 1;
  int interactions = parent_charge;
  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;
  int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
  int nn_stop = nn_first[level][node_id];
  int numnn = nn_stop - nn_start;
  uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
  uint32_t fn_stop = fn_first[level][node_id];
  uint32_t numfn = fn_stop - fn_start;
  uint32_t common_start = (node_id==0) ? 0 : common_first[level][node_id-1];
  uint32_t common_stop = common_first[level][node_id];
  uint32_t numcomm = common_stop - common_start;
  if(numnn>0){
    for(int i=nn_start; i<nn_stop; i++){
      interactions += expansion[level][nn_list[level][i]];
    }
  }
  if(numfn>0){
    for(uint32_t i=fn_start; i<fn_stop; i++){
      interactions += expansion[level][(uint)fn_list[level][i]];
    }
  }
  if(numcomm>0){
    for(uint32_t i=common_start; i<common_stop; i++){
      interactions += expansion[level+1][(uint)common_list[level][i]];
    }
  }
  if(isnotleaf>0 && level<tree_height){
    for(int i=children_start; i<children_stop; i++){
      pass &= verify_interactions_compressed(expansion, edges, nn_first, 
					     nn_list, fn_first, fn_list,
					     common_first, common_list,
					     i, level+1, interactions, N, tree_height);
    }
  }else{
    if(interactions == N){
      pass = 1;
    }else{
      pass = 0;
    }
  }
  return pass;
}

int verify_interactions_compressed_wrapper(int **expansion, int **edges, 
					   uint32_t **nn_first, int **nn_list, 
					   uint32_t **fn_first, int **fn_list, 
					   uint32_t **common_first, int **common_list,
					   int nnodes, int N, int tree_height){
  int pass = 1;
  for(int i=0; i<nnodes; i++){
    pass &= verify_interactions_compressed(expansion, edges, nn_first, 
					   nn_list, fn_first, fn_list,
					   common_first, common_list,
					   i, 0, 0, N, tree_height);
  }
  return pass;
}
