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
