#include "utils.h"

int verify_tree(int **expansions, int** c_count2,
		int** node_pointers2, int *leaf_populations,
		int node_id, int level){
  int charge = 0;
  int children_stop = c_count2[level][node_id];
  int children_start = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int isnotleaf = children_stop - children_start;
  if(isnotleaf==0){
    charge = leaf_populations[node_pointers2[level][node_id]];
  }else{
    for(int i=children_start; i<children_stop; i++){
      charge += verify_tree(expansions, c_count2, node_pointers2,
			    leaf_populations, i, level+1);
    }
  }
  expansions[level][node_id] = charge;
  return charge;
}
