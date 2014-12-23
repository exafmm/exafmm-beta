#include "utils.h"

void form_interaction_lists(int **node_codes, int **children_first,
			    int **node_codes2, int **children_first2, 
			    uint32_t (**restrict nn_count), 
			    uint32_t (**restrict clgs_count), 
			    uint32_t (**restrict common_count),
			    int (**restrict nn_link_list), 
			    int (**restrict clgs_link_list),
			    int (**restrict common_list),
			    int (*restrict nodes_per_level), 
			    int (*restrict nodes_per_level2), 
			    int height){
  double tmp_list_physical = 0;
  double tmp_list_workspace = 0;
  for(int i=0; i<height; i++){
    nn_count[i] = (uint32_t *)sakura_calloc(nodes_per_level[i], sizeof(uint32_t), 
					    "Counters for near neighbors");
    clgs_count[i] = (uint32_t *)sakura_calloc(nodes_per_level[i], sizeof(uint32_t), 
					      "Counters for far neighbors");
    common_count[i] = (uint32_t *)sakura_calloc(nodes_per_level[i], sizeof(uint32_t), 
						"Counters for the common neighbors");
  }
  interaction_list_formation(node_codes, children_first, node_codes2,
			     children_first2, nn_count,
			     clgs_count, common_count, 
			     nn_link_list, clgs_link_list, 
			     common_list, nodes_per_level, nodes_per_level2, 
			     height, &tmp_list_physical, &tmp_list_workspace);
}

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
