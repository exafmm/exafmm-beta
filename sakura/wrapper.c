#include "utils.h"

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
