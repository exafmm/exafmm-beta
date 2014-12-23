#include "utils.h"

void encodeParticles(int N, float * X, float * min, float *max, uint64_t *mcodes, int maxlev) {
  uint32_t *codes = (uint32_t *)sakura_malloc(N, DIM * sizeof(uint32_t), "Hash code array");
  int nbins = (1 << maxlev);
  start_timer();
  compute_quantization_codes_TL(codes, X, N, nbins, min, max);
  stop_timer("Quantization");
  start_timer();
  morton_encoding_T(mcodes, codes, N);
  stop_timer("Morton encoding");
  free(codes);
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
