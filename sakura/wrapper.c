#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "utils.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
/* additional includes */
#include <asm-generic/unistd.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>
#include <sys/ioctl.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

#define DIM 3
#define LDIM 12
#define NP3 64
#define NN 26
#define FN 37
#define CN 152
#define ML 1000000

void encodeParticles(int N, float * X, float * min, float *max, void *particle_codes, int maxlev) {
  uint64_t *mcodes = (uint64_t *)particle_codes;
  uint32_t *codes = (uint32_t *)sakura_malloc(N, DIM * sizeof(uint32_t), "Hash code array");
  int nbins = (1 << maxlev);
  start_timer();
  compute_quantization_codes_TL(codes, X, N, nbins, min, max);
  stop_timer("Quantization");
  start_timer();
  morton_encoding_T(mcodes, codes, N, maxlev);
  stop_timer("Morton encoding");
  free(codes);
}

void decomposeSpace(int N, void **particle_codes, 
		    uint32_t *permutation_vector, void *bin_rep, float **X,
		    int maxlev, int population_threshold, int dist) {
  uint64_t *mcodes = (uint64_t *)*particle_codes;
  uint64_t *scodes = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Code buffer array");
  uint32_t *bit_map = (uint32_t *)bin_rep;
  uint32_t *index = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), "Index vector");
  float *Y = NULL;
  start_timer();
  build_tree(Y, *X, mcodes, scodes, permutation_vector, 
	     index, bit_map, N, maxlev, maxlev, 
	     population_threshold, dist);
  stop_timer("Tree building");
  *particle_codes = (void*)scodes;
  free(mcodes);
  free(index);
}

int tree_formation(void *binrep, void *particle_codes, 
		   int *nodes_per_level, int **node_pointers, 
		   int **num_children, int **children_first, 
		   void **codes, int maxlevel, int N){
  uint32_t *bit_map = (uint32_t *)binrep;
  uint64_t *scodes = (uint64_t *)particle_codes;
  int **node_codes = (int **)codes;
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

void form_interaction_lists(int **node_codes, int **children_first,
			    int **node_codes2, int **children_first2, 
			    uint32_t (**restrict nn_count), 
			    uint32_t (**restrict clgs_count), 
			    uint32_t (**restrict common_count),
			    int (**restrict nn_link_list), 
			    int (**restrict clgs_link_list),
			    int (**restrict common_list),
			    int (*restrict common_stencil),
			    int (*restrict far_stencil),
			    int (*restrict near_stencil),
			    int (**restrict node_pointers), 
			    int (*restrict nodes_per_level), 
			    int (*restrict nodes_per_level2), 
			    int height, int height2, int N){
  double memory_count = 0;
  double workspace_memory = 0;
  double physical_memory = 0;
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
			     common_list, common_stencil, far_stencil, near_stencil,
			     node_pointers, nodes_per_level, nodes_per_level2, 
			     height, height2, N, &memory_count, &workspace_memory, 
			     &physical_memory, &tmp_list_physical, 
			     &tmp_list_workspace);
}

void verify_all(int **node_pointers, 
		int **node_pointers2, 
		int **children_first,
		int **children_first2, 
		int *nodes_per_level, int *nodes_per_level2,
		void *binrep, void *binrep2,
		int **clgs_link_list,
		int **nn_link_list,
		int **common_list,
		uint32_t **nn_count,
		uint32_t **clgs_count,
		uint32_t **common_count,
		int height, int height2, int N){
  uint32_t *bit_map = (uint32_t *)binrep;
  uint32_t *bit_map2 = (uint32_t *)binrep2;
  int **expansions = (int **)malloc(height2*sizeof(int *));
  for(int i=0; i<height2; i++){
    expansions[i] = (int *)sakura_malloc(nodes_per_level2[i],sizeof(int),"Node expansions");
  }
  int *leaf_populations = (int *)sakura_malloc(N, sizeof(int),"Leaf population array");
  leaf_populations[0:N] = 0;
  uint64_t numleaves = find_leaf_populations(leaf_populations, bit_map2, N);
  int charge = verify_tree_wrapper(expansions, children_first2, 
				   node_pointers2, leaf_populations, 
				   nodes_per_level2[0], 0, N);
  int ss = 0;
  for(int i=0; i<N; i++){
    ss += leaf_populations[i];
  }      
  printf("Tree %s\n", (charge) ? "PASS" : "FAIL");
  int pass = verify_interactions_compressed_wrapper(expansions, children_first, 
						    nn_count, nn_link_list, 
						    clgs_count, clgs_link_list,
						    common_count, common_list,
						    nodes_per_level[0], N, height);
  printf("List %s\n", (pass) ? "PASS" : "FAIL");
  uint64_t inter_list_edges = 0;
  uint64_t num_tree_nodes = 0;
  for(int i=0;i<height; i++){
    inter_list_edges += nn_count[i][nodes_per_level[i]-1] + 
      clgs_count[i][nodes_per_level[i]-1];
    num_tree_nodes += nodes_per_level[i];
  }
  printf("%-20s: %d\n", "Tree height", height);
  printf("%-20s: %lu\n", "Tree nodes", num_tree_nodes);
  printf("%-20s: %lu\n", "Tree leaves", numleaves);
  printf("%-20s: %lu\n", "Edges",inter_list_edges);
}

void relocateParticles(int N, float **X, uint32_t *permutation_vector){
  float *Y  = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle buffer");
  uint64_t working_memo = (uint64_t)N*LDIM*sizeof(float);
  float *tmp;
  start_timer();
  rearrange_dataTL(Y, *X, permutation_vector, N);
  stop_timer("Relocate particles");
  printf("%-20s:   %luMB\n", "Relocation work mem", working_memo / ML);
  tmp = Y;
  Y = *X;
  *X = tmp;
  free(Y);
}
