#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "utils.h"

int main(int argc, char** argv){
  int N = atoi(argv[1]);
  int dist = atoi(argv[2]);
  int population_threshold = atoi(argv[3]);
  int nworkers = atoi(argv[4]);
  __cilkrts_set_param("nworkers",argv[4]);
  printf("N = %d, T=%d\n", N, nworkers);
  start_timer();
  float *X = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle array");
  float *X2 = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle array");
  stop_timer("Data mem. alloc.");
  int maxlev = 20;
  start_timer();
  create_dataset_TL(X, N, dist);
  create_dataset_TL(X2, N, dist);
  stop_timer("Create data");
  start_timer();
  float min[DIM], max[DIM];
  for(int i=0; i<DIM; i++){
    min[i] = __sec_reduce_min(X[i:N:LDIM]);
    max[i] = __sec_reduce_max(X[i:N:LDIM]);
  }
  float min2[DIM], max2[DIM];
  for(int i=0; i<DIM; i++){
    min2[i] = __sec_reduce_min(X2[i:N:LDIM]);
    max2[i] = __sec_reduce_max(X2[i:N:LDIM]);
  }
  min[:] = MIN(min[:], min2[:]);
  max[:] = MAX(max[:], max2[:]);
  stop_timer("Box bounds");
  uint64_t *particle_codes = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Morton code array");
  uint32_t *bit_map = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map");
  uint32_t *permutation_vector = (uint32_t *)sakura_malloc(N, sizeof(uint32_t),"Permutation vector");
  uint64_t *particle_codes2 = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Morton code array");
  uint32_t *bit_map2 = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map");
  uint32_t *permutation_vector2 = (uint32_t *)sakura_malloc(N, sizeof(uint32_t),"Permutation vector");
  encodeParticles(N, X, min, max, particle_codes, maxlev);
  decomposeSpace(N, &particle_codes, permutation_vector, bit_map, maxlev, population_threshold);
  relocateParticles(N, &X, permutation_vector);
  encodeParticles(N, X2, min, max, particle_codes2, maxlev);
  decomposeSpace(N, &particle_codes2, permutation_vector2, bit_map2, maxlev, population_threshold);
  relocateParticles(N, &X2, permutation_vector2);

  int nodes_per_level[20];
  int **node_pointers = (int **)malloc(maxlev*sizeof(int *)); 
  int **num_children = (int **)malloc(maxlev*sizeof(int *)); 
  int **children_first = (int **)malloc(maxlev*sizeof(int *)); 
  int **node_codes = (int **)malloc(maxlev*sizeof(int *));
  int nodes_per_level2[20];
  int **node_pointers2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **num_children2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **children_first2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **node_codes2 = (int **)malloc(maxlev*sizeof(int *));
  int height = tree_formation(bit_map, particle_codes, 
			      nodes_per_level, node_pointers, 
			      num_children, children_first, 
			      node_codes, maxlev, N);
  int height2 = tree_formation(bit_map2, particle_codes2, 
			       nodes_per_level2, node_pointers2, 
			       num_children2, children_first2, 
			       node_codes2, maxlev, N);
  int **clgs_link_list = (int **)malloc(height*sizeof(int *)); 
  int **nn_link_list = (int **)malloc(height*sizeof(int *)); 
  int **common_list = (int **)malloc(height*sizeof(int *));
  uint32_t **nn_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
  uint32_t **clgs_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
  uint32_t **common_count = (uint32_t **)malloc(height*sizeof(uint32_t *));
  form_interaction_lists(node_codes, children_first, node_codes2, children_first2, 
			 nn_count, clgs_count, common_count,
			 nn_link_list, clgs_link_list, common_list,
			 nodes_per_level, nodes_per_level2, 
			 height);
  int **expansions = (int **)malloc(height2*sizeof(int *));
  int **interactions = (int **)malloc(height2*sizeof(int *));
  for(int i=0; i<height2; i++){
    expansions[i] = (int *)sakura_malloc(nodes_per_level2[i],sizeof(int),"Node expansions");
  }
  int *leaf_populations = (int *)sakura_malloc(N, sizeof(int),"Leaf population array");
  leaf_populations[0:N] = 0;
  uint64_t numleaves = find_leaf_populations(leaf_populations, bit_map2, N);
  int charge = verify_tree_wrapper(expansions, children_first2,
				   node_pointers2, leaf_populations,
				   nodes_per_level2[0], N);
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
  for(int i=0; i<height; i++){
    free(nn_count[i]);
    free(clgs_count[i]);
    free(clgs_link_list[i]);
    free(nn_link_list[i]);
    free(common_count[i]);
    free(common_list[i]);
  }

  free(nn_count);
  free(clgs_link_list);
  free(nn_link_list);
  free(common_count);
  free(common_list);
  for(int i=0; i<height; i++){
    free(node_pointers[i]);
    free(num_children[i]);
    free(children_first[i]);
    free(node_codes[i]);
  }
  for(int i=0; i<height2; i++){
    free(node_pointers2[i]);
    free(num_children2[i]);
    free(children_first2[i]);
    free(node_codes2[i]);
    free(expansions[i]);
    free(interactions[i]);
  }
  free(node_pointers);
  free(num_children);
  free(children_first);
  free(node_codes);
  free(node_pointers2);
  free(num_children2);
  free(children_first2);
  free(node_codes2);
  free(expansions);
  free(interactions);
  free(leaf_populations);
  free(bit_map);
  free(particle_codes);
  free(permutation_vector);
  free(bit_map2);
  free(particle_codes2);
  free(permutation_vector2);
  free(X);
  free(X2);
}
