/*

A simple program that demonstrate building of an octree using cilk.
This programm test the transpose version

author: Nikos Sismanis
date: Jul 2014

*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "timer.h"
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

void formInteractionStencil(int *common_stencil, int *far_stencil, int *near_stencil){
  int toy_parent[DIM];
  toy_parent[0] = 1; toy_parent[1] = 1; toy_parent[2] = 1;
  interaction_list_stencil(common_stencil, far_stencil, near_stencil, toy_parent);
}

void memalloc_encoding(void **mcodes, int N){
  *mcodes = sakura_malloc(N, sizeof(uint64_t), "Morton code array");
  uint64_t physical_memory = N*sizeof(uint64_t);
  printf("%-20s:   %luMB\n", "Encoding Phy mem", physical_memory / ML);
}

void memfree_encoding(void *mcodes){
  free(mcodes);
}

void encodeParticles(int N, float * X, float * min, float *max, void *particle_codes, int maxlev) {
  uint64_t *mcodes = (uint64_t *)particle_codes;
  uint32_t *codes = (uint32_t *)sakura_malloc(N, DIM * sizeof(uint32_t), "Hash code array");
  uint64_t working_memory = (uint64_t)DIM * N * sizeof(uint32_t);
  int nbins = (1 << maxlev);
  start_timer();
  compute_quantization_codes_TL(codes, X, N, nbins, min, max);
  stop_timer("Quantization");
  start_timer();
  morton_encoding_T(mcodes, codes, N, maxlev);
  stop_timer("Morton encoding");
  free(codes);
  printf("%-20s:   %luMB\n", "Encoding work mem", working_memory / ML);
}

void memalloc_decomposeSpace(uint32_t **permutation_vector, void **bit_map, int N){
  *bit_map = sakura_calloc(N, sizeof(uint32_t), "Bit map");
  uint64_t physical_memory = N*sizeof(uint32_t);
  *permutation_vector = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), 
						  "Permutation vector");
  physical_memory += N*sizeof(uint32_t);
  printf("%-20s:   %luMB\n", "Decom. Phy mem", physical_memory / ML);
}

void free_decomposeSpace(uint32_t *permutation_vector, void *bit_map ){
  free(permutation_vector);
  free(bit_map);
}

void decomposeSpace(int N, void **particle_codes, 
		    uint32_t *permutation_vector, void *bin_rep, float **X,
		    int maxlev, int population_threshold, int dist) {
  uint64_t *mcodes = (uint64_t *)*particle_codes;
  uint64_t *scodes = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), 
					       "Code buffer array");
  uint32_t *bit_map = (uint32_t *)bin_rep;
  uint64_t working_memory = N*sizeof(uint64_t);
  uint32_t *index = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), 
					      "Index vector");
  working_memory = N*sizeof(uint32_t);
  float *Y;
  if(N <= SMALLTH){
    Y = (float*)sakura_malloc(N, LDIM*sizeof(float), "Particle buffer");
    working_memory = (uint64_t)N*LDIM*sizeof(float);
  }
  start_timer();
  build_tree(Y, *X, mcodes, scodes, permutation_vector, 
	     index, bit_map, N, maxlev, maxlev, 
	     population_threshold, dist);
  stop_timer("Tree building");
  printf("%-20s:   %luMB\n", "Decomp. work mem", working_memory / ML);
  *particle_codes = (void*)scodes;
  if(N <= SMALLTH){
    float * tmp;
    tmp = Y;
    Y = *X;
    *X = tmp;
    free(Y);
  }
  free(mcodes);
  free(index);
}

#ifndef LIBRARY

int tree_formation(void *binrep, void *particle_codes, 
		   int *nodes_per_level, int **node_pointers, 
		   int **num_children, int **children_first, 
		   void **codes, int maxlevel, int N){
  uint64_t physical_mem = 0;
  uint32_t *bit_map = (uint32_t *)binrep;
  uint64_t *scodes = (uint64_t *)particle_codes;

#ifndef MORTON_ONLY
  int **node_codes = (int **)codes;
#else
  uint64_t **node_codes = (uint64_t **)codes;
#endif

  int *nodes_block_first = (int*)sakura_malloc(NP3*maxlevel, sizeof(int), "Node block");
    start_timer();
    int height = count_bins_bitmap_wrapper(nodes_per_level, 
					   nodes_block_first, 
					   bit_map, N, maxlevel);
    stop_timer("Count nodes");
    for(int i=0; i<height; i++){
      node_pointers[i] = (int *)sakura_malloc(nodes_per_level[i], sizeof(int), 
					      "Array of pointers to data");
      physical_mem += nodes_per_level[i]*sizeof(int);
      num_children[i] = (int *)sakura_malloc(nodes_per_level[i],sizeof(int), 
					     "Number of children array");
      physical_mem += nodes_per_level[i]*sizeof(int);
      children_first[i] = (int *)sakura_calloc(nodes_per_level[i], sizeof(int), 
					       "Children first array");
      physical_mem += nodes_per_level[i]*sizeof(int);
#ifndef MORTON_ONLY
      node_codes[i] = (int *)sakura_malloc(3*nodes_per_level[i], sizeof(int), 
					   "Node hash codes");
      physical_mem += 3*nodes_per_level[i]*sizeof(int);
#else
      node_codes[i] = (uint64_t *)sakura_malloc(nodes_per_level[i], sizeof(uint64_t),
						"New code array");
      physical_mem += nodes_per_level[i]*sizeof(int);
#endif 
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
    printf("%-20s:   %luMB\n", "Tree phys mem", physical_mem / ML);
    free(nodes_block_first);
    return(height);
}

#endif


#ifdef MORTON_ONLY
void form_interaction_lists(uint64_t **node_codes, int **children_first,
			    uint64_t **node_codes2, int **children_first2, 
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

#else
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
#endif

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
  printf("%-20s:   %luMB\n", "Inter. list phys mem", (uint64_t)physical_memory / ML);

}

void free_interaction_list_memo(uint32_t **nn_count, uint32_t **clgs_count, 
				uint32_t **common_count, int **clgs_link_list,
				int **nn_link_list, int **common_list,
				int height){

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
}

void free_tree_struct(int **node_pointers, int **num_children,
		      int **children_first, void **node_codes, int height){


  for(int i=0; i<height; i++){
    free(node_pointers[i]);
    free(num_children[i]);
    free(children_first[i]);
    free(node_codes[i]);
  }

  free(node_pointers);
  free(num_children);
  free(children_first);
  free(node_codes);

}


void generate_interaction_stencil(int *common_stencil, int *far_stencil, 
				  int *near_stencil){

  int toy_parent[DIM];
  toy_parent[0] = 1; toy_parent[1] = 1; toy_parent[2] = 1;
  interaction_list_stencil(common_stencil, far_stencil, near_stencil, toy_parent);

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
	expansions[i] = (int *)sakura_malloc(nodes_per_level2[i],sizeof(int), 
					     "Node expansions");
      }
      int *leaf_populations = (int *)sakura_malloc(N, sizeof(int), 
						   "Leaf population array");
      leaf_populations[0:N] = 0;
      uint64_t numleaves = find_leaf_populations(leaf_populations, bit_map2, N);
      int charge = verify_tree_wrapper(expansions, children_first2, 
				       node_pointers2, leaf_populations, 
				       nodes_per_level2[0], 0, N);
      int ss = 0;
      for(int i=0; i<N; i++){
	ss += leaf_populations[i];
      }      
      printf("Population: %d\n", ss);

      printf("Tree %s\n", (charge) ? "PASS" : "FAIL");


      /* Interaction list verification */

#ifdef NO_SYMBOLIC 

      int pass = verify_interactions_compressed_so_symbolic_wrapper(expansions, 
								    children_first, 
								    nn_link_list, 
								    nn_count,
								    clgs_link_list, 
								    clgs_count,
								    common_list, 
								    common_count,
								    nodes_per_level[0], N, 
								    height);
      
#else

      int pass = verify_interactions_compressed_wrapper(expansions, children_first, 
							nn_count, nn_link_list, 
							clgs_count, clgs_link_list,
							common_count, common_list,
							nodes_per_level[0], N, height);
#endif

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

void verify_dense(uint32_t *bit_map, int *near_stencil, 
		  int *far_stencil, int *common_stencil, 
		  int *nodes_per_level,
		  int N, int height){
      int **expansions = (int **)malloc(height*sizeof(int *));
      for(int i=0; i<height; i++){
	expansions[i] = (int *)sakura_malloc(nodes_per_level[i],sizeof(int), 
					     "Node expansions");
      }
      int *leaf_populations = (int *)sakura_malloc(N, sizeof(int), 
						   "Leaf population array");
      leaf_populations[0:N] = 0;
      find_leaf_populations(leaf_populations, bit_map, N);
      int charge = verify_tree_dense_wrapper(expansions, 
					     bit_map, leaf_populations, 
					     height, N);
      printf("Tree %s\n", (charge) ? "PASS" : "FAIL");
      int pass = verify_interactions_compressed_dense_wrapper(expansions, 
							      near_stencil, 
							      far_stencil, 
							      common_stencil,
							      N, height);
      printf("List %s\n", (pass) ? "PASS" : "FAIL");
      free(leaf_populations);
      for(int i=0; i<height; i++){
	free(expansions[i]);
      }
      free(expansions);
}

/* extenral funtion to be linked with exafmm */
#ifdef LIBRARY
void decomposeSpace(int N, void ** particle_codes, void ** particle_codes_buffer, 
		    uint32_t ** permutation_vector, uint32_t ** index, 
		    int maxlev) {
  float *X, *Y;
  uint32_t *tmp_id;
  //uint32_t *bit_map = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map 1");
  uint32_t *bit_map;

  uint64_t **mcodes = (uint64_t **)particle_codes;
  uint64_t **scodes = (uint64_t **)particle_codes_buffer;
  uint64_t *tmp;
    build_tree(Y, X, mcodes[0], scodes[0], permutation_vector[0],
	       index[0], bit_map, N, maxlev, maxlev, 64, 1);
    
    
#ifdef LIBRARY
    int st = (int)ceil((float)maxlev / (float)2);
    
    if((st & 1) == 0){
      tmp = mcodes[0];
      mcodes[0] = scodes[0];
      scodes[0] = tmp; 
      
      tmp_id = index[0];
      index[0] = permutation_vector[0];
      permutation_vector[0] = tmp_id;
    }
#endif
  
}
#endif

void decomposeSpacePermute(int N, float * Y, float * X, uint32_t * keys,
                           uint32_t *permutation_vector, int maxlev){

  bin_sort_dense_singlepass(Y, X, keys, permutation_vector, N, maxlev);

}


#ifdef DENSE
void decomposeSpacePermute(int N, float * Y, float * X, uint32_t * mcodes, 
			   uint32_t * scodes, uint32_t * permutation_vector, 
			   uint32_t * index, int maxlev, int population_threshold){
  uint32_t *bit_map = NULL;
  bin_sort_radix6_bitmap_small(Y, X, mcodes, scodes, permutation_vector,
			       index, bit_map,
			       N, 3*(maxlev-2),
			       0, 0, population_threshold);


}

#endif


void relocateParticles(int N, float **X, uint32_t *permutation_vector){
  float *Y  = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle buffer");
  uint64_t working_memo = (uint64_t)N*LDIM*sizeof(float);
  float *tmp;
  start_timer();
  if(N >= SMALLTH){
    rearrange_dataTL(Y, *X, permutation_vector, N);
  }
  stop_timer("Relocate particles");
  printf("%-20s:   %luMB\n", "Relocation work mem", working_memo / ML);
  tmp = Y;
  Y = *X;
  *X = tmp;
  free(Y);
}

#ifdef LIBRARY
void relocateParticles(int N, float * X, float * Y, uint32_t * permutation_vector) {
  cilk_spawn rearrange_dataTL(Y, X, permutation_vector, N);
}
#endif

#ifndef LIBRARY
/* main function */

#if 0
int main(int argc, char** argv){

  // Time counting variables 
  struct timeval startwtime, endwtime;
  struct Body *data, *permuted_data, *data2, *permuted_data2;

  //const bool printNow = false;
  const bool printNow = true;

  double memory_count = 0;
  double workspace_memory = 0;
  double workspace_memory_base = 0;
  double physical_memory = 0;
  double physical_memory_base = 0;
  double base_memory = 0;

  if (argc != 6) {
    printf("Usage: %s N dist pop rep P\n"
	   " where\n"
	   " N    : number of points\n"
	   " dist : distribution code (1-3)\n"
	   " pop  : population threshold\n"
	   " rep  : repetitions\n"
	   " P    : number of threads.\n", argv[0]);
    return (1);
  }

  // Input command line arguments
  int N = atoi(argv[1]); // Number of points
  int dist = atoi(argv[2]); // Distribution identifier 
  int population_threshold = atoi(argv[3]);
  int repeat = atoi(argv[4]);
  int nworkers = atoi(argv[5]);
  
  __cilkrts_set_param("nworkers",argv[5]); // the number of working threads
  printf("N = %d, T=%d\n", N, nworkers);

  start_timer();
  data = (Body *)sakura_malloc(N, sizeof(Body), "Particle array 1");
  permuted_data = (Body *)sakura_malloc(N, sizeof(Body), "Ordered particle data array 1");
  float *X = (float *)data;
  float *Y = (float *)permuted_data;
  uint32_t *mcodes = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), "Morton code array 1");
  uint32_t *scodes = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), "Sorted Morton code array 1");
  uint32_t *permutation_vector = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), 
							   "Permutation vector 1");
  uint32_t *index = (uint32_t *)sakura_malloc(N, sizeof(uint32_t), "Index vector 1");
  stop_timer("Memory allocation");

  /* Generate a 3-dimensional data distribution */
  start_timer();
  create_dataset_TL(X, N, dist);
  stop_timer("Create data");

  int lv = (int)ceil(log2f((double)ceil( (double) N / (double) population_threshold )) / (double)3);
  int maxlev = lv; // Maximum level of the tree

  /* Compute the bounds of the computational cube */

  start_timer();
  float min[DIM], max[DIM]; // the minimum and maximum in each dimension 
  for(int i=0; i<DIM; i++){
    min[i] = __sec_reduce_min(X[i:N:LDIM]);
    max[i] = __sec_reduce_max(X[i:N:LDIM]);
  }
  stop_timer("Box bounds");

  /* Independant iteractions */
  float * time = new float [repeat];
  for (int it=0; it<repeat; it++) {
    start_timer();
    //getKey2(N, X, min, max, mcodes, maxlev);
    stop_timer("Morton key");

    for(int i=0; i<N; i++){
      index[i] = i;
    }

    start_timer();
    //radixSort2(N, mcodes, scodes, permutation_vector, index, maxlev);
    stop_timer("Radix sort");
    start_timer();
    cilk_sync;
    stop_timer("Permutation");
    int nodes_per_level[20]; // Assume a maximum of 20 levels for the tree
    int height =  maxlev;
    for(int i=0; i<height; i++){
      nodes_per_level[i] = ( 1 << (3*(i+1)) );
    }
    
    int **node_codes = (int **)malloc(height*sizeof(int *)); // Cell array with the hash code of each node 
    int **node_pointers = (int **)malloc(height*sizeof(int *)); // Cell array with pointer to the data
    int **num_children = (int **)malloc(height*sizeof(int *)); // Cell array with the number of children of each node
    int **children_first = (int **)malloc(height*sizeof(int *)); // The position of the first child of each node.     
    start_timer();
    int **clgs_link_list = (int **)malloc(height*sizeof(int *)); // The link-list of the colleagues
    int **nn_link_list = (int **)malloc(height*sizeof(int *)); // The link-list of the neighbors
    int **common_list = (int **)malloc(height*sizeof(int *));
    uint32_t **nn_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); // Cell array that holds the number of neighbors each nodes has
    uint32_t **clgs_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); // Cell array that holds the number of colleagues each node has
    uint32_t **common_count = (uint32_t **)malloc(height*sizeof(uint32_t *));
    int common_stencil[DIM*CN] = {0};
    int far_stencil[8*DIM*FN] = {0};
    int near_stencil[8*DIM*NN] = {0};
    stop_timer("Malloc list");
    double datare_time = 0;
    double tmp_list_physical = 0;
    double tmp_list_workspace = 0;
      interaction_list_formation(node_codes, children_first, node_codes,
				 children_first, nn_count,
				 clgs_count, common_count,
				 nn_link_list, clgs_link_list, 
				 common_list, common_stencil, far_stencil, near_stencil,
				 node_pointers, nodes_per_level, nodes_per_level,
				 height, height, N, &memory_count, &workspace_memory, 
				 &physical_memory, &tmp_list_physical, 
				 &tmp_list_workspace);
      int **expansions = (int **)malloc(height*sizeof(int *));
      int **interactions = (int **)malloc(height*sizeof(int *));
      
      for(int i=0; i<height; i++){
	expansions[i] = (int *)sakura_malloc(nodes_per_level[i], sizeof(int), 
					     "Node expansions");
	interactions[i] = (int *)sakura_calloc(nodes_per_level[i], sizeof(int), 
					       "Node interactions");
      }
    
      int *leaf_populations = (int *)sakura_malloc(N, sizeof(int), 
						   "Leaf population array");
      leaf_populations[0:N] = 0;

      uint32_t * bit_map = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map 1");
      find_leaf_populations(leaf_populations, bit_map, N);

      int ss = 0;
      for(int i=0; i<N; i++){
	ss += leaf_populations[i];
      }      
      printf("sum = %d\n", ss);

      int charge = verify_tree_dense_wrapper(expansions, 
					     bit_map, leaf_populations, 
					     height, N);
      free(bit_map);
      
      printf("Tree %s\n", (charge) ? "PASS" : "FAIL");
      int pass = verify_interactions_compressed_dense_wrapper(expansions, 
							      near_stencil, 
							      far_stencil, 
							      common_stencil,
							      N, height);
      printf("List %s\n", (pass) ? "PASS" : "FAIL");
      
      for(int i=0; i<height; i++){
	free(expansions[i]);
	free(interactions[i]);
      }
      free(interactions);
      free(expansions);
      free(children_first);
      free(num_children);
      free(node_codes);
      free(node_pointers);
      free(nn_count);
      free(clgs_count);
      free(clgs_link_list);
      free(nn_link_list);
  }
  free(X);
  free(Y);
  free(scodes);
  free(permutation_vector);
  free(index);
}
#endif
#endif
