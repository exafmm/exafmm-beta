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

int main(int argc, char** argv){

  // Time counting variables 
  struct Body *data, *permuted_data, *data2, *permuted_data2;

  //const bool printNow = false;
  const bool printNow = true;

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

  /* Allocate memory for the particles */  
  start_timer();
  float *X = (float *)sakura_malloc(N, LDIM*sizeof(float), 
				    "Particle array");
  uint64_t particle_memory = (uint64_t) N*LDIM*sizeof(float);
  float *X2 = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle array");
  particle_memory += (uint64_t) N*LDIM*sizeof(float);

  stop_timer("Data mem. alloc.");

  printf("%-20s:   %luMB\n", "Particle mem", particle_memory / ML);

  /* Set up the parameters */
#ifdef SMALL_DENSE
  int maxlev = 5; // Maximum level of the tree
  int maxheight = 5; // Maximum height of the tree (extra control for debugin)
  int nbins = (1 << maxlev); // maximum number of boxes at the leaf level
#elif DENSE
  int lv = (int)ceil(log2f((double)ceil( (double) N / (double) population_threshold )) / (double)3);

  int maxlev = lv; // Maximum level of the tree
  int maxheight = lv; // Maximum height of the tree (extra control for debugin)
  int nbins = (1 << maxlev); // maximum number of boxes at the leaf level
#else
  int maxlev = 20; // Maximum level of the tree
  int maxheight = 20; // Maximun height of the tree (extra control for debugin)
  int nbins = (1 << maxlev); // maximum number of boxes at the leaf level
#endif


  /* Independant iteractions */
  for (int it=0; it<repeat; it++) {
    
    printf("Run: %d\n", it);
    
    /* Generate a 3-dimensional data distribution */
    start_timer();
    create_dataset_TL(X, N, dist);
    create_dataset_TL(X2, N, dist);
    stop_timer("Create data");
    
    /* Compute the bounds of the computational cube */
    
    start_timer();
    float min[DIM], max[DIM]; // the minimum and maximum in each dimension 
    
    //space_bounds(min, max, X, N);
    
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
    
#ifdef SMALL_DENSE
    uint16_t *particle_codes;
    uint16_t *bit_map;
    uint16_t *particle_codes2;
    uint16_t *bit_map2;
#elif DENSE
    uint32_t *particle_codes;
    uint32_t *bit_map;
    uint32_t *particle_codes2;
    uint32_t *bit_map2;
#else
    uint64_t *particle_codes;
    uint32_t *bit_map;
    uint64_t *particle_codes2;
    uint32_t *bit_map2;
#endif
    uint32_t *permutation_vector;
    uint32_t *permutation_vector2;
    memalloc_encoding((void**)&particle_codes, N);
#ifdef INTERLEAVE
    start_timer();
#endif

    encodeParticles(N, X, min, 
		    max, particle_codes, 
		    maxlev);

#ifdef INTERLEAVE
    stop_timer("Encoding");
#endif

    /* Spatial decomposition */
    memalloc_decomposeSpace(&permutation_vector, (void**)&bit_map, N);
   
#ifdef INTERLEAVE
    start_timer();
#endif
    
    decomposeSpace(N, (void **)&particle_codes, 
		   permutation_vector, (void*)bit_map, &X,
		   maxlev, population_threshold, dist);

#ifdef INTERLEAVE
    stop_timer("Tree building");
#endif

#ifndef SMALL_DENSE
#ifdef INTERLEAVE
cilk_spawn
#endif
    relocateParticles(N, &X, permutation_vector);
#endif

#ifdef INTERLEAVE
    start_timer();
#endif

    memalloc_encoding((void**)&particle_codes2, N);
    encodeParticles(N, X2, min, max, particle_codes2, maxlev);
    memalloc_decomposeSpace(&permutation_vector2, (void**)&bit_map2, N);
    decomposeSpace(N, (void**)&particle_codes2, 
		   permutation_vector2, (void*)bit_map2, &X2,
		   maxlev, population_threshold, dist);
#ifndef SMALL_DENSE
#ifdef INTERLEAVE
    cilk_sync;
    cilk_spawn 
#endif
    relocateParticles(N, &X2, permutation_vector2);
#endif

    /* Tree data structure. Parent to children connection */
    /* data structure */

#ifndef DENSE
    int nodes_per_level[20];
    int **node_pointers = (int **)malloc(maxlev*sizeof(int *)); 
    int **num_children = (int **)malloc(maxlev*sizeof(int *)); 
    int **children_first = (int **)malloc(maxlev*sizeof(int *)); 
#ifndef MORTON_ONLY
    int **node_codes = (int **)malloc(maxlev*sizeof(int *));
#else
    uint64_t **node_codes = (uint64_t **)malloc(maxlev*sizeof(uint64_t *)); 
#endif

    int nodes_per_level2[20];
    int **node_pointers2 = (int **)malloc(maxlev*sizeof(int *)); 
    int **num_children2 = (int **)malloc(maxlev*sizeof(int *)); 
    int **children_first2 = (int **)malloc(maxlev*sizeof(int *)); 
#ifndef MORTON_ONLY
    int **node_codes2 = (int **)malloc(maxlev*sizeof(int *));
#else
    uint64_t **node_codes2 = (uint64_t **)malloc(maxlev*sizeof(uint64_t *)); 
#endif

#ifndef DESNE
    int height = tree_formation((void *)bit_map, (void *)particle_codes, 
				nodes_per_level, node_pointers, 
				num_children, children_first, 
				(void**)node_codes, maxlev, N);
    int height2 = tree_formation((void *)bit_map2, (void *)particle_codes2, 
				 nodes_per_level2, node_pointers2, 
				 num_children2, children_first2, 
				 (void**)node_codes2, maxlev, N);
#endif

#endif
      
    /* The interaction list structure */

#ifdef DENSE
    int common_stencil[DIM*CN] = {0};
    int far_stencil[8*DIM*FN] = {0};
    int near_stencil[8*DIM*NN] = {0};

    int nodes_per_level[20]; // Assume a maximum of 20 levels for the tree
    int height =  maxlev;
    
    for(int i=0; i<height; i++){
      nodes_per_level[i] = ( 1 << (3*(i+1)) );
    }


#else
    int **clgs_link_list = (int **)malloc(height*sizeof(int *)); 
    int **nn_link_list = (int **)malloc(height*sizeof(int *)); 
    int **common_list = (int **)malloc(height*sizeof(int *));
    uint32_t **nn_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
    uint32_t **clgs_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
    uint32_t **common_count = (uint32_t **)malloc(height*sizeof(uint32_t *));
#endif


#ifndef DENSE


    form_interaction_lists(node_codes, children_first,
			   node_codes2, children_first2, 
			   nn_count, 
			   clgs_count, 
			   common_count,
			   nn_link_list, 
			   clgs_link_list,
			   common_list,
			   NULL,
			   NULL,
			   NULL,
			   node_pointers, 
			   nodes_per_level, 
			   nodes_per_level2, 
			   height, 
			   height2, 
			   N);

#else

    generate_interaction_stencil(common_stencil, far_stencil, 
				 near_stencil);
#endif


#ifdef INTERLEAVE
    cilk_sync;
    stop_timer("Relocation and formation");
#endif 


#ifdef DENSE 
    printf("Tree height: %d\n", maxlev);
#else
    printf("Tree height: %d\n", height);
#endif

    /* Verification part */
#ifndef DENSE

    verify_all(node_pointers, 
	       node_pointers2,
	       children_first, 
	       children_first2,
	       nodes_per_level, nodes_per_level2,
	       bit_map, bit_map2,
	       clgs_link_list,
	       nn_link_list,
	       common_list,
	       nn_count,
	       clgs_count,
	       common_count,
	       height, height2, N);

#else 
    verify_dense(bit_map, near_stencil, 
		 far_stencil, common_stencil, 
		 nodes_per_level,
		 N, maxheight);

#endif


#ifndef DENSE


    free_interaction_list_memo(nn_count, clgs_count, 
			       common_count, clgs_link_list,
			       nn_link_list, common_list,
			       height);
    free_tree_struct(node_pointers, num_children,
		     children_first, (void **)node_codes, height);


 free_tree_struct(node_pointers2, num_children2,
		  children_first2, (void **)node_codes2, height2);
#endif


 free(bit_map);
 free(particle_codes);
 free(permutation_vector);
 free(bit_map2);
 free(particle_codes2);
 free(permutation_vector2);
  }
  free(X);
  free(X2);
}








