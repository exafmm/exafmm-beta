#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/time.h>

#define DIM 3
#define LDIM 12
#define NP3 64
#define uint uint32_t

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))

static double get_time(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)(tv.tv_sec+tv.tv_usec*1e-6);
}

static double timer;

static void start_timer(){
  timer = get_time();
}

static void stop_timer(char *event_name){
  printf("%-20s: %8.4lf s\n",event_name, get_time()-timer);
}

static void* sakura_malloc(size_t items, size_t size, char *message){
  void *ptr = malloc(items*size);
  if(ptr == 0){
    printf("Out of memory at %s\n", message);
  }
  return ptr;
}

static void* sakura_calloc(size_t items, size_t size, char* message){
  void *ptr = calloc(items, size);
  if(ptr == 0){
    printf("Out of memory %s\n", message);
  }
  return ptr;
}

void create_dataset_TL(float *X, int N, int dist);
void compute_quantization_codes_TL(uint32_t (*restrict codes), float (*restrict X), 
				   int N, int nbins, float (*restrict min), 
				   float (*restrict max));
void morton_encoding_T(uint64_t (*restrict mcodes), uint32_t (*restrict codes), int N);
void decode_morton_code(int *x, int *y, int *z, uint64_t mcode);
int count_bins_bitmap_wrapper(int *nodes_per_level, int *node_block_first,
			      uint32_t *bit_map, int N, int L);
void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					int (**restrict node_codes), 
					int (*restrict nodes_block_first),
					uint32_t (*restrict bit_map), 
					uint64_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev);
void first_child_position_wrapper(int **children_first, 
				  int **num_children, 
				  int *nodes_per_level, 
				  int L);
void build_tree(uint64_t (*restrict zcodes), 
		uint64_t (*restrict codes), 
		uint32_t (*restrict pointIds), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold);
void rearrange_dataTL(float (*restrict Y), float (*restrict X), 
		      uint (*restrict Ids), int N);
void interaction_list_formation(int **node_codes, int **children_first,
				int **node_codes2,int **children_first2, 
				uint32_t (**restrict nn_count), 
				uint32_t (**restrict clgs_count), 
				uint32_t (**restrict common_count), 
				int (**restrict nn_link_list), 
				int (**restrict clgs_link_list), 
				int (**restrict common_list),
				int (*restrict nodes_per_level), 
				int (*restrict nodes_per_level2), 
				int height,
			        double *interaction_list_physical, 
				double *interaction_list_workspace);
uint64_t find_leaf_populations(int *populations, uint32_t* bit_map, int N);
int verify_tree_wrapper(int **expansions, int** edges, 
			int** node_pointers, int *leaf_populations, 
			int nnodes, int N);
int verify_interactions_wrapper_iterative_singlearray(int **expansion, int **interactions,
						      int **edges,
						      uint32_t **nn_first, int **nn_list,
						      uint32_t **fn_first, int **fn_list,
						      uint32_t **common_first,
						      int **common_list,
						      int *nodes_per_level, int N,
						      int tree_height);
void encodeParticles(int N, float * X, float * min, 
		     float *max, uint64_t *particle_codes, int maxlev);
void decomposeSpace(int N, uint64_t **particle_codes,
		    uint32_t *permutation_vector, uint32_t *bit_map,
		    int maxlev, int population_threshold);
void relocateParticles(int N, float **X, uint32_t *permutation_vector);
int tree_formation(uint32_t *bit_map, uint64_t *particle_codes,
		   int *nodes_per_level, int **node_pointers, 
		   int **num_chuildren, int **children_first, 
		   int **codes, int maxlevel, int N);
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
			    int height);

