#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/time.h>

#define DIM 3
#define LDIM 12
#define NP3 64

struct Body{
  float x[3];
  float SRC;
  int IBODY;
  int IRANC;
  float WEIGHT;
  float T;
  float TRE[4];
};

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))

namespace {
  double get_time(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return double(tv.tv_sec+tv.tv_usec*1e-6);
  }

  double timer;

  void start_timer(){
    timer = get_time();
  }

  void stop_timer(char *event_name){
    printf("%-20s: %8.4lf s\n",event_name, get_time()-timer);
  }

  void* sakura_malloc(size_t items, size_t size, char *message){
    void *ptr = malloc(items*size);
    if(ptr == 0){
      printf("Out of memory at %s\n", message);
    }
    return ptr;
  }

  void* sakura_calloc(size_t items, size_t size, char* message){
    void *ptr = calloc(items, size);
    if(ptr == 0){
      printf("Out of memory %s\n", message);
    }
    return ptr;
  }
}

void create_dataset_TL(float *X, int N, int dist);
void compute_quantization_codes_TL(uint32_t (*restrict codes), float (*restrict X), 
				   int N, int nbins, float (*restrict min), 
				   float (*restrict max));
void morton_encoding_T(uint64_t (*restrict mcodes), 
		       uint32_t (*restrict codes), 
		       int N, int max_level);
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
void build_tree(float *Y, float *X, uint64_t (*restrict zcodes), 
		uint64_t (*restrict codes), 
		uint32_t (*restrict pointIds), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist);
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
				int (*restrict common_stencil),
				int (*restrict far_stencil),
				int (*restrict near_stencil),
				int (**restrict node_pointers), 
				int (*restrict nodes_per_level), 
				int (*restrict nodes_per_level2), 
				int height, int height2, int N, 
				double *memory_count, double *workspace_memory, 
				double *physical_memory, double *interaction_list_physical, 
				double *interaction_list_workspace);
uint64_t find_leaf_populations(int *populations, uint32_t* bit_map, int N);
int verify_tree_wrapper(int **expansions, int** edges, 
			int** node_pointers, int *leaf_populations, 
			int nnodes, int levels, int N);
int verify_interactions_compressed_wrapper(int **expansion, int **edges, 
					   uint **nn_first, int **nn_list, 
					   uint **fn_first, int **fn_list, 
					   uint **common_first, int **common_list,
					   int nnodes, int N, int tree_height);
void encodeParticles(int N, float * X, float * min, 
		     float *max, void *particle_codes, 
		     int maxlev);
void decomposeSpace(int N, void **particle_codes, 
		    uint32_t *permutation_vector, void *bit_map, float **X,
		    int maxlev, int population_threshold, int dist);
void relocateParticles(int N, float **X, uint32_t *permutation_vector);
int tree_formation(void *binrep, void *particle_codes,
		   int *nodes_per_level, int **node_pointers, 
		   int **num_chuildren, int **children_first, 
		   void **codes, int maxlevel, int N);
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
			    int height, int height2, int N);
void verify_all(int **node_pointers, int **node_pointers2, 
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
		int height, int height2, int N);

