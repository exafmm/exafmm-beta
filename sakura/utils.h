#include <stdint.h>

#define SMALLTH 500000

struct Body{
  float x[3];
  float SRC;
  int IBODY;
  int IRANC;
  float WEIGHT;
  float T;
  float TRE[4];
};

struct node{
  float mpexp; // A variable to hold the multipole expansion                                                                    
  int leaf;
  struct node *children[64]; // pointers to the children                                                                         
};

struct tree_node{
  int codes[3];
  int edge;
  int num_children;
};

struct fmmls{
  uint64_t *trg;
  uint64_t *src;
  fmmls *next[8];
};

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))

/* general functions */
void* sakura_malloc(size_t items, size_t size, char *message);
void* sakura_calloc(size_t items, size_t size, char *message);

/* datasets.c file */
void create_dataset_TL(float *X, int N, int dist);

/* quantization.c file */
#ifdef SMALL_DENSE
void quantize_and_encoding(uint16_t (*restrict mcodes),float (*restrict X), int N,
			   int nbins, float (*restrict min), float (*restrict max));
#endif

#ifdef SMALL_DENSE
void compute_quantization_codes_TL(uint16_t (*restrict codes), float (*restrict X), 
				   int N, int nbins, float (*restrict min), 
				   float (*restrict max));
#else
void compute_quantization_codes_TL(uint32_t (*restrict codes), float (*restrict X), 
				   int N, int nbins, float (*restrict min), 
				   float (*restrict max));
#endif

void space_bounds(float (*restrict min), float (*restrict max), 
		  float (*restrict X), int N);

/* build_tree_cilk.c file */
#ifdef SMALL_DENSE
void morton_encoding_T(uint16_t (*restrict mcodes), 
		       uint16_t (*restrict codes), 
		       int N, int max_level);
#elif DENSE
void morton_encoding_T(uint32_t (*restrict mcodes), 
		       uint32_t (*restrict codes), 
		       int N, int max_level);
#else
void morton_encoding_T(uint64_t (*restrict mcodes), 
		       uint32_t (*restrict codes), 
		       int N, int max_level);
#endif


void scan_colleagues(uint32_t *C, uint32_t *Y, uint32_t *X, int N);

void bin_sort_dense_singlepass(float *Y, float *X, uint32_t *keys, 
			       uint32_t *permutation_vector, int N, int maxLevel);

void bin_sort_radix6_bitmap(uint64_t (*restrict zcodes), uint64_t (*restrict codes), 
			    uint32_t (*restrict pointIds), uint32_t (*restrict index),
			    uint32_t (*restrict bit_map),
			    int N, int sft, int lv, int stop, 
			    int population_threshold);

#ifdef SMALL_DENSE
void bin_sort_radix6_bitmap_small(float (*restrict Y), float (*restrict X),
                                  uint16_t (*restrict zcodes), uint16_t (*restrict codes),
                                  uint32_t (*restrict pointIds), 
				  uint32_t (*restrict index),
                                  uint16_t (*restrict bit_map),
                                  int N, int sft, int lv, int stop,
                                  int population_threshold);
#elif DENSE
  void bin_sort_radix6_bitmap_small(float (*restrict Y), float (*restrict X),
                                    uint32_t (*restrict zcodes), 
				    uint32_t (*restrict codes),
                                    uint32_t (*restrict pointIds), 
				    uint32_t (*restrict index),
                                    uint32_t (*restrict bit_map),
                                    int N, int sft, int lv, int stop,
                                    int population_threshold);
#else
    void bin_sort_radix6_bitmap_small(float (*restrict Y), float (*restrict X),
				      uint64_t (*restrict zcodes), 
				      uint64_t (*restrict codes),
				      uint32_t (*restrict pointIds), 
				      uint32_t (*restrict index),
				      uint32_t (*restrict bit_map),
				      int N, int sft, int lv, int stop,
				      int population_threshold);
#endif



#ifdef SMALL_DENSE
int count_bins_bitmap_wrapper(int *nodes_per_level, int *node_block_first,
			       uint16_t *bit_map, int N, int L);
#else
int count_bins_bitmap_wrapper(int *nodes_per_level, int *node_block_first,
			      uint32_t *bit_map, int N, int L);
#endif

#ifdef SMALL_DENSE
void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					int (**restrict node_codes),
					int (*restrict nodes_block_first),
					uint16_t (*restrict bit_map), 
					uint16_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev);
#else
#ifndef MORTON_ONLY
void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					int (**restrict node_codes), 
					int (*restrict nodes_block_first),
					uint32_t (*restrict bit_map), 
					uint64_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev);
#else
void parent_children_connection_wrapper(int (**restrict node_pointers), 
					int (**restrict num_children),
					uint64_t (**restrict node_codes),
					int (*restrict nodes_block_first),
					uint32_t (*restrict bit_map), 
					uint64_t (*restrict leaf_morton_codes),
					int N, int L, int maxL, int maxlev);
#endif
#endif
void first_child_position_wrapper(int **children_first, 
				  int **num_children, 
				  int *nodes_per_level, 
				  int L);
void scan_neighbors(uint *nn_memory, 
		    uint *nn_first, uint *nn_count, 
		    int N);
void bin_sort_radix6_bitmap_plummer(uint long *zcodes, uint long* codes, 
				    uint *pointIds, uint* index,
				    uint *bit_map,
				    int N, int sft, int lv, int stop, 
				    int population_threshold);
void bin_sort_radix6_bitmap_old(uint long *zcodes, uint long* codes, 
				uint *pointIds, uint* index,
				uint *bit_map,
				int N, int sft, int lv, int stop, 
				int population_threshold);

void bin_sort_radix6_bitmap_small(float (*restrict Y), float (*restrict X), 
				  uint64_t (*restrict zcodes), 
				  uint64_t (*restrict codes), 
				  uint32_t (*restrict pointIds), 
				  uint32_t (*restrict index), 
				  uint32_t (*restrict bit_map),
				  int N, 
				  int sft, int lv, int stop, 
				  int population_threshold);


#ifdef SMALL_DENSE
void build_tree(float *Y, float *X, uint16_t (*restrict zcodes), 
		uint16_t (*restrict codes), 
		uint32_t (*restrict pointIds), uint32_t (*restrict index),
		uint16_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist);
#elif DENSE
void build_tree(float *Y, float *X, uint32_t (*restrict zcodes), 
		uint32_t (*restrict codes), 
		uint32_t (*restrict pointIds), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist);
#else
void build_tree(float *Y, float *X, uint64_t (*restrict zcodes), 
		uint64_t (*restrict codes), 
		uint32_t (*restrict pointIds), uint32_t (*restrict index),
		uint32_t (*restrict bit_map),
		int N, int maxlev, int maxheight, 
		int population_threshold, int dist);
#endif


/* data_rearrange.c file*/
void rearrange_dataTL(float (*restrict Y), float (*restrict X), 
		      uint (*restrict Ids), int N);

/* cilk_linked_list.c file */
#ifndef MORTON_ONLY
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
#else
void interaction_list_formation(uint64_t **node_codes, int **children_first,
				uint64_t **node_codes2,int **children_first2, 
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
#endif

void interaction_list_stencil(int *common_stencil, int *far_stencil, 
			      int *near_stencil, int *parent_code);

/* Verification */
int check_index(unsigned int *index, int N);
int verify_tree_wrapper(int **expansions, int** edges, 
			int** node_pointers, int *leaf_populations, 
			int nnodes, int levels, int N);
int verify_interactions_wrapper(int **expansion, int **edges, 
				uint **nn_first, int **nn_list, 
				uint **fn_first, int **fn_list, 
				int nnodes, int N, int tree_height);
#ifdef SMALL_DENSE
uint64_t find_leaf_populations(int *populations, uint16_t* bit_map, int N);
#else
uint64_t find_leaf_populations(int *populations, uint32_t* bit_map, int N);
#endif

#ifdef SMALL_DENSE
uint64_t find_leaf_populations(float *Y, int *populations, uint16_t* bit_map, int N);
#else
uint64_t find_leaf_populations(float *Y, int *populations, uint32_t* bit_map, int N);
#endif


void cumsum(uint* X, int N);
int verify_interactions_symetric_wrapper(int **expansion, int** interactions, 
					 int **edges, 
					 uint **nn_first, int **nn_list, 
					 uint **fn_first, int **fn_list, 
					 int nnodes, int N);
int verify_interactions_compressed_wrapper(int **expansion, int **edges, 
					   uint **nn_first, int **nn_list, 
					   uint **fn_first, int **fn_list, 
					   uint **common_first, int **common_list,
					   int nnodes, int N, int tree_height);
int verify_interactions_symetric_wrapper_compressed(int **expansion, int** interactions, 
						    int **edges, 
						    uint **nn_first, int **nn_list, 
						    uint **fn_first, int **fn_list, 
						    uint **common_first, int **common_list,
						    int nnodes, int N);
void interaction_list_compressed_expanded_driver(int **clgs_link_list, 
						 uint **clgs_count,
						 int **nn_link_list, 
						 uint**nn_count,
						 int **common_list, 
						 uint **common_count,
						 int **target_tree_nodes, 
						 int **target_tree_edges,
						 int **source_tree_nodes, 
						 int **source_tree_edges,
						 int *nodes_per_level_target, 
						 int *nodes_per_level_source, 
						 int maxlev, int operation);
int verify_tree_dense_wrapper(int **expansions, 
			      uint32_t *pointers2data, int *leaf_populations, 
			      int leaf_level, int N);
int verify_interactions_compressed_dense_wrapper(int **expansion, 
						 int *near_stencil, 
						 int *far_stencil, 
						 int *common_stencil,
						 int N, 
						 int leaf_level);
int verify_interactions_compressed_so_symbolic_wrapper(int **expansion, 
						       int **edges, 
						       int **nn_list, 
						       uint32_t **nn_count,
						       int **fn_list, 
						       uint32_t **fn_count,
						       int **common_list, 
						       uint32_t **common_count,
						       int nnodes, int N, 
						       int tree_height);
 int verify_interactions_wrapper_iterative(int **expansion, int **interactions,
					   int **edges, 
					   uint32_t **nn_first, int **nn_list, 
					   uint32_t **fn_first, int **fn_list, 
					   int *nodes_per_level, int N, 
					   int tree_height);
 int verify_interactions_wrapper_iterative_singlearray(int **expansion, int **interactions,
						      int **edges, 
						      uint32_t **nn_first, int **nn_list, 
						      uint32_t **fn_first, int **fn_list, 
						      int *nodes_per_level, int N, 
						      int tree_height);
void generate_interaction_stencil(int *common_stencil, int *far_stencil, 
				  int *near_stencil);

/* Wrappers */
void memalloc_encoding(void **mcodes, int N);

void encodeParticles(int N, float * X, float * min, 
		     float *max, void *particle_codes, 
		     int maxlev);

void memalloc_decomposeSpace(uint32_t **permutation_vector, void **bit_map, int N);

void decomposeSpace(int N, void **particle_codes, 
		    uint32_t *permutation_vector, void *bit_map, float **X,
		    int maxlev, int population_threshold, int dist);
void relocateParticles(int N, float **X, uint32_t *permutation_vector);
int tree_formation(void *binrep, void *particle_codes,
		   int *nodes_per_level, int **node_pointers, 
		   int **num_chuildren, int **children_first, 
		   void **codes, int maxlevel, int N);
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
			    int height, int height2, int N);
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
			    int height, int height2, int N);
#endif
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
void verify_dense(uint32_t *bit_map, int *near_stencil, 
		  int *far_stencil, int *common_stencil, 
		  int *nodes_per_level,
		  int N, int height);
void free_interaction_list_memo(uint32_t **nn_count, uint32_t **clgs_count, 
				uint32_t **common_count, int **clgs_link_list,
				int **nn_link_list, int **common_list,
				int height);
void free_tree_struct(int **node_pointers, int **num_children,
		      int **children_first, void **node_codes, int height);

