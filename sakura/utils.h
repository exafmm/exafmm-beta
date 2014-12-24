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

void create_dataset(float *X, int N, int dist);
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
uint64_t find_leaf_populations(int *populations, uint32_t* bit_map, int N);
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
int upward_pass(int **expansions, int** c_count2,
		int** node_pointers2, int *leaf_populations,
		int node_id, int level);
