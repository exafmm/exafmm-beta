#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"
#include "utils.h"
#include <string.h>

#define NN 26
#define FN 37
#define CN 152

#ifndef LDIM
#define LDIM 12
#endif
#define DIM 3

void decode_morton_code(int *x, int *y, int *z, uint64_t mcode);

uint32_t mortonEncode_magicbits(uint32_t x, uint32_t y, uint32_t z);

int cmpfunc(const void *a, const void *b){
   return ( *(int *)a - *(int *)b );
}

int check_index(uint32_t *index, int N){

  /* sort the permutation vector */
  qsort(index, N, sizeof(uint32_t), cmpfunc);

  /* Check if all indexes are present in the input vector */
  int count = 0;
  for(int i=0; i<N; i++){
    count += (index[i] == i);

  }

  return(count == N);

}

#ifdef SMALL_DENSE
uint64_t find_leaf_populations(int *populations, uint16_t *bit_map, int N){
#else
uint64_t find_leaf_populations(int *populations, uint32_t *bit_map, int N){
#endif

  int pointer = 0;
  uint64_t numleaves = 0;  

  for (int i=1; i<N; i++){

    if(bit_map[i]>0){
      int pop = i - pointer;
      populations[pointer] = pop;
      pointer = i;
      numleaves++;
    }

  }

  populations[pointer] = N - pointer;
  return(numleaves);
}

#ifdef SMALL_DENSE
uint64_t find_leaf_populations(float *Y, int *populations, uint16_t *bit_map, int N){
#else
uint64_t find_leaf_populations(float *Y,int *populations, uint32_t *bit_map, int N){
#endif

  int counter = 1;
  int pointer = 0;
  uint64_t numleaves = 0;

  for (int i=1; i<N; i++){

    if(bit_map[i]>0){
      int pop = counter;
      populations[pointer] = pop;
      //counter = Y[i*LDIM + 3];
      counter = 0;
      pointer = i;
      numleaves++;
    }
    counter += Y[i*LDIM+3];
  }

  populations[pointer] = counter;
  return(numleaves);
}



void cumsum(uint32_t *X, int N){

  for(int i=1; i<N; i++){
    X[i] += X[i-1];
  }

}

int verify_tree(int **expansions, int** edges, 
		int** node_pointers, int *leaf_populations, 
		int node_id, int level){
  
  int charge = 0;

  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];

  int isnotleaf = children_stop - children_start;

  if(isnotleaf==0){
    charge = leaf_populations[node_pointers[level][node_id]];
  }
  else{

    for(int i=children_start; i<children_stop; i++){
      charge += verify_tree(expansions, edges, node_pointers, 
			    leaf_populations, i, level+1); 
    }
    
  }

  expansions[level][node_id] = charge;

  return charge;

}

int verify_tree_wrapper(int **expansions, int **edges, 
			int **node_pointers, int *leaf_populations, 
			int nnodes, int levels, int N){

  int charge = 0;
  int pass = 0;

  for(int i=0; i<nnodes; i++){
    charge += verify_tree(expansions, edges, node_pointers, leaf_populations, i, 0); 
  }

  if(charge == N){
    pass = 1;
  }
  else{
    printf("charge: %d\n", charge);
  }
  
  return pass;
}

/* Dense tree verification */
int verify_tree_dense(int **expansions, 
		      uint32_t *pointers2data, int *leaf_populations, 
		      uint32_t node_id, int level, int leaf_level, int N){

  int charge = 0;

  if(level == leaf_level-1){ // The node is a leaf
    //charge = leaf_populations[pointers2data[node_id]];
    int stop = ( node_id+1 == ( 1<<(3*leaf_level) ) ) ? N : pointers2data[node_id+1];
    int start = pointers2data[node_id];

    charge = stop - start;

  }
  else{ // the node is not a leaf

    // loop over the children
    for(int i=0; i<8; i++){
      uint32_t child_id = 8*node_id + i;
      charge += verify_tree_dense(expansions, pointers2data, 
				  leaf_populations, child_id, 
				  level+1, leaf_level, N);

    }
  }

  expansions[level][node_id] = charge;
  return charge;
}


int verify_tree_dense_wrapper(int **expansions, 
			      uint32_t *pointers2data, int *leaf_populations, 
			      int leaf_level, int N){

  int charge = 0;
  int pass = 0;

  for(int i=0; i<8; i++){
    charge += verify_tree_dense(expansions, 
				pointers2data, leaf_populations,
				i, 0, leaf_level, N);
  }

  if(charge == N){
    pass = 1;
  }
  else{
    printf("charge: %d\n", charge);
  }

  return pass;

}


/* Interaction list verification - Recursive version */
int verify_interaction(int **expansion, int **edges, uint32_t **nn_first, 
		       int **nn_list, uint32_t **fn_first, int **fn_list, 
		       int node_id, int level, int parent_charge, int N, int tree_height){

  int pass = 1;
  int interactions = parent_charge;

  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;

  int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
  int nn_stop = nn_first[level][node_id];
  int numnn = nn_stop - nn_start;
  
  uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
  uint32_t fn_stop = fn_first[level][node_id];
  uint32_t numfn = fn_stop - fn_start;  


    if(numnn>0){

    for(int i=nn_start; i<nn_stop; i++){
      interactions += expansion[level][nn_list[level][i]];
    }

  }

  if(numfn>0){

    for(uint32_t i=fn_start; i<fn_stop; i++){
      interactions += expansion[level][(uint32_t)fn_list[level][i]];
    }

  }

  if(isnotleaf>0 && level<tree_height){

    for(int i=children_start; i<children_stop; i++){
      
      pass &= verify_interaction(expansion, edges, nn_first, 
				nn_list, fn_first, fn_list, 
				i, level+1, interactions, N, tree_height);
      
    }
  }
  else{
    //interactions += expansion[level][node_id];
    if(interactions == N){
      pass = 1;
    }
    else{
      pass = 0;
      //printf("inter: %d\n", interactions);
    }
  }

  return pass;

}


int verify_interactions_wrapper(int **expansion, int **edges, 
				uint32_t **nn_first, int **nn_list, 
				uint32_t **fn_first, int **fn_list, 
				int nnodes, int N, int tree_height){

  int pass = 1;


  /* First pass */
  for(int i=0; i<nnodes; i++){
    pass &= verify_interaction(expansion, edges, nn_first, 
			      nn_list, fn_first, fn_list, 
			      i, 0, 0, N, tree_height);
  }

  return pass;

}


/* Interaction list verification - Iterative version*/
int verify_interaction_iterative(int **expansion, int **interactions, 
				 int **edges, uint32_t **nn_first, 
				 int **nn_list, uint32_t **fn_first, int **fn_list, 
				 int *nodes_per_level, int N, int tree_height){
  int pass = 1;

  for(int level=0; level<tree_height; level++){
    int nbins = nodes_per_level[level];

    for(int node_id=0; node_id<nbins; node_id++){

      // Find the position of the children
      int children_stop = edges[level][node_id];
      int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
      // If there are no children the node is a leaf
      int isnotleaf = children_stop - children_start; 
      

      // Find the index and the number of the near neighbors
      // number of near neighbors is 0 the node dose not have near neighbor 
      // interactions
      int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
      int nn_stop = nn_first[level][node_id];
      int numnn = nn_stop - nn_start;
      

      // Find the index and the number of far neighbor interactions
      uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
      uint32_t fn_stop = fn_first[level][node_id];
      uint32_t numfn = fn_stop - fn_start;


      if(numnn>0){ // Simulate P2P interactions	
	for(int i=nn_start; i<nn_stop; i++){
	  interactions[level][node_id] += expansion[level][nn_list[level][i]];
	}	
      }
      
      if(numfn>0){ // Simulate M2M interactions
	for(uint32_t i=fn_start; i<fn_stop; i++){
	  interactions[level][node_id] += expansion[level][(uint32_t)fn_list[level][i]];
	}	
      }      

      if(isnotleaf>0 && level<tree_height){ // If this is not a leaf node	
	for(int i=children_start; i<children_stop; i++){ // Simulate M2L interactions
	  interactions[level+1][i] += interactions[level][node_id];
	}
      }
      else{
	if(interactions[level][node_id] == N){
	  pass &= 1;
	}
	else{
	  pass &= 0;
	}
      }
      
      
    }
    
  }

  return(pass);

}

 int verify_interactions_wrapper_iterative(int **expansion, int **interactions,
					   int **edges, 
					   uint32_t **nn_first, int **nn_list, 
					   uint32_t **fn_first, int **fn_list, 
					   int *nodes_per_level, int N, 
					   int tree_height){

  /* First pass */
  int pass = verify_interaction_iterative(expansion, interactions, 
				 edges, nn_first, 
				 nn_list, fn_first, fn_list, 
				 nodes_per_level, N, tree_height);

  return pass;

}

/* Interaction list verification - Iterative version - Single array */
int verify_interaction_iterative_singlearray(int **expansion, int *interactions, 
					     int **edges, uint32_t **nn_first, 
					     int **nn_list, uint32_t **fn_first, 
					     int **fn_list, 
					     int *nodes_per_level, int N, int tree_height){
  int pass = 1;

  int level = 0; // level counter
  for(int bb=0; bb<nodes_per_level[tree_height-1]; bb++){

    if(bb>=nodes_per_level[level]){
      level++;
    }
    int offset = (level==0) ? 0 : nodes_per_level[level-1];
    int node_id = bb - offset;

    // Find the position of the children
    int children_stop = edges[level][node_id];
    int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
    // If there are no children the node is a leaf
    int isnotleaf = children_stop - children_start; 

    // Find the index and the number of the near neighbors
    // number of near neighbors is 0 the node dose not have near neighbor 
    // interactions
    int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
    int nn_stop = nn_first[level][node_id];
    int numnn = nn_stop - nn_start;

    // Find the index and the number of far neighbor interactions
    uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
    uint32_t fn_stop = fn_first[level][node_id];
    uint32_t numfn = fn_stop - fn_start;
    
    if(numnn>0){ // Simulate P2P interactions	
      for(int i=nn_start; i<nn_stop; i++){
	interactions[bb] += expansion[level][nn_list[level][i]];
      }	
    }
    
    if(numfn>0){ // Simulate M2M interactions
      for(uint32_t i=fn_start; i<fn_stop; i++){
	interactions[bb] += expansion[level][(uint32_t)fn_list[level][i]];
      }	
    }      

    if(isnotleaf>0 && level<tree_height){ // If this is not a leaf node	

      int nl_offset = nodes_per_level[level]; 

      for(int i=children_start; i<children_stop; i++){ // Simulate M2L interactions
	interactions[i+nl_offset] += interactions[bb];
      }
    }
    else{
      if(interactions[bb] == N){
	pass &= 1;
      }
      else{
	pass &= 0;
      }
    }    
    
  }

  return(pass);
}


 int verify_interactions_wrapper_iterative_singlearray(int **expansion, int **interactions,
						      int **edges, 
						      uint32_t **nn_first, int **nn_list, 
						      uint32_t **fn_first, int **fn_list, 
						      int *nodes_per_level, int N, 
						      int tree_height){

   int *nodes_sum = (int *)malloc(tree_height*sizeof(int));
   nodes_sum[0] = nodes_per_level[0];
   for(int i=1; i<tree_height; i++){
     nodes_sum[i] = nodes_sum[i-1] + nodes_per_level[i]; 
   }

   int *tmp_interactions = (int *)calloc(nodes_sum[tree_height-1],sizeof(int));

  /* First pass */
  int pass = verify_interaction_iterative_singlearray(expansion, tmp_interactions, 
						     edges, nn_first, 
						     nn_list, fn_first, fn_list, 
						     nodes_sum, N, tree_height);

  return pass;

  free(nodes_sum);
  free(tmp_interactions);
}

/* Symmetric interaction list verification - Recursive */
void verify_interactions_symetric(int **expansion, int **interactions,
				  int **edges, uint32_t **nn_first, 
				  int **nn_list, uint32_t **fn_first, int **fn_list, 
				  int node_id, int level, int parent_charge, int N){


  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;

  int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
  int nn_stop = nn_first[level][node_id];
  int numnn = nn_stop - nn_start;
  
  uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
  uint32_t fn_stop = fn_first[level][node_id];
  uint32_t numfn = fn_stop - fn_start;


  if(numnn>0){
    
    for(int i=nn_start; i<nn_stop; i++){
      interactions[level][node_id] += expansion[level][nn_list[level][i]];
      interactions[level][nn_list[level][i]] += expansion[level][node_id];
    }
    
  }


  if(numfn>0){

    for(int i=fn_start; i<fn_stop; i++){
      interactions[level][node_id] += expansion[level][fn_list[level][i]];
      interactions[level][fn_list[level][i]] += expansion[level][node_id];
    }

  }

  if(isnotleaf>0){

    for(int i=children_start; i<children_stop; i++){
      
      verify_interactions_symetric(expansion, interactions,edges, nn_first, 
				   nn_list, fn_first, fn_list, 
				   i, level+1, 0, 
				   N);
      
    }
  }
  else{
    interactions[level][node_id] += expansion[level][node_id];
  }


}

int verify_interactions_symetric_secondpass(int **interactions,
					    int **edges, int node_id, 
					    int level, int parent_charge, int N){

  int pass = 1;
  int inter = parent_charge + interactions[level][node_id];

  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;


  if(isnotleaf>0){
    for(int i=children_start; i<children_stop; i++){
      
      pass &= verify_interactions_symetric_secondpass(interactions, edges, 
						     i, level+1, inter, 
						     N);
    }
  }
  else{

    if(inter == N){
      pass = 1;
    }
    else{
      //printf("level: %d, inter = %d\n", level, inter);
      pass = 0;
    }
  }

  return pass;
}

int verify_interactions_symetric_wrapper(int **expansion, int** interactions, 
					 int **edges, 
					 uint32_t **nn_first, int **nn_list, 
					 uint32_t **fn_first, int **fn_list, 
					 int nnodes, int N){

  int pass = 1;
  
  /* First pass */
  for(int i=0; i<nnodes; i++){
    verify_interactions_symetric(expansion, interactions, edges, nn_first, 
				 nn_list, fn_first, fn_list, 
				 i, 0, 0, N);
  }

  /* Second Pass */
  for(int i=0; i<nnodes; i++){
    pass &= verify_interactions_symetric_secondpass(interactions, edges, 
						   i, 0, 0, N);
  }

  return pass;

}

int verify_interactions_compressed(int **expansion, int **edges, uint32_t **nn_first, 
				   int **nn_list, uint32_t **fn_first, int **fn_list, 
				   uint32_t **common_first, int **common_list,
				   int node_id, int level, int parent_charge, 
				   int N, int tree_height){


  int pass = 1;
  int interactions = parent_charge;

  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;

  int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
  int nn_stop = nn_first[level][node_id];
  int numnn = nn_stop - nn_start;
  
  uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
  uint32_t fn_stop = fn_first[level][node_id];
  uint32_t numfn = fn_stop - fn_start;

  uint32_t common_start = (node_id==0) ? 0 : common_first[level][node_id-1];
  uint32_t common_stop = common_first[level][node_id];
  uint32_t numcomm = common_stop - common_start;


  if(numnn>0){

    for(int i=nn_start; i<nn_stop; i++){
      interactions += expansion[level][nn_list[level][i]];
    }

  }

  if(numfn>0){

    for(uint32_t i=fn_start; i<fn_stop; i++){
      interactions += expansion[level][(uint)fn_list[level][i]];
    }

  }

  if(numcomm>0){ // Add the interactions of the common neighbors
    
    for(uint32_t i=common_start; i<common_stop; i++){
      interactions += expansion[level+1][(uint)common_list[level][i]];
    }
    
  }

  if(isnotleaf>0 && level<tree_height){

    for(int i=children_start; i<children_stop; i++){
      
      pass &= verify_interactions_compressed(expansion, edges, nn_first, 
					    nn_list, fn_first, fn_list,
					    common_first, common_list,
					    i, level+1, interactions, N, tree_height);
      
    }
  }
  else{
    //interactions += expansion[level][node_id];
    if(interactions == N){
      pass = 1;
    }
    else{
      pass = 0;
      //printf("inter: %d\n", interactions);
    }
  }

  return pass;

}


int verify_interactions_compressed_wrapper(int **expansion, int **edges, 
					   uint32_t **nn_first, int **nn_list, 
					   uint32_t **fn_first, int **fn_list, 
					   uint32_t **common_first, int **common_list,
					   int nnodes, int N, int tree_height){

  int pass = 1;

  for(int i=0; i<nnodes; i++){
    pass &= verify_interactions_compressed(expansion, edges, nn_first, 
					  nn_list, fn_first, fn_list,
					  common_first, common_list,
					  i, 0, 0, N, tree_height);
  }
  
  return pass;

}


void verify_interactions_symetric_compressed(int **expansion, int **interactions,
					     int **edges, uint32_t **nn_first, 
					     int **nn_list, uint32_t **fn_first, 
					     int **fn_list,
					     uint32_t **common_first, int **common_list,
					     int node_id, int level, int parent_charge, 
					     int N){


  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;

  int nn_start = (node_id==0)? 0 : nn_first[level][node_id-1];
  int nn_stop = nn_first[level][node_id];
  int numnn = nn_stop - nn_start;
  
  uint32_t fn_start = (node_id==0) ? 0 : fn_first[level][node_id-1];
  uint32_t fn_stop = fn_first[level][node_id];
  uint32_t numfn = fn_stop - fn_start;

  uint32_t common_start = (node_id == 0) ? 0 : common_first[level][node_id-1];
  uint32_t common_stop = common_first[level][node_id];
  uint32_t numcomm = common_stop - common_start;

  if(numnn>0){
    
    for(int i=nn_start; i<nn_stop; i++){
      interactions[level][node_id] += expansion[level][nn_list[level][i]];
      interactions[level][nn_list[level][i]] += expansion[level][node_id];
    }
    
  }

  if(numfn>0){

    for(int i=fn_start; i<fn_stop; i++){
      interactions[level][node_id] += expansion[level][fn_list[level][i]];
      interactions[level][fn_list[level][i]] += expansion[level][node_id];
    }

  }

  if(numcomm>0){

    for(int i=common_start; i<common_stop; i++){
      for(int j=children_start; j<children_stop; j++){
	interactions[level+1][j] += expansion[level+1][common_list[level][i]];
	interactions[level+1][common_list[level][i]] += expansion[level+1][j];
      }
    }

  }

  if(isnotleaf>0){

    for(int i=children_start; i<children_stop; i++){
      
      verify_interactions_symetric_compressed(expansion, interactions,edges, nn_first, 
					      nn_list, fn_first, fn_list, 
					      common_first, common_list,
					      i, level+1, 0, 
					      N);
      
    }
  }
  else{
    interactions[level][node_id] += expansion[level][node_id];
  }

}

int verify_interactions_symetric_wrapper_compressed(int **expansion, int** interactions, 
						    int **edges, 
						    uint32_t **nn_first, int **nn_list, 
						    uint32_t **fn_first, int **fn_list, 
						    uint32_t **common_first, 
						    int **common_list,
						    int nnodes, int N){

  int pass = 1;
  
  /* First pass */
  for(int i=0; i<nnodes; i++){
    verify_interactions_symetric_compressed(expansion, interactions, edges, nn_first, 
					    nn_list, fn_first, fn_list, 
					    common_first, common_list,
					    i, 0, 0, N);
  }

  /* Second Pass */
  for(int i=0; i<nnodes; i++){
    pass &= verify_interactions_symetric_secondpass(interactions, edges, 
						   i, 0, 0, N);
  }

  return pass;

}

int verify_interactions_compressed_no_symbolic(int **expansion, 
					       int **edges,
					       int **nn_list,
					       uint32_t **nn_count,
					       int **fn_list,
					       uint32_t **fn_count,
					       int **common_list,
					       uint32_t **common_count,
					       int node_id, int level, 
					       int parent_charge, 
					       int N){

  int pass = 1;
  int interactions = parent_charge;

  int children_stop = edges[level][node_id];
  int children_start = (node_id==0) ? 0 : edges[level][node_id-1];
  int isnotleaf = children_stop - children_start;

  int num_nn = nn_count[level][node_id];
  int num_fn = fn_count[level][node_id];
  int num_cn = common_count[level][node_id];

  if(num_nn>0){ // If this is a leaf bin
    int nn_start = node_id*NN;
    int nn_stop = node_id*NN + num_nn;
    for(int i=nn_start; i<nn_stop; i++){
      int ii = nn_list[level][i];
      interactions += expansion[level][ii];
    }
  }

  if(num_fn){ 
    int fn_start = node_id*FN;
    int fn_stop = node_id*FN + num_fn;
    for(uint32_t i=fn_start; i<fn_stop; i++){
      int ii = fn_list[level][i];
      interactions += expansion[level][ii];
    }
  }

  if(isnotleaf > 0){
    int common_start = num_cn*CN;
    int common_stop = num_cn*CN + num_cn;
    for(uint32_t i=common_start; i<common_stop; i++){
      int ii = common_list[level][i];
      interactions += expansion[level+1][ii];
    }
  }

  if(isnotleaf > 0){

    for(int i=children_start; i<children_stop; i++){
      
      pass &= verify_interactions_compressed_no_symbolic(expansion,
							 edges,
							 nn_list,
							 nn_count,
							 fn_list,
							 fn_count,
							 common_list,
							 common_count,
							 i, level+1, interactions, 
							 N);
      
    }
  }
  else{
    //interactions += expansion[level][node_id];
    if(interactions == N){
      pass = 1;
      //printf("node: %d\n", node_id);
    }
    else{
      if(node_id<100){ 
	//printf("node: %d, level: %d, inter: %d\n", node_id, level, interactions);
      }
    }
  }

  return pass;

}

int verify_interactions_compressed_so_symbolic_wrapper(int **expansion, 
						       int **edges, 
						       int **nn_list, 
						       uint32_t **nn_count,
						       int **fn_list, 
						       uint32_t **fn_count,
						       int **common_list, 
						       uint32_t **common_count,
						       int nnodes, int N, 
						       int tree_height){

  int pass = 1;

  for(int i=0; i<nnodes; i++){
    pass &= verify_interactions_compressed_no_symbolic(expansion, 
						       edges,
						       nn_list,
						       nn_count,
						       fn_list,
						       fn_count,
						       common_list,
						       common_count,
						       i, 0, 0, N);

  }

  return pass;

}

#ifndef SMALL_DENSE
int verify_interactions_compressed_dense(int **expansion, 
					 int *near_stencil, 
					 int *far_stencil,
					 int *common_stencil,
					 int node_id, int level, int parent_charge, 
					 int N, int leaf_level){

  int pass = 1;
  int interactions = parent_charge;
  int bound = 1 << (level+1);

  int node[DIM] = {0};
  int parent[DIM] = {0};
  int code[DIM] = {0};

  int children_stop = (node_id+1)*8;
  int children_start = (node_id)*8;

  decode_morton_code(&node[0], &node[1], &node[2], node_id);

  parent[0] = 2 * (node[0] / 2);
  parent[1] = 2 * (node[1] / 2);
  parent[2] = 2 * (node[2] / 2);

  uint64_t loc = node_id & 0x07;

  //printf("%d: %d\n", level, bound);
  
  for(int i=0; i<FN; i++){ // Find the far neighbors
    code[0] = parent[0] + far_stencil[loc*DIM*FN + i*DIM + 0];
    if(code[0] >= 0 && code[0] < bound){
      code[1] = parent[1] + far_stencil[loc*DIM*FN + i*DIM + 1];
      if(code[1] >= 0 && code[1] < bound){
	code[2] = parent[2] + far_stencil[loc*DIM*FN + i*DIM + 2];
	if(code[2] >= 0 && code[2] < bound){
	  int neighbor = mortonEncode_magicbits(code[0], code[1], code[2]); 
	  interactions += expansion[level][neighbor];

	}
      }
    }
  }
  
  
  if(level == (leaf_level-1)){ // If leaves add up the near neighbor interactions

    for(int i=0; i<NN; i++){

      code[0] = parent[0] + near_stencil[loc*DIM*NN + i*DIM + 0];
      if(code[0] >= 0 && code[0] < bound){
	code[1] = parent[1] + near_stencil[loc*DIM*NN + i*DIM + 1];
	if(code[1] >= 0 && code[1] < bound){
	  code[2] = parent[2] + near_stencil[loc*DIM*NN + i*DIM + 2];
	  if(code[2] >= 0 && code[2] < bound){

	    int neighbor = mortonEncode_magicbits(code[0], code[1], code[2]); 
	    interactions += expansion[level][neighbor];

	  }
	}
      }

    }
  }
  
  
  if(level < (leaf_level-1) ){ // Take into account the common list

    for(int i=0; i<CN; i++){

      code[0] = 2*node[0] + common_stencil[i*DIM + 0];
      if(code[0] >= 0 && code[0] < 2*bound){
	code[1] = 2*node[1] + common_stencil[i*DIM + 1];
	if(code[1] >= 0 && code[1] < 2*bound){
	  code[2] = 2*node[2] + common_stencil[i*DIM + 2];
	  if(code[2] >= 0 && code[2] < 2*bound){
	    int neighbor = mortonEncode_magicbits(code[0], code[1], code[2]); 
	    interactions += expansion[level+1][neighbor];

	  }
	}
      }

    }

  }
  
  //printf("level: %d, node %d, leaf: %d\n", level, node_id, leaf_level);

  if(level<leaf_level-1){
    for(int i=children_start; i<children_stop; i++){
      pass &= verify_interactions_compressed_dense(expansion, 
						   near_stencil, 
						   far_stencil,
						   common_stencil,
						   i, level+1, interactions, 
						   N, leaf_level);
    }

  }
  else{
    interactions += expansion[level][node_id];
    if(interactions == N){
      pass = 1;
      //printf("Passed\n");
    }
    else{
      pass = 0;
      //printf("inter: %d, %d, %d, %d\n", interactions, node[0], node[1], node[2]);
      //printf("parent: %d, %d, %d\n", parent[0], parent[1], parent[2]);
      
    }
  }

  return pass;

}

int verify_interactions_compressed_dense_wrapper(int **expansion, 
						 int *near_stencil, 
						 int *far_stencil, 
						 int *common_stencil,
						 int N, 
						 int leaf_level){
 
  int pass = 1;

  for(int i=0; i<8; i++){
    pass &= verify_interactions_compressed_dense(expansion,
						 near_stencil, 
						 far_stencil,
						 common_stencil,
						 i, 0, 0, N, leaf_level);
  }
  
  return pass;

}
#endif

