/* 
This file contain the functions used in for the formation of the FMM linked list

author: Nikos Sismanis
date: Jul 2014

*/

#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"
#include "timer.h"
//#include "trace.h"
#include "utils.h"
#include <string.h>
#include <stdint.h>

#define MASK 0x07

#define MASK_X 0x04 
#define MASK_Y 0x02
#define MASK_Z 0x01

#define DIM 3

#define NN 26
#define FN 37
#define CN 152

#define X_3 0x09
#define Y_3 0x12
#define Z_3 0x24

#define X_2 0x08
#define Y_2 0x10
#define Z_2 0x20

#define FMASK 0x0FFFFFFF
#define mask_X 0x6DB6DB6DB6DB6DB6
#define mask_Y 0xDB6DB6DB6DB6DB6D
#define mask_Z 0xB6DB6DB6DB6DB6DB


typedef void (*base_function)(uint32_t*, int*, int, int); 

void decode_morton_code(int*, int*, int*, uint64_t);
uint64_t mortonEncode_magicbits(uint32_t , uint32_t , uint32_t);

void increase_counter(uint32_t (*restrict count), int (*restrict link_list), int target, int source){
  count[target]++;	
}

void store_pointer(uint32_t (*restrict count), int (*restrict link_list), int target, int source){
  uint32_t cursor = count[target];
  link_list[cursor] = source;
  cursor++;
  count[target] = cursor;
}

void interaction_list(int **clgs_link_list, 
		      uint32_t **clgs_count, 
		      int **nn_link_list, 
		      uint32_t **nn_count,
		      int **target_tree_nodes, 
		      int **target_tree_edges,
		      int **source_tree_nodes, 
		      int **source_tree_edges,
		      int target, int source, int level_target, 
		      int level_source, int maxlev, base_function seleced_op);

void triangular_co_split_interaction(int **clgs_link_list, 
				     uint32_t **clgs_count, 
				     int **nn_link_list, 
				     uint32_t **nn_count,
				     int **target_tree_nodes, int** target_tree_edges,
				     int **source_tree_nodes, int** source_tree_edges,
				     int target, int source, int level_target, 
				     int level_source, int maxlev, 
				     int source_children_stop, base_function selected_op){

  for(int i=target; i<source_children_stop; i++){
    interaction_list(clgs_link_list, clgs_count, 
		     nn_link_list, nn_count,
		     target_tree_nodes, target_tree_edges,
		     source_tree_nodes, source_tree_edges,
		     target, i, level_target, 
		     level_source+1, maxlev, selected_op);
    
  }
  
}


void interaction_list(int **clgs_link_list, 
		      uint32_t **clgs_count, 
		      int **nn_link_list, 
		      uint32_t **nn_count,
		      int **target_tree_nodes, 
		      int **target_tree_edges,
		      int **source_tree_nodes, 
		      int **source_tree_edges,
		      int target, int source, int level_target, 
		      int level_source, int maxlev, base_function selected_op){

  // Check if target and sources are leaves
  int target_children_stop = target_tree_edges[level_target][target];
  int target_children_start = (target==0) ? 0 : target_tree_edges[level_target][target-1];
  int source_children_stop = source_tree_edges[level_source][source];
  int source_children_start = (source == 0) ? 0 : source_tree_edges[level_source][source-1];
  
  if(level_target==level_source){ // The target and source are at the same level

    int Target_code_X = target_tree_nodes[level_target][DIM*target];
    int Target_code_Y = target_tree_nodes[level_target][DIM*target+1];
    int Target_code_Z = target_tree_nodes[level_target][DIM*target+2];
    int Source_code_X = source_tree_nodes[level_source][DIM*source];
    int Source_code_Y = source_tree_nodes[level_source][DIM*source+1];
    int Source_code_Z = source_tree_nodes[level_source][DIM*source+2];

    // compute the difference between dimension codes
    int abs_diff_x = abs(Target_code_X - Source_code_X);
    int abs_diff_y = abs(Target_code_Y - Source_code_Y);
    int abs_diff_z = abs(Target_code_Z - Source_code_Z);
    
    int are_neighbors = (abs_diff_x <= 1) && (abs_diff_y <= 1) && (abs_diff_z <= 1);
    int are_same_bin = (abs_diff_x == 0) && (abs_diff_y == 0) && (abs_diff_z == 0);
    int src_isnot_leaf = (((source_children_stop - source_children_start) > 0) && 
			  (level_source < maxlev));
    int trg_isnot_leaf = (((target_children_stop - target_children_start) > 0) && 
			  (level_target < maxlev));

    if(are_same_bin && trg_isnot_leaf && src_isnot_leaf){
	for(int i=target_children_start; i<target_children_stop; i++){
	  cilk_spawn interaction_list(clgs_link_list, 
				      clgs_count,
				      nn_link_list, 
				      nn_count,
				      target_tree_nodes, target_tree_edges,
				      source_tree_nodes, source_tree_edges,
				      i, source, level_target+1,
				      level_source, maxlev, selected_op);
	}
      //cilk_sync;

    }
    else if(are_neighbors && !are_same_bin){
       
      if(!src_isnot_leaf || !trg_isnot_leaf){
	selected_op(nn_count[level_target], 
		    nn_link_list[level_target], target, source);
      }
      else{

	  for(int i=target_children_start; i<target_children_stop; i++){ // recursive call  
	    cilk_spawn interaction_list(clgs_link_list, 
					clgs_count, 
					nn_link_list, 
					nn_count,
					target_tree_nodes, target_tree_edges,
					source_tree_nodes, source_tree_edges,
					i, source, level_target+1, 
					level_source, maxlev, selected_op);
	    
	  }
	  //cilk_sync;
      }
    }
    else{
      selected_op(clgs_count[level_target], 
		  clgs_link_list[level_target], target, source);
    }

  }
  else{ // The nodes are not at the same level

    if(level_source > level_target){ // The target is larger
	for(int i=target_children_start; i<target_children_stop; i++){ // recursive call  
	  cilk_spawn interaction_list(clgs_link_list, 
				      clgs_count, 
				      nn_link_list, 
				      nn_count,
				      target_tree_nodes, target_tree_edges,
				      source_tree_nodes, source_tree_edges,
				      i, source, level_target+1, 
				      level_source, maxlev, selected_op);
	  
	}
      //cilk_sync;
    }
    else{ // The source is larger
      for(int i=source_children_start; i<source_children_stop; i++){
	interaction_list(clgs_link_list, 
			 clgs_count, 
			 nn_link_list, 
			 nn_count,
			 target_tree_nodes, target_tree_edges,
			 source_tree_nodes, source_tree_edges,
			 target, i, level_target, 
			 level_source+1, maxlev, selected_op);	
      }
    }
    
  }

}

void interaction_list_wrapper(int **clgs_link_list,
			      uint32_t **clgs_count,
			      int **nn_link_list,
			      uint32_t **nn_count,
			      int **target_tree_nodes, 
			      int **target_tree_edges,
			      int **source_tree_nodes, 
			      int **source_tree_edges,
			      int target, int *nodes_per_level,
			      int level, int maxlev, base_function selected_op){

  for(int i=0; i<nodes_per_level[level]; i++){
    interaction_list(clgs_link_list, 
		     clgs_count,
		     nn_link_list,
		     nn_count,
		     target_tree_nodes, target_tree_edges,
		     source_tree_nodes, source_tree_edges,
		     target, i, level, 
		     level, maxlev, selected_op);
  }
}


void interaction_list_classical(int **clgs_link_list, 
				uint32_t **clgs_count,
				int **nn_link_list, 
				uint32_t **nn_count,
				int **target_tree_nodes, 
				int **target_tree_edges,
				int **source_tree_nodes, 
				int **source_tree_edges,
				int *nodes_per_level_target, 
				int *nodes_per_level_source, 
				int maxlev, int operation){


  int level = 0;

  base_function use_function = (operation==0) ? &increase_counter : &store_pointer;

  for(int i=0; i<nodes_per_level_target[level]; i++){
    
    cilk_spawn interaction_list_wrapper(clgs_link_list, 
					clgs_count,
					nn_link_list, 
					nn_count,
					target_tree_nodes, 
					target_tree_edges,
					source_tree_nodes, 
					source_tree_edges,
					i, nodes_per_level_source,
					level, maxlev-1, use_function);
    
  }
  cilk_sync;
  //writeTracer();
}

/* pran prefix of the colleagues count. 
   The purpose of this is to compute the required 
   memory for the colleagues list and find the position of 
   the first colleague for each tree node */
 void scan_colleagues(uint32_t (*restrict clgs_memory), uint32_t (*restrict clgs_first), 
		      int N){

  int offset = 0;
  for(int i=0; i<N; i++){
    int ss = clgs_first[i];
    clgs_first[i] = offset;
    offset += ss;
  }
  clgs_memory[0] = offset;

}


void interaction_list_compressed(int **clgs_link_list, 
				 uint32_t **clgs_count, 
				 int **nn_link_list, 
				 uint32_t **nn_count,
				 int **common_list, 
				 uint32_t **common_count,
				 int **target_tree_nodes, 
				 int **target_tree_edges,
				 int **source_tree_nodes, 
				 int **source_tree_edges,
				 int target, int source, int level_target, 
				 int level_source, int maxlev, 
				 base_function selected_op);


void triangular_co_split_interaction_compressed(int **clgs_link_list, 
						uint32_t **clgs_count, 
						int **nn_link_list, 
						uint32_t **nn_count,
						int **common_list,
						uint32_t **common_count,
						int **target_tree_nodes, 
						int** target_tree_edges,
						int **source_tree_nodes, 
						int** source_tree_edges,
						int target, int source, int level_target, 
						int level_source, int maxlev, 
						int source_children_stop, 
						base_function selected_op){

  for(int i=target; i<source_children_stop; i++){
    interaction_list_compressed(clgs_link_list, clgs_count, 
				nn_link_list, nn_count,
				common_list, common_count,
				target_tree_nodes, target_tree_edges,
				source_tree_nodes, source_tree_edges,
				target, i, level_target, 
				level_source+1, maxlev, selected_op);
    
  }
  
}

 
void intermidiate_call(int **clgs_link_list, 
		       uint32_t **clgs_count, 
		       int **nn_link_list, 
		       uint32_t **nn_count,
		       int **common_list, 
		       uint32_t **common_count,
		       int **target_tree_nodes, 
		       int **target_tree_edges,
		       int **source_tree_nodes, 
		       int **source_tree_edges,
		       int source, int level_target, 
		       int level_source, int maxlev, 
		       base_function selected_op,
		       int target_children_start,
		       int target_children_stop){


  for(int i=target_children_start; i<target_children_stop; i++){ // recursive call  
    interaction_list_compressed(clgs_link_list, 
				clgs_count, 
				nn_link_list, 
				nn_count,
				common_list, 
				common_count,
				target_tree_nodes, target_tree_edges,
				source_tree_nodes, source_tree_edges,
				i, source, level_target+1, 
				level_source, maxlev, 
				selected_op);
  }





}


void interaction_list_compressed(int **clgs_link_list, 
				 uint32_t **clgs_count, 
				 int **nn_link_list, 
				 uint32_t **nn_count,
				 int **common_list, 
				 uint32_t **common_count,
				 int **target_tree_nodes, 
				 int **target_tree_edges,
				 int **source_tree_nodes, 
				 int **source_tree_edges,
				 int target, int source, int level_target, 
				 int level_source, int maxlev, 
				 base_function selected_op){ 
  
   
  // Check if target and sources are leaves
  int target_children_stop = target_tree_edges[level_target][target];
  int target_children_start = (target==0) ? 0 : target_tree_edges[level_target][target-1];
  int source_children_stop = source_tree_edges[level_source][source];
  int source_children_start = (source == 0) ? 0 : source_tree_edges[level_source][source-1];
  

  if(level_target == level_source){

    int Target_code_X = target_tree_nodes[level_target][DIM*target];
    int Target_code_Y = target_tree_nodes[level_target][DIM*target+1];
    int Target_code_Z = target_tree_nodes[level_target][DIM*target+2];
    int Source_code_X = source_tree_nodes[level_source][DIM*source];
    int Source_code_Y = source_tree_nodes[level_source][DIM*source+1];
    int Source_code_Z = source_tree_nodes[level_source][DIM*source+2];

    // compute the difference between dimension codes
    int abs_diff_x = abs(Target_code_X - Source_code_X);
    int abs_diff_y = abs(Target_code_Y - Source_code_Y);
    int abs_diff_z = abs(Target_code_Z - Source_code_Z);
    
    int are_neighbors = (abs_diff_x <= 1) && (abs_diff_y <= 1) && (abs_diff_z <= 1);
    int are_same_bin = (abs_diff_x == 0) && (abs_diff_y == 0) && (abs_diff_z == 0);
    int src_isnot_leaf = (((source_children_stop - source_children_start) > 0) && 
			  (level_source < maxlev));
    int trg_isnot_leaf = (((target_children_stop - target_children_start) > 0) && 
			  (level_target < maxlev));
    if(are_same_bin && trg_isnot_leaf && src_isnot_leaf){
      cilk_for(int i=target_children_start; i<target_children_stop; i++){
	interaction_list_compressed(clgs_link_list, 
				      clgs_count,
				      nn_link_list, 
				      nn_count,
				      common_list, 
				      common_count,
				      target_tree_nodes, target_tree_edges,
				      source_tree_nodes, source_tree_edges,
				      i, source, level_target+1,
				      level_source, maxlev, selected_op);
	}
    } else if(are_neighbors && !are_same_bin){
      if(!src_isnot_leaf || !trg_isnot_leaf){
	selected_op(nn_count[level_target], 
		    nn_link_list[level_target], target, source);
      } else {
	  for(int i=source_children_start; i<source_children_stop; i++){ // recursive call  
	    interaction_list_compressed(clgs_link_list, 
					clgs_count, 
					nn_link_list, 
					nn_count,
					common_list, 
					common_count,
					target_tree_nodes, target_tree_edges,
					source_tree_nodes, source_tree_edges,
					target, i, level_target, 
					level_source+1, maxlev, selected_op);
	  }
      }
    } else {
      selected_op(clgs_count[level_target], 
		  clgs_link_list[level_target], target, source);
    }
  } else {
    if(level_source > level_target){ // The target is larger
      int Target_code_X = target_tree_nodes[level_target][DIM*target];
      int Target_code_Y = target_tree_nodes[level_target][DIM*target+1];
      int Target_code_Z = target_tree_nodes[level_target][DIM*target+2];
      int Source_code_X = source_tree_nodes[level_source][DIM*source];
      int Source_code_Y = source_tree_nodes[level_source][DIM*source+1];
      int Source_code_Z = source_tree_nodes[level_source][DIM*source+2];
      
      int diff_x = 2*Target_code_X - Source_code_X;
      int diff_y = 2*Target_code_Y - Source_code_Y;
      int diff_z = 2*Target_code_Z - Source_code_Z;
      
      int is_common_neighbor = (diff_x == 2 || diff_x == -3) ||
	(diff_y == 2 || diff_y == -3) ||
	(diff_z == 2 || diff_z == -3);
      
      
      if(is_common_neighbor){ // If S is a commpon neighbor of T
	
	selected_op(common_count[level_target], 
		    common_list[level_target], target, source);
		
      }
      
      else{

	intermidiate_call(clgs_link_list, 
			  clgs_count, 
			  nn_link_list, 
			  nn_count,
			  common_list, 
			  common_count,
			  target_tree_nodes, target_tree_edges,
			  source_tree_nodes, source_tree_edges,
			  source, level_target, 
			  level_source, maxlev, 
			  selected_op,
			  target_children_start,
			  target_children_stop);

      }
    }
    else{ // The source is larger
      for(int i=source_children_start; i<source_children_stop; i++){
	interaction_list_compressed(clgs_link_list, 
				    clgs_count, 
				    nn_link_list, 
				    nn_count,
				    common_list, 
				    common_count,
				    target_tree_nodes, target_tree_edges,
				    source_tree_nodes, source_tree_edges,
				    target, i, level_target, 
				    level_source+1, maxlev, selected_op);	
      }
    }


  }

}

int are_near_neighbors(uint64_t child_mcode, 
		       uint64_t parent_mcode ){

  uint64_t parent_l = parent_mcode << 3;

  uint64_t mask_child_x = ( child_mcode | mask_X );
  uint64_t mask_child_y = ( child_mcode | mask_Y );
  uint64_t mask_child_z = ( child_mcode | mask_Z );


  uint64_t mask_parent_x = ( parent_l | mask_X );
  uint64_t mask_parent_y = ( parent_l | mask_Y );
  uint64_t mask_parent_z = ( parent_l | mask_Z );

  
  int is_common_neighbor = ( (( mask_child_x + X_2 ) | mask_X)  == mask_parent_x ) || 
    ( (( mask_parent_x + X_3 ) | mask_X) == mask_child_x ) || ( (( mask_child_y + Y_2 ) | mask_Y) == mask_parent_y) ||
    ( (( mask_parent_y + Y_3 ) | mask_Y) == mask_child_y ) || ( (( mask_child_z + Z_2 ) | mask_Z) == mask_parent_z ) ||
    ( (( mask_parent_z + Z_3 ) | mask_Z)== mask_child_z);
  
  return(is_common_neighbor);

}

int are_near_neighbors( int (*restrict child),  
			int (*restrict parent) ){

  
  int diff_x = 2*parent[0] - child[0];
  int diff_y = 2*parent[1] - child[1];
  int diff_z = 2*parent[2] - child[2];
  
  int is_common_neighbor = (diff_x == 2 || diff_x == -3) ||
    (diff_y == 2 || diff_y == -3) ||
    (diff_z == 2 || diff_z == -3);
  

  /*
  if(is_common_neighbor != is_common_neighbor1){
    printf("problem %d, %d\n", is_common_neighbor, is_common_neighbor1);
    int code1[DIM] = {0};
    int code2[DIM] = {0};
    decode_morton_code(&code1[0], &code1[1], &code1[2], mask_parent_x);
    decode_morton_code(&code2[0], &code2[1], &code2[2], mask_child_x);
    printf("mas_parent_x: %lx, mask_child_x: %lx\n", mask_parent_x | mask_X , mask_child_x | mask_X);
    printf("parent_x: %d %d, %d %d, child_x %d %d, %d %d - %lx\n", code1[0], 2*parent[0], code1[1], code1[2], code2[0], child[0], code2[1], code2[2], (mask_X));

  }
  */
  
  return(is_common_neighbor);

}


#ifndef MORTON_ONLY
void interaction_list_compressed_expanded(int (**restrict clgs_link_list), 
					  uint32_t (**restrict clgs_count), 
					  int (**restrict nn_link_list), 
					  uint32_t (**restrict nn_count),
					  int (**restrict common_list), 
					  uint32_t (**restrict common_count),
					  int **target_tree_nodes, 
					  int **target_tree_edges,
					  int **source_tree_nodes, 
					  int **source_tree_edges,
					  int target, int source, int level_target, 
					  int level_source, int maxlev, 
					  base_function selected_op){
#else
void interaction_list_compressed_expanded(int (**restrict clgs_link_list), 
					  uint32_t (**restrict clgs_count), 
					  int (**restrict nn_link_list), 
					  uint32_t (**restrict nn_count),
					  int (**restrict common_list), 
					  uint32_t (**restrict common_count),
					  uint64_t **target_tree_nodes,
					  int **target_tree_edges,
					  uint64_t **source_tree_nodes,
					  int **source_tree_edges,
					  int target, int source, int level_target, 
					  int level_source, int maxlev, 
					  base_function selected_op){
#endif
   

  int have_neighbors_T[8] = {0};
  int have_neighbors_S[8] = {0};

  int target_children_stop = target_tree_edges[level_target][target];
  int target_children_start = (target==0) ? 0 : target_tree_edges[level_target][target-1];
  int source_children_stop = source_tree_edges[level_source][source];
  int source_children_start = (source == 0) ? 0 : source_tree_edges[level_source][source-1];

  int nchild_source = source_children_stop - source_children_start;
  int nchild_target = target_children_stop - target_children_start;

#ifndef MORTON_ONLY
  int Target_code_X = target_tree_nodes[level_target][DIM*target];
  int Target_code_Y = target_tree_nodes[level_target][DIM*target+1];
  int Target_code_Z = target_tree_nodes[level_target][DIM*target+2];
  int Source_code_X = source_tree_nodes[level_source][DIM*source];
  int Source_code_Y = source_tree_nodes[level_source][DIM*source+1];
  int Source_code_Z = source_tree_nodes[level_source][DIM*source+2];

  int abs_diff_x = abs(Target_code_X - Source_code_X);
  int abs_diff_y = abs(Target_code_Y - Source_code_Y);
  int abs_diff_z = abs(Target_code_Z - Source_code_Z);


#else

  uint64_t Target_mcode = target_tree_nodes[level_target][target];
  uint64_t Source_mcode = source_tree_nodes[level_source][source];

#endif
  
#ifndef MORTON_ONLY
  int are_same_bin = (abs_diff_x == 0) && (abs_diff_y == 0) && (abs_diff_z == 0);
#else
  int are_same_bin = Target_mcode == Source_mcode;
#endif
  int src_isnot_leaf = (((source_children_stop - source_children_start) > 0) && 
			(level_source < maxlev));
  int trg_isnot_leaf = (((target_children_stop - target_children_start) > 0) && 
			(level_target < maxlev));
  
  if(are_same_bin && trg_isnot_leaf && src_isnot_leaf){
    cilk_for(int i=target_children_start; i<target_children_stop; i++){
      for(int j=source_children_start; j<source_children_stop; j++){
	interaction_list_compressed_expanded(clgs_link_list, 
					     clgs_count,
					     nn_link_list, 
					     nn_count,
					     common_list, 
					     common_count,
					     target_tree_nodes, target_tree_edges,
					     source_tree_nodes, source_tree_edges,
					     i, j, level_target+1,
					     level_source+1, maxlev, selected_op);
      }
    }
    
  }

  else if((!src_isnot_leaf || !trg_isnot_leaf) ){
#ifdef NO_SYMBOLIC

    selected_op(nn_count[level_target], 
		&nn_link_list[level_target][target*NN], target, source);
#else
    selected_op(nn_count[level_target], 
		nn_link_list[level_target], target, source);
#endif
  }
  else{

    int offset = target_children_start;
    for(int i=target_children_start; i<target_children_stop; i++){
      
#ifndef MORTON_ONLY
      have_neighbors_T[i-offset] = are_near_neighbors(&target_tree_nodes[level_target+1][i*DIM],
						      &source_tree_nodes[level_source][DIM*source]);
#else      
      
      have_neighbors_T[i-offset] = are_near_neighbors(target_tree_nodes[level_target+1][i],
						      Source_mcode);
      
#endif

    }

    offset = source_children_start;
    int num_nn_S = 0;
    for(int i=source_children_start; i<source_children_stop; i++){
      
#ifndef MORTON_ONLY
      have_neighbors_S[i-offset] = are_near_neighbors(&source_tree_nodes[level_source+1][i*DIM], 
						      &target_tree_nodes[level_target][target*DIM]);
#else 
      
      have_neighbors_S[i-offset] = are_near_neighbors(source_tree_nodes[level_source+1][i], 
						      Target_mcode);
      
#endif
      if(have_neighbors_S[i-offset]){
#ifdef NO_SYMBOLIC
	selected_op(common_count[level_target], 
		    &common_list[level_target][target*CN], target, i);
#else
	selected_op(common_count[level_target], 
		    common_list[level_target], target, i);

#endif
	num_nn_S++;
      }
      
    }
    
    offset = target_children_start;
    for(int i=target_children_start; i<target_children_stop; i++){
      
      if(!have_neighbors_T[i-offset]){

	for(int j=source_children_start, m=0; j<source_children_stop; j++, m++){

	  if(!have_neighbors_S[m]){
	    interaction_list_compressed_expanded(clgs_link_list, 
						 clgs_count,
						 nn_link_list, 
						 nn_count,
						 common_list, 
						 common_count,
						 target_tree_nodes, target_tree_edges,
						 source_tree_nodes, source_tree_edges,
						 i, j, level_target+1,
						 level_source+1, maxlev, selected_op);
	    
	  }

	}

      }
      else{
	for(int j=source_children_start, m=0; j<source_children_stop; j++, m++){
	  if(!have_neighbors_S[m]){
#ifdef NO_SYMBOLIC
	    selected_op(clgs_count[level_target+1], 
			&clgs_link_list[level_target+1][i*FN], i, j);

#else
	    selected_op(clgs_count[level_target+1], 
			clgs_link_list[level_target+1], i, j);
#endif

	  }
	}


      }

      
    }
    
  }


}


void interaction_list_wrapper_compressed(int **clgs_link_list, 
					 uint32_t **clgs_count,
					 int **nn_link_list, 
					 uint32_t **nn_count,
					 int **common_list, 
					 uint32_t **common_count,
					 int **target_tree_nodes, 
					 int **target_tree_edges,
					 int **source_tree_nodes, 
					 int **source_tree_edges,
					 int target, int *nodes_per_level,
					 int level, int maxlev, base_function selected_op){
  for(int i=0; i<nodes_per_level[level]; i++){
    interaction_list_compressed(clgs_link_list, 
				clgs_count, 
				nn_link_list, 
				nn_count,
				common_list, 
				common_count,
				target_tree_nodes, target_tree_edges,
				source_tree_nodes, source_tree_edges,
				target, i, level, 
				level, maxlev, selected_op);
  }

}

void interaction_list_compressed_driver(int **clgs_link_list, 
					uint32_t **clgs_count,
					int **nn_link_list, 
					uint32_t**nn_count,
					int **common_list, 
					uint32_t **common_count,
					int **target_tree_nodes, 
					int **target_tree_edges,
					int **source_tree_nodes, 
					int **source_tree_edges,
					int *nodes_per_level_target, 
					int *nodes_per_level_source, 
					int maxlev, int operation){


  int level = 0;

  base_function use_function = (operation==0) ? &increase_counter : &store_pointer;

  cilk_for(int i=0; i<nodes_per_level_target[level]; i++){
    
    interaction_list_wrapper_compressed(clgs_link_list, 
					clgs_count,
					nn_link_list, 
					nn_count,
					common_list, 
					common_count,
					target_tree_nodes, 
					target_tree_edges,
					source_tree_nodes, 
					source_tree_edges,
					i, nodes_per_level_source,
					level, maxlev-1, use_function);
    
  }
  //writeTracer();
}

#ifndef MORTON_ONLY

 void interaction_list_compressed_expanded_driver(int (**restrict clgs_link_list), 
						  uint32_t (**restrict clgs_count),
						  int (**restrict nn_link_list), 
						  uint32_t (**restrict nn_count),
						  int (**restrict common_list), 
						  uint32_t (**restrict common_count),
						  int **target_tree_nodes, 
						  int **target_tree_edges,
						  int **source_tree_nodes, 
						  int **source_tree_edges,
						  int (*restrict nodes_per_level_target), 
						  int (*restrict nodes_per_level_source), 
						  int maxlev, int operation){
#else
 void interaction_list_compressed_expanded_driver(int (**restrict clgs_link_list), 
						  uint32_t (**restrict clgs_count),
						  int (**restrict nn_link_list), 
						  uint32_t (**restrict nn_count),
						  int (**restrict common_list), 
						  uint32_t (**restrict common_count),
						  uint64_t **target_tree_nodes,
						  int **target_tree_edges,
						  uint64_t **source_tree_nodes,
						  int **source_tree_edges,
						  int (*restrict nodes_per_level_target), 
						  int (*restrict nodes_per_level_source), 
						  int maxlev, int operation){
#endif


  int level = 0;
  
  base_function use_function = (operation==0) ? &increase_counter : &store_pointer;
  
  cilk_for(int i=0; i<nodes_per_level_target[level]; i++){
    for(int j=0; j<nodes_per_level_source[level]; j++){
	interaction_list_compressed_expanded(clgs_link_list, 
					     clgs_count, 
					     nn_link_list, 
					     nn_count,
					     common_list, 
					     common_count,
					     target_tree_nodes, target_tree_edges,
					     source_tree_nodes, source_tree_edges,
					     i, j, level, 
					     level, maxlev-1, use_function);

    }    
    
  }
}
  

void interaction_list_dense_wrapper(int (**restrict clgs_link_list), 
				    int (**restrict nn_link_list), 
				    int (**restrict common_list),
				    int (*restrict nodes_per_level),
				    int maxlev){
  /*
  for(int level=0; level<maxlev; level++){
    cilk_for(int i=0; i<nodes_per_level[level]; i++){
      interaction_list_compressed_dense(clgs_link_list, 
						   nn_link_list, 
						   common_list,
						   i, level, 
						   maxlev);
  
    }
  }
  */

}

int cmpfunc2 (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}


int near_neighbor_roule(int (*restrict candidate),  
			int (*restrict child) ){

  int abs_diff_x = abs(candidate[0] - child[0]);
  int abs_diff_y = abs(candidate[1] - child[1]);
  int abs_diff_z = abs(candidate[2] - child[2]);

  int are_neighbors = (abs_diff_x<=1) && (abs_diff_y<=1) && (abs_diff_z<=1);
  int are_same_bin = (abs_diff_x==0) && (abs_diff_y==0) && (abs_diff_z==0);

  return(are_neighbors && !are_same_bin);
}


int same_bin(int (*restrict candidate),  
	     int (*restrict child) ){

  int abs_diff_x = abs(candidate[0] - child[0]);
  int abs_diff_y = abs(candidate[1] - child[1]);
  int abs_diff_z = abs(candidate[2] - child[2]);

  int are_same_bin = (abs_diff_x==0) && (abs_diff_y==0) && (abs_diff_z==0);

  return(are_same_bin);

} 

void interaction_list_stencil(int *common_stencil, int *far_stencil, 
			      int *near_stencil, int *parent_code){ 


  int mask[DIM] = {0};
  int candidate[DIM] = {0};
  int child[DIM] = {0};
  int common_count = 0;
  int nn_count[8] = {0};
  int far_count[8] = {0};
  int child_count = 0;

  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      for(int k=0; k<6; k++){

	mask[0] = i-2; mask[1] = j-2; mask[2] = k-2;
	candidate[0] = 2*parent_code[0] + mask[0]; 
	candidate[1] = 2*parent_code[1] + mask[1];
	candidate[2] = 2*parent_code[2] + mask[2];

	if(are_near_neighbors( candidate, parent_code )){ // Record far neighbors
	  common_stencil[common_count*DIM + 0] = mask[0];
	  common_stencil[common_count*DIM + 1] = mask[1];
	  common_stencil[common_count*DIM + 2] = mask[2];
	  common_count++;
	}
	else{

	  child_count = 0;
	  for(int m=0; m<2; m++){
	    for(int p=0; p<2; p++){
	      for(int q=0; q<2; q++){

		child[0] = 2*parent_code[0] + q;
		child[1] = 2*parent_code[1] + p;
		child[2] = 2*parent_code[2] + m;

		if(near_neighbor_roule(candidate, child)){ // check if the boxes are near neighbors

		  near_stencil[child_count*NN*DIM + DIM*nn_count[child_count] + 0] = mask[0];
		  near_stencil[child_count*NN*DIM + DIM*nn_count[child_count] + 1] = mask[1];
		  near_stencil[child_count*NN*DIM + DIM*nn_count[child_count] + 2] = mask[2];
		  nn_count[child_count]++;
		}
		else if(!same_bin(candidate, child)){ // check id the boxes are far neighbors
		  far_stencil[child_count*FN*DIM + DIM*far_count[child_count] + 0] = mask[0];
		  far_stencil[child_count*FN*DIM + DIM*far_count[child_count] + 1] = mask[1];
		  far_stencil[child_count*FN*DIM + DIM*far_count[child_count] + 2] = mask[2];
		  far_count[child_count]++;
		}
		child_count++;

	      }
	    }
	  }


	}

      }
    }
  }

}


 int executeStencilNN(int *nn_codes, int *node_code, int *stencil, int periodic, int lv){

   int parent[DIM];
   int bound = (1<<lv);

   int nn_count = 0;

   parent[0] = 2 * (node_code[0] / 2);
   parent[1] = 2 * (node_code[1] / 2);
   parent[2] = 2 * (node_code[2] / 2);

   int lc = (node_code[0] & 1) | ( (node_code[1] & 1) << 1 ) | ( (node_code[2] & 1) << 2 );

   if(periodic){
     for(int i=0; i<NN; i++){
       nn_codes[i*DIM + 0] = parent[0] + stencil[lc*DIM*NN + i*DIM + 0];
       nn_codes[i*DIM + 1] = parent[1] + stencil[lc*DIM*NN + i*DIM + 1];
       nn_codes[i*DIM + 2] = parent[2] + stencil[lc*DIM*NN + i*DIM + 2];
     }

     nn_count = NN;
   }
   else{
    for(int i=0; i<NN; i++){
      nn_codes[nn_count*DIM + 0] = parent[0] + stencil[lc*DIM*NN + i*DIM + 0];
      if(nn_codes[nn_count*DIM + 0] >=0 && nn_codes[nn_count*DIM + 0] < bound){

        nn_codes[nn_count*DIM + 1] = parent[1] + stencil[lc*DIM*NN + i*DIM + 1];
        if(nn_codes[nn_count*DIM + 1] >= 0 && nn_codes[nn_count*DIM + 1] < bound){

          nn_codes[nn_count*DIM + 2] = parent[2] + stencil[lc*DIM*NN + i*DIM + 2];
          if(nn_codes[nn_count*DIM + 2] >= 0 && nn_codes[nn_count*DIM + 2] < bound){

            nn_count++;
          }
        }
      }
    }

  }

  return(nn_count);

 }


 int executeStencilFN(int *fn_codes, int *node_code, int *fn_stencil, int *cn_stencil, int periodic, int lv){

   int parent[DIM];
   int bound = (1 << lv);
   int fn_count = 0;
   int cn_count = 0;

   parent[0] = 2 * (node_code[0] / 2);
   parent[1] = 2 *(node_code[1] / 2);
   parent[2] = 2 *(node_code[2] / 2);

   int lc = (node_code[0] & 1) | ( (node_code[1] & 1) << 1 ) | ( (node_code[2] & 1) << 2 );

   if(periodic){

     for(int i=0; i<CN; i++){
       fn_codes[i*DIM + 0] = parent[0] + cn_stencil[i*DIM + 0];
       fn_codes[i*DIM + 1] = parent[1] + cn_stencil[i*DIM + 1];
       fn_codes[i*DIM + 2] = parent[2] + cn_stencil[i*DIM + 2];
     }

     for(int i=0; i<FN; i++){
       fn_codes[CN*DIM + i*DIM + 0] = parent[0] + fn_stencil[lc*DIM*FN + i*DIM + 0];
       fn_codes[CN*DIM + i*DIM + 1] = parent[1] + fn_stencil[lc*DIM*FN + i*DIM + 1];
       fn_codes[CN*DIM + i*DIM + 2] = parent[2] + fn_stencil[lc*DIM*FN + i*DIM + 2];
     }
     fn_count = FN;
     cn_count = CN;
   }
   else{

     for(int i=0; i<CN; i++){
       fn_codes[cn_count*DIM + 0] = parent[0] + cn_stencil[i*DIM + 0];
       if(fn_codes[cn_count*DIM + 0] >= 0 && fn_codes[cn_count*DIM + 0] < bound){

	 fn_codes[cn_count*DIM + 1] = parent[1] + cn_stencil[i*DIM + 1];
	 if(fn_codes[cn_count*DIM + 1] >= 0 && fn_codes[cn_count*DIM + 1] < bound){

	   fn_codes[cn_count*DIM + 2] = parent[2] + cn_stencil[i*DIM + 2];
	   if(fn_codes[cn_count*DIM + 2] >= 0 && fn_codes[cn_count*DIM + 2] < bound){
	     cn_count++;
	   }
	 }
       }
     }
     for(int i=0; i<FN; i++){
       fn_codes[cn_count*DIM + fn_count*DIM + 0] = parent[0] + fn_stencil[lc*DIM*FN + i*DIM + 0];
       if(fn_codes[cn_count*DIM + fn_count*DIM + 0] >= 0 && fn_codes[cn_count*DIM + fn_count*DIM + 0] < bound){

	 fn_codes[cn_count*DIM + fn_count*DIM + 1] = parent[1] + fn_stencil[lc*DIM*FN + i*DIM + 1];
	 if(fn_codes[cn_count*DIM + fn_count*DIM + 1] >= 0 && fn_codes[cn_count*DIM + fn_count*DIM + 1] < bound){

	   fn_codes[cn_count*DIM + fn_count*DIM + 2] = parent[2] + fn_stencil[lc*DIM*FN + i*DIM + 2];
	   if(fn_codes[cn_count*DIM + fn_count*DIM + 2] >= 0 && fn_codes[cn_count*DIM + fn_count*DIM + 2] < bound){
	     fn_count++;
	   }
	 }
       }
     }

   }

   return(fn_count + cn_count);
 }



#ifndef MORTON_ONLY
void interaction_list_formation(int **node_codes, int **children_first,
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
				int height, int height2, int N, 
				double *memory_count, double *workspace_memory, 
				double *physical_memory, 
				double *interaction_list_physical, 
				double *interaction_list_workspace){
#else
void interaction_list_formation(uint64_t **node_codes, int **children_first,
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
				int height, int height2, int N, 
				double *memory_count, double *workspace_memory,
				double *physical_memory, 
				double *interaction_list_physical, 
				double *interaction_list_workspace){

#endif

  struct timeval startwtime, endwtime;
  const bool printNow = true;
  /* Start the link list formation */

  uint32_t *clgs_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t)); // Array that holds the amount of memory per level required by the colleague linked list 
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);

  uint32_t *nn_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t)); // Array that holds the amount of memory per level required for the neighbors linked list
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);

  uint32_t *common_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t));
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);

  int operation = 0; // First part, symbolic pre-processing    

#ifndef DENSE
#ifndef NO_SYMBOLIC 


 #ifndef INTERLEAVE
   start_timer();
 #endif
   operation = 0; // First part, symbolic pre-processing    

   interaction_list_compressed_expanded_driver(clgs_link_list, 
					       clgs_count,
					       nn_link_list, 
					       nn_count, 
					       common_list, 
					       common_count,
					       node_codes, children_first,
					       node_codes2, children_first2,
					       nodes_per_level, 
					       nodes_per_level2, height, operation);    

 #ifndef INTERLEAVE
   stop_timer("Count neighbors");
 #endif


 #endif 
#endif
   // Find out how much memory we need for the targets linked list

#ifndef DENSE
#ifdef NO_SYMBOLIC

   for(int i=0; i<height; i++){
     clgs_memory_per_level[i] = FN * nodes_per_level[i];
     common_memory_per_level[i] = CN * nodes_per_level[i];
     nn_memory_per_level[i] = NN * nodes_per_level[i];
   }

#else

#ifndef INTERLEAVE
  start_timer();
#endif    
  
  for(int i=0; i<height; i++){
    cilk_spawn scan_colleagues(&clgs_memory_per_level[i], 
			       clgs_count[i], 
			       nodes_per_level[i]);
    cilk_spawn scan_colleagues(&nn_memory_per_level[i],
			       nn_count[i],
			       nodes_per_level[i]);
    cilk_spawn scan_colleagues(&common_memory_per_level[i],
			       common_count[i],
			       nodes_per_level[i]);
    
    
  }
  cilk_sync;
  
#ifndef INTERLEAVE
  stop_timer("Scan colleagues");
#endif
  
  
#endif
#endif

#ifndef DENSE
  // Allocate the required memory for the linked list  
  for(int i=0; i<height; i++){
    clgs_link_list[i] = (int *)sakura_calloc(clgs_memory_per_level[i], sizeof(int), 
					     "Far neighbors interaction list");
    memory_count[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)clgs_memory_per_level[i]*sizeof(int);

#ifdef DENSE
    clgs_link_list[i][0:clgs_memory_per_level[i]] = -1;
#endif
    int pp = (nn_memory_per_level[i]>0)? nn_memory_per_level[i] : 1000;
    nn_link_list[i] = (int *)sakura_calloc(pp, sizeof(int), 
					   "Near neighbors interaction list");
    memory_count[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)clgs_memory_per_level[i]*sizeof(int);

    common_list[i] = (int *)sakura_calloc(common_memory_per_level[i], sizeof(int), 
					  "Common neighbors interaction list");
    memory_count[0] += (double)common_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)common_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)common_memory_per_level[i]*sizeof(int);

#ifdef DENSE
    common_list[i][0:common_memory_per_level[i]] = -1;
#endif
  }
  nn_link_list[height-1][0:nn_memory_per_level[height-1]] = -1;

#endif


#ifdef DENSE
  /*
  interaction_list_compressed_dense(clgs_link_list, 
				    nn_link_list, 
				    common_list,
				    0, height-2, 
				    height);
  */

#ifndef INTERLEAVE
  start_timer();
#endif    

  printf("Generate the stencil\n");

  int toy_parent[DIM];
  toy_parent[0] = 1; toy_parent[1] = 1; toy_parent[2] = 1;
  interaction_list_stencil(common_stencil, far_stencil, near_stencil, toy_parent);

  

#ifndef INTERLEAVE
  stop_timer("Link list");
#endif    

#else


  // Form the link-list
#ifndef INTERLEAVE
  start_timer();
#endif    

  operation = 1; // Second part, formation of the interaction list   

  interaction_list_compressed_expanded_driver(clgs_link_list, 
					      clgs_count,
					      nn_link_list, 
					      nn_count, 
					      common_list, 
					      common_count,
					      node_codes, children_first,
					      node_codes2, children_first2,
					      nodes_per_level, 
					      nodes_per_level2, height, operation);    

#ifndef INTERLEAVE
  stop_timer("Link list");
#endif    
  
#endif

  // free the memory assosiated with the linked list
  free(nn_memory_per_level);
  free(clgs_memory_per_level);

}

