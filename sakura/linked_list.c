#include "cilk/cilk.h"
#include "utils.h"

typedef void (*base_function)(uint32_t*, int*, int, int); 

void increase_counter(uint32_t (*restrict count), int (*restrict link_list), int target, int source){
  count[target]++;	
}

void store_pointer(uint32_t (*restrict count), int (*restrict link_list), int target, int source){
  uint32_t cursor = count[target];
  link_list[cursor] = source;
  cursor++;
  count[target] = cursor;
}

void scan_colleagues(uint32_t (*restrict clgs_memory), uint32_t (*restrict clgs_first), int N){
  int offset = 0;
  for(int i=0; i<N; i++){
    int ss = clgs_first[i];
    clgs_first[i] = offset;
    offset += ss;
  }
  clgs_memory[0] = offset;
}

int are_near_neighbors( int (*restrict child),  
			int (*restrict parent) ){
  int diff_x = 2*parent[0] - child[0];
  int diff_y = 2*parent[1] - child[1];
  int diff_z = 2*parent[2] - child[2];
  int is_common_neighbor = (diff_x == 2 || diff_x == -3) ||
    (diff_y == 2 || diff_y == -3) ||
    (diff_z == 2 || diff_z == -3);
  return(is_common_neighbor);
}

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
  int have_neighbors_T[8] = {0};
  int have_neighbors_S[8] = {0};
  int target_children_stop = target_tree_edges[level_target][target];
  int target_children_start = (target==0) ? 0 : target_tree_edges[level_target][target-1];
  int source_children_stop = source_tree_edges[level_source][source];
  int source_children_start = (source == 0) ? 0 : source_tree_edges[level_source][source-1];
  int nchild_source = source_children_stop - source_children_start;
  int nchild_target = target_children_stop - target_children_start;
  int Target_code_X = target_tree_nodes[level_target][DIM*target];
  int Target_code_Y = target_tree_nodes[level_target][DIM*target+1];
  int Target_code_Z = target_tree_nodes[level_target][DIM*target+2];
  int Source_code_X = source_tree_nodes[level_source][DIM*source];
  int Source_code_Y = source_tree_nodes[level_source][DIM*source+1];
  int Source_code_Z = source_tree_nodes[level_source][DIM*source+2];
  int abs_diff_x = abs(Target_code_X - Source_code_X);
  int abs_diff_y = abs(Target_code_Y - Source_code_Y);
  int abs_diff_z = abs(Target_code_Z - Source_code_Z);
  int are_same_bin = (abs_diff_x == 0) && (abs_diff_y == 0) && (abs_diff_z == 0);
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
  }else if((!src_isnot_leaf || !trg_isnot_leaf) ){
    selected_op(nn_count[level_target], 
		nn_link_list[level_target], target, source);
  }else{
    int offset = target_children_start;
    for(int i=target_children_start; i<target_children_stop; i++){
      have_neighbors_T[i-offset] = are_near_neighbors(&target_tree_nodes[level_target+1][i*DIM],
						      &source_tree_nodes[level_source][DIM*source]);
    }
    offset = source_children_start;
    int num_nn_S = 0;
    for(int i=source_children_start; i<source_children_stop; i++){
      have_neighbors_S[i-offset] = are_near_neighbors(&source_tree_nodes[level_source+1][i*DIM], 
						      &target_tree_nodes[level_target][target*DIM]);
      if(have_neighbors_S[i-offset]){
	selected_op(common_count[level_target], 
		    common_list[level_target], target, i);
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
	    selected_op(clgs_count[level_target+1], 
			clgs_link_list[level_target+1], i, j);
	  }
	}
      }
    }
  }
}

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
  struct timeval startwtime, endwtime;
  const bool printNow = true;
  uint32_t *clgs_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t));
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);
  uint32_t *nn_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t));
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);
  uint32_t *common_memory_per_level = (uint32_t *)malloc(height*sizeof(uint32_t));
  memory_count[0] += (double)height*sizeof(uint32_t);
  workspace_memory[0] += (double)height*sizeof(uint32_t);
  interaction_list_workspace[0] += (double)height*sizeof(uint32_t);
  int operation = 0;
  start_timer();
  operation = 0;
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
  stop_timer("Count neighbors");
  start_timer();
  for(int i=0; i<height; i++){
    cilk_spawn scan_colleagues(&clgs_memory_per_level[i], clgs_count[i], nodes_per_level[i]);
    cilk_spawn scan_colleagues(&nn_memory_per_level[i], nn_count[i], nodes_per_level[i]);
    cilk_spawn scan_colleagues(&common_memory_per_level[i], common_count[i], nodes_per_level[i]);
  }
  cilk_sync;
  stop_timer("Scan colleagues");
  for(int i=0; i<height; i++){
    clgs_link_list[i] = (int *)sakura_calloc(clgs_memory_per_level[i], sizeof(int), "Far neighbors interaction list");
    memory_count[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    int pp = (nn_memory_per_level[i]>0)? nn_memory_per_level[i] : 1000;
    nn_link_list[i] = (int *)sakura_calloc(pp, sizeof(int), "Near neighbors interaction list");
    memory_count[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)clgs_memory_per_level[i]*sizeof(int);
    common_list[i] = (int *)sakura_calloc(common_memory_per_level[i], sizeof(int), "Common neighbors interaction list");
    memory_count[0] += (double)common_memory_per_level[i]*sizeof(int);
    physical_memory[0] += (double)common_memory_per_level[i]*sizeof(int);
    interaction_list_physical[0] += (double)common_memory_per_level[i]*sizeof(int);
  }
  nn_link_list[height-1][0:nn_memory_per_level[height-1]] = -1;

  start_timer();
  operation = 1;
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
  stop_timer("Link list");
  free(nn_memory_per_level);
  free(clgs_memory_per_level);
}
