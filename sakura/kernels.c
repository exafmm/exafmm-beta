void upward_pass(int **expansions, int** c_count2,
		int** node_pointers2, int *leaf_populations,
		int node_id, int level){
  int children_stop = c_count2[level][node_id];
  int children_start = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int isnotleaf = children_stop - children_start;
  if(isnotleaf==0){ // P2M
    expansions[level][node_id] = leaf_populations[node_pointers2[level][node_id]];
  }else{ // M2M
    for(int i=children_start; i<children_stop; i++){
      upward_pass(expansions, c_count2, node_pointers2, leaf_populations, i, level+1);
      expansions[level][node_id] += expansions[level+1][i];
    }
  }
}
