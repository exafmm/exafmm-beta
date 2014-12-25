void upward_pass(int **expansions, int** c_count2, int** node_pointers2,
		 int *leaf_populations2, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int c_end = c_count2[level][node_id];
  int c_size = c_end - c_begin;
  if(c_size==0){ // P2M
    expansions[level][node_id] = leaf_populations2[node_pointers2[level][node_id]];
  }else{ // M2M
    for(int i=c_begin; i<c_end; i++){
      upward_pass(expansions, c_count2, node_pointers2, leaf_populations2, i, level+1);
      expansions[level][node_id] += expansions[level+1][i];
    }
  }
}

void downward_pass(int **expansions, int** c_count, int** node_pointers,
		 int *leaf_populations, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
  int c_end = c_count[level][node_id];
  int c_size = c_end - c_begin;
  if(c_size==0){ // L2P
  }else{ // L2L
    for(int i=c_begin; i<c_end; i++){
      expansions[level][node_id] += expansions[level+1][i];
      downward_pass(expansions, c_count, node_pointers, leaf_populations, i, level+1);
    }
  }
}
