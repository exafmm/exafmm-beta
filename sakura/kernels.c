#include "utils.h"

void upward_pass(float *X2, int **expansions, int **node_codes2, int** c_count2, int** node_pointers2,
		 int *leaf_populations2, float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int c_end = c_count2[level][node_id];
  int c_size = c_end - c_begin;
  if(c_size==0){ // P2M
    expansions[level][node_id] = leaf_populations2[node_pointers2[level][node_id]];
    int l_begin = node_pointers2[level][node_id];
    int l_size = leaf_populations2[l_begin];
    int l_end = l_begin + l_size;
    int nbins = 1 << (level + 1);
    float Xnode[3], dX[3];
    for(int d=0; d<3; d++){
      float node_dim = (Xmax[d] - Xmin[d]) / nbins;
      Xnode[d] = node_dim * (node_codes2[level][3*node_id+d] + .5) + Xmin[d];
    }
    for(int i=l_begin; i<l_end; i++){
      for(int d=0; d<3; d++) {
	dX[d] = X2[LDIM*i+d] - Xnode[d];
	float radius = (Xmax[d] - Xmin[d]) / nbins / 2;
	if(dX[d]>radius) printf("%d %d %f: %f\n",i,d,X2[LDIM*i+d],Xnode[d]);
      }
    }
  }else{ // M2M
    for(int i=c_begin; i<c_end; i++){
      upward_pass(X2, expansions, node_codes2, c_count2, node_pointers2, leaf_populations2,
		  Xmin, Xmax, i, level+1);
      expansions[level][node_id] += expansions[level+1][i];
    }
  }
}

void downward_pass(float *X, int **interactions, int **node_codes, int** c_count, int** node_pointers,
		   int *leaf_populations, float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
  int c_end = c_count[level][node_id];
  int c_size = c_end - c_begin;
  if(c_size==0){ // L2P
  }else{ // L2L
    for(int i=c_begin; i<c_end; i++){
      interactions[level+1][i] += interactions[level][node_id];
      downward_pass(X, interactions, node_codes, c_count, node_pointers, leaf_populations,
		    Xmin, Xmax, i, level+1);
    }
  }
}
