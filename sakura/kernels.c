#include <cilk/cilk.h>
#include <math.h>
#include "utils.h"
#include "core.h"

void upward_pass(float *X2, float (**Multipole)[MTERM], int **node_codes2,
		 int** c_count2, int** node_pointers2, int *leaf_populations2,
		 float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int c_end = c_count2[level][node_id];
  int c_size = c_end - c_begin;
  float ranges[DIM], dX[DIM];
  ranges[:] = fabs(Xmax[0:DIM] - Xmin[0:DIM]);
  ranges[:] *= 1.00001;
  float range = __sec_reduce_max(ranges[:]);
  if(c_size==0){ // P2M
    int l_begin = node_pointers2[level][node_id];
    int l_size = leaf_populations2[l_begin];
    int l_end = l_begin + l_size;
    int nbins = 1 << (level + 1);
    float qstep = range / nbins;
    float Xnode[DIM];
    for(int d=0; d<DIM; d++){
      Xnode[d] = qstep * (node_codes2[level][3*node_id+d] + .5) + Xmin[d];
    }
    for(int i=l_begin; i<l_end; i++){
      for(int d=0; d<DIM; d++) dX[d] = Xnode[d] - X2[LDIM*i+d];
      float M[MTERM];
      M[0] = X2[LDIM*i+3];
      powerM(M,dX);
      for(int m=0; m<MTERM; m++) Multipole[level][node_id][m] += M[m];
    }
  }else{ // M2M
    int nbins = 1 << (level + 1);
    float step_p = range / nbins;
    float step_c = range / nbins / 2;
    float Xp[DIM], Xc[DIM];
    for(int d=0; d<DIM; d++){
      Xp[d] = step_p * (node_codes2[level][3*node_id+d] + .5) + Xmin[d];
    }
    for(int i=c_begin; i<c_end; i++){
      upward_pass(X2, Multipole, node_codes2, c_count2, node_pointers2, leaf_populations2,
		  Xmin, Xmax, i, level+1);
      for(int d=0; d<DIM; d++){
	Xc[d] = step_c * (node_codes2[level+1][3*i+d] + .5) + Xmin[d];
	dX[d] = Xp[d] - Xc[d];
      }
      float M[MTERM], C[LTERM];
      C[0] = 1;
      powerM(C,dX);
      for(int m=0; m<MTERM; m++) {
	M[m] = Multipole[level+1][i][m];
	Multipole[level][node_id][m] += C[m] * M[0];
      }
      M2MSum(Multipole[level][node_id],C,M);
    }
  }
}

void evaluation(float (**Multipole)[MTERM], float (**Local)[LTERM],
		int *nodes_per_level, int **node_codes, int **node_codes2,
		int **c_count, int **n_list, uint32_t **n_count,
		int **f_list, uint32_t **f_count, int **s_list, uint32_t **s_count, int height){
  float Xtarget[DIM], Xsource[DIM], dX[3];
  for(int level=0; level<height; level++){
    for(int node_id=0; node_id<nodes_per_level[level]; node_id++){
      int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
      int c_end = c_count[level][node_id];
      int n_begin = (node_id==0)? 0 : n_count[level][node_id-1];
      int n_end = n_count[level][node_id];
      int f_begin = (node_id==0) ? 0 : f_count[level][node_id-1];
      int f_end = f_count[level][node_id];
      int s_begin = (node_id==0) ? 0 : s_count[level][node_id-1];
      int s_end = s_count[level][node_id];
      for(int i=n_begin; i<n_end; i++){ // P2P
	Local[level][node_id][0] += Multipole[level][n_list[level][i]][0];
      }
      for(int i=f_begin; i<f_end; i++){ // M2L
	Local[level][node_id][0] += Multipole[level][f_list[level][i]][0];
      }
      for(int i=c_begin; i<c_end; i++){
	for(int j=s_begin; j<s_end; j++){ // M2L
	  Local[level+1][i][0] += Multipole[level+1][s_list[level][j]][0];
	}
      }
    }
  }
}

void downward_pass(float *X, float (**Local)[LTERM], int **node_codes,
		   int** c_count, int** node_pointers, int *leaf_populations,
		   float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
  int c_end = c_count[level][node_id];
  int c_size = c_end - c_begin;
  if(c_size==0){ // L2P
  }else{ // L2L
    for(int i=c_begin; i<c_end; i++){
      Local[level+1][i][0] += Local[level][node_id][0];
      downward_pass(X, Local, node_codes, c_count, node_pointers, leaf_populations,
		    Xmin, Xmax, i, level+1);
    }
  }
}
