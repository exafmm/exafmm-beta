#include <cilk/cilk.h>
#include <math.h>
#include "utils.h"
#include "core.h"

void upward_pass(float *X2, double ***Multipole, int **node_codes2,
		 int** c_count2, int** node_pointers2, int *leaf_populations2,
		 float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count2[level][node_id-1];
  int c_end = c_count2[level][node_id];
  int c_size = c_end - c_begin;
  double ranges[DIM], dX[DIM];
  ranges[:] = fabs(Xmax[0:DIM] - Xmin[0:DIM]);
  ranges[:] *= 1.00001;
  double range = __sec_reduce_max(ranges[:]);
  if(c_size==0){ // P2M
    int l_begin = node_pointers2[level][node_id];
    int l_size = leaf_populations2[l_begin];
    int l_end = l_begin + l_size;
    int nbins = 1 << (level + 1);
    double qstep = range / nbins;
    double Xnode[DIM];
    for(int d=0; d<DIM; d++){
      Xnode[d] = qstep * (node_codes2[level][3*node_id+d] + .5) + Xmin[d];
    }
    for(int i=l_begin; i<l_end; i++){
      for(int d=0; d<DIM; d++) dX[d] = Xnode[d] - X2[LDIM*i+d];
      double M[MTERM];
      M[0] = X2[LDIM*i+3];
      powerM(M,dX);
      for(int m=0; m<MTERM; m++) Multipole[level][node_id][m] += M[m];
    }
  }else{ // M2M
    int nbins = 1 << (level + 1);
    double step_p = range / nbins;
    double step_c = range / nbins / 2;
    double Xp[DIM], Xc[DIM];
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
      double M[MTERM], C[LTERM];
      C[0] = 1;
      powerM(C,dX);
      for(int m=0; m<MTERM; m++){
	M[m] = Multipole[level+1][i][m];
	Multipole[level][node_id][m] += C[m] * M[0];
      }
      M2MSum(Multipole[level][node_id],C,M);
    }
  }
}

void evaluation(float *X, float *X2, float *TRG, double ***Multipole, double ***Local,
		int *nodes_per_level, int **node_pointers, int **node_codes, int *leaf_populations,
		int **node_pointers2, int **node_codes2, int *leaf_populations2,
		int **c_count, int **n_list, uint32_t **n_count,
		int **f_list, uint32_t **f_count, int **s_list, uint32_t **s_count,
		float *Xmin, float *Xmax, int height){
  double Xtarget[DIM], Xsource[DIM], ranges[DIM], dX[DIM];
  ranges[:] = fabs(Xmax[0:DIM] - Xmin[0:DIM]);
  ranges[:] *= 1.00001;
  double range = __sec_reduce_max(ranges[:]);
  for(int level=0; level<height; level++){
    int nbins = 1 << (level + 1);
    double qstep = range / nbins;
    cilk_for(int node_id=0; node_id<nodes_per_level[level]; node_id++){
      int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
      int c_end = c_count[level][node_id];
      int n_begin = (node_id==0)? 0 : n_count[level][node_id-1];
      int n_end = n_count[level][node_id];
      int f_begin = (node_id==0) ? 0 : f_count[level][node_id-1];
      int f_end = f_count[level][node_id];
      int s_begin = (node_id==0) ? 0 : s_count[level][node_id-1];
      int s_end = s_count[level][node_id];
      int i_begin = node_pointers[level][node_id];
      int i_size = leaf_populations[i_begin];
      int i_end = i_begin + i_size;
      for(int n=n_begin; n<n_end; n++){ // P2P
	int j_begin = node_pointers2[level][n_list[level][n]];
	int j_size = leaf_populations2[j_begin];
	int j_end = j_begin + j_size;
	for(int i=i_begin; i<i_end; i++) {
	  for(int j=j_begin; j<j_end; j++) {
	    for(int d=0; d<DIM; d++) dX[d] = X[LDIM*i+d] - X2[LDIM*j+d];
	    double R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	    double invR2 = 1.0 / R2;
	    if( R2 == 0 ) invR2 = 0;
	    double invR = X2[LDIM*j+3] * sqrt(invR2);
	    double invR3 = invR2 * invR;
	    TRG[4*i+0] += invR;
	    TRG[4*i+1] -= dX[0] * invR3;
	    TRG[4*i+2] -= dX[1] * invR3;
	    TRG[4*i+3] -= dX[2] * invR3;
#ifdef TEST
	    if(i==i_begin) Local[level][node_id][0]++;
#endif
	  }
	}
      }
      double L[LTERM];
      for(int l=0; l<LTERM; l++) L[l] = 0;
      for(int i=f_begin; i<f_end; i++){ // M2L
	for(int d=0; d<DIM; d++){
	  Xtarget[d] = qstep * (node_codes[level][3*node_id+d] + .5) + Xmin[d];
	  Xsource[d] = qstep * (node_codes2[level][3*f_list[level][i]+d] + .5) + Xmin[d];
	  dX[d] = Xtarget[d] - Xsource[d];
	}
	double invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
	double invR  = sqrt(invR2);
	double C[LTERM], M[MTERM];
	getCoef(C,dX,invR2,invR);
	for(int m=0; m<MTERM; m++){
	  M[m] = Multipole[level][f_list[level][i]][m];
#ifdef TEST
	  if(m>0) M[m] = 0;
#endif
	}
#ifdef TEST
	C[0] = 1;
#endif
	M2LSum(L,C,M);
      }
      for(int l=0; l<LTERM; l++) Local[level][node_id][l] += L[l];
      for(int i=c_begin; i<c_end; i++){
	for(int l=0; l<LTERM; l++) L[l] = 0;
	for(int j=s_begin; j<s_end; j++){ // M2L
	  for(int d=0; d<DIM; d++){
	    Xtarget[d] = .5 * qstep * (node_codes[level+1][3*i+d] + .5) + Xmin[d];
	    Xsource[d] = .5 * qstep * (node_codes2[level+1][3*s_list[level][j]+d] + .5) + Xmin[d];
	    dX[d] = Xtarget[d] - Xsource[d];
	  }
	  double invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
	  double invR  = sqrt(invR2);
	  double C[LTERM], M[MTERM];
	  getCoef(C,dX,invR2,invR);
	  for(int m=0; m<MTERM; m++){
	    M[m] = Multipole[level+1][s_list[level][j]][m];
#ifdef TEST
	    if(m>0) M[m] = 0;
#endif
	  }
#ifdef TEST
	  C[0] = 1;
#endif
	  M2LSum(L,C,M);
	}
	for(int l=0; l<LTERM; l++) Local[level+1][i][l] += L[l];
      }
    }
  }
}

void downward_pass(float *X, float *TRG, double ***Local, int **node_codes,
		   int** c_count, int** node_pointers, int *leaf_populations,
		   float *Xmin, float *Xmax, int node_id, int level){
  int c_begin = (node_id==0) ? 0 : c_count[level][node_id-1];
  int c_end = c_count[level][node_id];
  int c_size = c_end - c_begin;
  double ranges[DIM], dX[DIM];
  ranges[:] = fabs(Xmax[0:DIM] - Xmin[0:DIM]);
  ranges[:] *= 1.00001;
  double range = __sec_reduce_max(ranges[:]);
  if(c_size==0){ // L2P
    int l_begin = node_pointers[level][node_id];
    int l_size = leaf_populations[l_begin];
    int l_end = l_begin + l_size;
    int nbins = 1 << (level + 1);
    double qstep = range / nbins;
    double Xnode[DIM];
    for(int d=0; d<DIM; d++){
      Xnode[d] = qstep * (node_codes[level][3*node_id+d] + .5) + Xmin[d];
    }
    double L[LTERM];
    for(int l=0; l<LTERM; l++) L[l] = Local[level][node_id][l];
    for(int i=l_begin; i<l_end; i++){
      for(int d=0; d<DIM; d++) dX[d] = X[LDIM*i+d] - Xnode[d];
      double C[LTERM];
      C[0] = 1;
      powerL(C,dX);
      for(int d=0; d<4; d++) TRG[4*i+d] += L[d];
      for(int l=1; l<LTERM; l++) TRG[4*i+0] += C[l] * L[l];
      L2PSum(&TRG[4*i],C,L);
    }
  }else{ // L2L
    int nbins = 1 << (level + 1);
    double step_p = range / nbins;
    double step_c = range / nbins / 2;
    double Xp[DIM], Xc[DIM];
    for(int d=0; d<DIM; d++){
      Xp[d] = step_p * (node_codes[level][3*node_id+d] + .5) + Xmin[d];
    }
    for(int i=c_begin; i<c_end; i++){
      for(int d=0; d<DIM; d++){
	Xc[d] = step_c * (node_codes[level+1][3*i+d] + .5) + Xmin[d];
	dX[d] = Xp[d] - Xc[d];
      }
      double C[LTERM];
      C[0] = 1;
      powerL(C,dX);
      for(int l=0; l<LTERM; l++) Local[level+1][i][l] += Local[level][node_id][l];
#ifndef TEST
      for(int l=1; l<LTERM; l++) Local[level+1][i][0] += C[l] * Local[level][node_id][l];
#endif
      L2LSum(Local[level+1][i],C,Local[level][node_id]);      
      downward_pass(X, TRG, Local, node_codes, c_count, node_pointers, leaf_populations,
		    Xmin, Xmax, i, level+1);
    }
  }
}
