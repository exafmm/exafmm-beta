#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <math.h>
#include "utils.h"

int main(int argc, char** argv){
  int N = atoi(argv[1]);
  int dist = atoi(argv[2]);
  int population_threshold = atoi(argv[3]);
  int nworkers = atoi(argv[4]);
  __cilkrts_set_param("nworkers",argv[4]);
  printf("N = %d, T=%d\n", N, nworkers);
  start_timer();
  float *X = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle array");
  float *X2 = (float *)sakura_malloc(N, LDIM*sizeof(float), "Particle array");
  float *TRG = (float *)sakura_malloc(N, 4*sizeof(float), "Result array");
  float *TRG2 = (float *)sakura_malloc(N, 4*sizeof(float), "Result array");
  stop_timer("Data mem. alloc.");
  int maxlev = 20;
  start_timer();
#if 0
  create_dataset(X, N, dist);
  create_dataset(X2, N, dist);
#else
  for(int i=0; i<N; i++){
    for(int d=0; d<3; d++){
      X[LDIM*i+d] = X2[LDIM*i+d] = i;
    }
    for(int d=0; d<4; d++){
      TRG[4*i+d] = TRG2[4*i+d] = 0;
    }
  }
  X[LDIM*0+3] = X2[LDIM*0+3] = 1;
  X[LDIM*1+3] = X2[LDIM*1+3] = 0;
  X[LDIM*2+3] = X2[LDIM*2+3] = 0;
  X[LDIM*3+3] = X2[LDIM*3+3] = 1;
#endif
  stop_timer("Create data");
  start_timer();
  float Xmin[DIM], Xmax[DIM];
  for(int i=0; i<DIM; i++){
    Xmin[i] = __sec_reduce_min(X[i:N:LDIM]);
    Xmax[i] = __sec_reduce_max(X[i:N:LDIM]);
  }
  float Xmin2[DIM], Xmax2[DIM];
  for(int i=0; i<DIM; i++){
    Xmin2[i] = __sec_reduce_min(X2[i:N:LDIM]);
    Xmax2[i] = __sec_reduce_max(X2[i:N:LDIM]);
  }
  Xmin[:] = MIN(Xmin[:], Xmin2[:]);
  Xmax[:] = MAX(Xmax[:], Xmax2[:]);
  stop_timer("Box bounds");
  uint64_t *particle_codes = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Morton code array");
  uint32_t *bit_map = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map");
  uint32_t *permutation_vector = (uint32_t *)sakura_malloc(N, sizeof(uint32_t),"Permutation vector");
  uint64_t *particle_codes2 = (uint64_t *)sakura_malloc(N, sizeof(uint64_t), "Morton code array");
  uint32_t *bit_map2 = (uint32_t *)sakura_calloc(N, sizeof(uint32_t), "Bit map");
  uint32_t *permutation_vector2 = (uint32_t *)sakura_malloc(N, sizeof(uint32_t),"Permutation vector");
  encodeParticles(N, X, Xmin, Xmax, particle_codes, maxlev);
  decomposeSpace(N, &particle_codes, permutation_vector, bit_map, maxlev, population_threshold);
  relocateParticles(N, &X, permutation_vector);
  encodeParticles(N, X2, Xmin, Xmax, particle_codes2, maxlev);
  decomposeSpace(N, &particle_codes2, permutation_vector2, bit_map2, maxlev, population_threshold);
  relocateParticles(N, &X2, permutation_vector2);

  int nodes_per_level[20];
  int **node_pointers = (int **)malloc(maxlev*sizeof(int *)); 
  int **num_children = (int **)malloc(maxlev*sizeof(int *)); 
  int **c_count = (int **)malloc(maxlev*sizeof(int *)); 
  int **node_codes = (int **)malloc(maxlev*sizeof(int *));
  int nodes_per_level2[20];
  int **node_pointers2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **num_children2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **c_count2 = (int **)malloc(maxlev*sizeof(int *)); 
  int **node_codes2 = (int **)malloc(maxlev*sizeof(int *));
  int height = tree_formation(bit_map, particle_codes, 
			      nodes_per_level, node_pointers, 
			      num_children, c_count, 
			      node_codes, maxlev, N);
  int height2 = tree_formation(bit_map2, particle_codes2, 
			       nodes_per_level2, node_pointers2, 
			       num_children2, c_count2, 
			       node_codes2, maxlev, N);
  int **n_list = (int **)malloc(height*sizeof(int *)); 
  int **f_list = (int **)malloc(height*sizeof(int *)); 
  int **s_list = (int **)malloc(height*sizeof(int *));
  uint32_t **n_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
  uint32_t **f_count = (uint32_t **)malloc(height*sizeof(uint32_t *)); 
  uint32_t **s_count = (uint32_t **)malloc(height*sizeof(uint32_t *));
  form_interaction_lists(node_codes, c_count, node_codes2, c_count2, 
			 n_count, f_count, s_count, n_list, f_list, s_list,
			 nodes_per_level, nodes_per_level2, height);
  double ***Multipole = (double ***)malloc(height2*sizeof(double**));
  double ***Local = (double ***)malloc(height*sizeof(double**));
  for(int l=0; l<height2; l++){
    Multipole[l] = (double **)malloc(nodes_per_level2[l]*sizeof(double*));
    for (int i=0; i<nodes_per_level2[l]; i++) {
      Multipole[l][i] = (double *)sakura_calloc(MTERM,sizeof(double),"Multipole Multipole");
    }
  }
  for(int l=0; l<height; l++){
    Local[l] = (double **)malloc(nodes_per_level[l]*sizeof(double*));
    for (int i=0; i<nodes_per_level[l]; i++) {
      Local[l][i] = (double *)sakura_calloc(LTERM,sizeof(double),"Local Multipole");
    }
  }
  int *leaf_populations = (int *)sakura_calloc(N, sizeof(int),"Leaf population array");
  uint64_t numleaves = find_leaf_populations(leaf_populations, bit_map, N);
  int *leaf_populations2 = (int *)sakura_calloc(N, sizeof(int),"Leaf population array");
  uint64_t numleaves2 = find_leaf_populations(leaf_populations2, bit_map2, N);
  int charge = 0;
  start_timer();
  for(int i=0; i<nodes_per_level2[0]; i++){
    upward_pass(X2, Multipole, node_codes2, c_count2, node_pointers2, leaf_populations2,
		Xmin, Xmax, i, 0);
    charge += Multipole[0][i][0];
  }
#ifdef TEST
  printf("Tree %s\n", (charge == N ? "PASS" : "FAIL"));
#endif
  stop_timer("Upward pass");
  start_timer();
  evaluation(X, X2, TRG, Multipole, Local, nodes_per_level, node_pointers, node_codes, leaf_populations,
	     node_pointers2, node_codes2, leaf_populations2,
	     c_count, n_list, n_count, f_list, f_count, s_list, s_count,
	     Xmin, Xmax, height);
  stop_timer("Evaluation");
  start_timer();
  for(int i=0; i<nodes_per_level[0]; i++){
    downward_pass(X, TRG, Local, node_codes, c_count, node_pointers, leaf_populations,
		  Xmin, Xmax, i, 0);
  }
  stop_timer("Downward pass");
  int node_id = nodes_per_level[height-1] - 1;
#ifdef TEST
  printf("List %s\n", ((Local[height-1][node_id][0] - N) < 0.1 ? "PASS" : "FAIL"));
#endif
  uint64_t inter_list_edges = 0;
  uint64_t num_tree_nodes = 0;
  for(int i=0;i<height; i++){
    inter_list_edges += n_count[i][nodes_per_level[i]-1] +
      f_count[i][nodes_per_level[i]-1];
    num_tree_nodes += nodes_per_level[i];
  }
  printf("%-20s: %d\n", "Tree height", height);
  printf("%-20s: %lu\n", "Tree nodes", num_tree_nodes);
  printf("%-20s: %lu\n", "Tree leaves", numleaves);
  printf("%-20s: %lu\n", "Edges",inter_list_edges);

#ifndef TEST
  double dX[DIM];
  for(int i=0; i<4; i++){
    for(int j=0; j<N; j++){
      for(int d=0; d<3; d++) dX[d] = X[LDIM*i+d] - X2[LDIM*j+d];
      double R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
      double invR2 = 1.0 / R2;
      if( R2 == 0 ) invR2 = 0;
      double invR = X2[LDIM*j+3] * sqrt(invR2);
      double invR3 = invR2 * invR;
      TRG2[4*i+0] += invR;
      TRG2[4*i+1] -= dX[0] * invR3;
      TRG2[4*i+2] -= dX[1] * invR3;
      TRG2[4*i+3] -= dX[2] * invR3;
    }
    printf("%d %f %f\n",i,TRG[4*i],TRG2[4*i]);
  }
#endif
  
  for(int i=0; i<height; i++){
    free(n_count[i]);
    free(f_count[i]);
    free(s_count[i]);
    free(n_list[i]);
    free(f_list[i]);
    free(s_list[i]);
  }

  free(n_count);
  free(f_count);
  free(s_count);
  free(n_list);
  free(f_list);
  free(s_list);
  for(int i=0; i<height; i++){
    free(node_pointers[i]);
    free(num_children[i]);
    free(c_count[i]);
    free(node_codes[i]);
    free(Local[i]);
  }
  for(int i=0; i<height2; i++){
    free(node_pointers2[i]);
    free(num_children2[i]);
    free(c_count2[i]);
    free(node_codes2[i]);
    free(Multipole[i]);
  }
  free(node_pointers);
  free(num_children);
  free(c_count);
  free(node_codes);
  free(node_pointers2);
  free(num_children2);
  free(c_count2);
  free(node_codes2);
  free(Multipole);
  free(Local);
  free(leaf_populations);
  free(bit_map);
  free(particle_codes);
  free(permutation_vector);
  free(leaf_populations2);
  free(bit_map2);
  free(particle_codes2);
  free(permutation_vector2);
  free(X);
  free(X2);
  free(TRG);
  free(TRG2);
}
