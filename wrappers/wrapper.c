#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void FMM_Init(int images, int thread, double theta, double cutoff, int verbose );
void FMM_Partition( int* ni, int nimax, int* res_index, double* x, double* q, double* v, double* cycle);
void FMM_FMM( int ni, int *nj, int* res_index, double *x, double *q, double *p, double *f, double* cycle);
void FMM_Finalize();

void fmm_init_wrapper_(int *images, int * thread, double *theta, double *cutoff, int *verbose ){
  FMM_Init( *images, *thread, *theta, *cutoff, *verbose);
}
void fmm_partition_wrapper_( 
    int* ni, int *nimax, int* res_index, 
    double* x, double* q, double* v, double* cycle) {
  int nimax2 = *nimax ;
  FMM_Partition(ni, nimax2, res_index, x, q, v, cycle);
}

void fmm_fmm_wrapper_(
    int *ni, int *nj, int* res_index, 
    double *x, double *q, double *p, double *f, double* cycle){
  int ni2 = *ni;
  FMM_FMM(ni2, nj, res_index, x, q, p, f, cycle); 
}

void fmm_finalize_(){
  FMM_Finalize();
}
