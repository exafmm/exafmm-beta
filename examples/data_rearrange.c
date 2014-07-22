/*

This file contain all the functions used for data rearrangement based on the octree indexing

author: Nikos Sismanis
date: Jul 2014  

*/

#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"

#define DIM 3

#ifndef LDIM
#define LDIM 12
#endif

#define NP 128

void relocate(float *Y, float *X, uint *index, int N){

  Y[0:N] = X[index[0:N]];
}

void relocate2(float *Y, float *X, uint *index, int N){

  
  for(int i=0; i<N; i++){
    Y[i] = X[index[i]];
  }
  
}

void relocateT(float *Y, float *X, uint *index, int N){

  for(int i=0; i<N; i++){
    //Y[i*DIM] = X[index[i]*DIM];
    //Y[i*DIM] = X[index[i]*DIM + 1];
    //Y[i*DIM] = X[index[i]*DIM + 2];
    Y[i*DIM:DIM] = X[index[i]*DIM:DIM];
  }

}

void relocateTL(float *Y, float *X, uint *index, int N){

  for(int i=0; i<N; i++){
    //Y[i*LDIM] = X[index[i]*LDIM];
    //Y[i*LDIM] = X[index[i]*LDIM + 1];
    //Y[i*LDIM] = X[index[i]*LDIM + 2];
    Y[i*LDIM:LDIM] = X[index[i]*LDIM:LDIM];
    //printf("%d: %d\n",i, i*LDIM);
  }

}


void relocate_par(float *Y, float *X, uint *index, int N){

  int M = N / NP;

  for(int j=0; j<NP; j++){
    cilk_spawn relocate2(Y, X, &index[j*M], M);
  }
  cilk_sync;
}


void rearrange_data(float *Y, float *X, uint *Ids, int N){


  cilk_for(int i=0; i<DIM; i++){
    relocate_par(&Y[i*N], &X[i*N], Ids, N);
  }

}


void rearrange_dataT(float *Y, float *X, uint *Ids, int N){

  int M = N / NP;

  for(int i=0; i<NP; i++){
    cilk_spawn relocateT(&Y[i*M*DIM], X, &Ids[i*M], M);
  }

}

void rearrange_dataT_nlev(float *Y, float *X, uint *Ids, int N, int nbins, uint *bsizes, uint *baccum_sizes){


  for(int i=0; i<nbins; i++){
    rearrange_dataT(&Y[DIM*baccum_sizes[i]], &X[DIM*baccum_sizes[i]], &Ids[DIM*baccum_sizes[i]], bsizes[i]);
  }

}

void rearrange_dataTL(float *Y, float *X, uint *Ids, int N){

  int M = N / NP;

  for(int i=0; i<NP-1; i++){
    cilk_spawn relocateTL(&Y[i*M*LDIM], X, &Ids[i*M], M);
  }
  relocateTL(&Y[(NP-1)*M*LDIM], X, &Ids[(NP-1)*M], N-(NP-1)*M);
}

