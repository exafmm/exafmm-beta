#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"

#define LDIM 12
#define NP 128

void permuteBlock(float *Y, float *X, uint *index, int N){
  for(int i=0; i<N; i++){
    Y[i*LDIM:LDIM] = X[index[i]*LDIM:LDIM];
  }
}

void permute(float *Y, float *X, uint *index, int N){
  int M = N / NP;
  for(int i=0; i<NP-1; i++){
    cilk_spawn permuteBlock(&Y[i*M*LDIM], X, &index[i*M], M);
  }
  permuteBlock(&Y[(NP-1)*M*LDIM], X, &index[(NP-1)*M], N-(NP-1)*M);
}

