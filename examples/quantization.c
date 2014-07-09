/*
This file contains the files used for the computation of morton codes 

author: Nikos Sismanis
date: Jul 2014

*/

#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_min.h>
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"

#define DIM 3
#define LDIM 12
#define NP 128

uint compute_code(float x, float low, float step){

  return floor((x - low) / step);

}

void encode(uint *codes, float *X, float low, float step, int N){

  codes[0:N] = compute_code(X[0:N], low, step);
}

void encodeT(uint *codes, float *X, float low, float step, int N){


}

void quantize(uint *codes, float *X, float low, float step, int N){

  //int M = (int)ceil((float)N / (float)NP);
  int M = N / NP;

  for(int i=0; i<NP; i++){
      cilk_spawn encode(&codes[i*M], &X[i*M], low, step, M);
  }
  cilk_sync;

}


void quantizeT(uint *codes, float *X, float* low, float* step, int N){

  cilk_for(int i=0; i<N; i++){
    codes[i*DIM:DIM] = compute_code(X[i*DIM:DIM], low[0:DIM], step[0:DIM]);
  }

}

void quantizeTL(uint *codes, float *X, float* low, float* step, int N){

  cilk_for(int i=0; i<N; i++){
    codes[i*DIM:DIM] = compute_code(X[i*LDIM:DIM], low[0:DIM], step[0:DIM]);
  }

}

void compute_quantization_codes(uint* codes, float *X, int N, int nbins){

  float min[DIM], max[DIM], range[DIM], qstep[DIM];

  /* Compute the boundaries */ 
  cilk_for(int i=0; i<DIM; i++){ 
    max[i] = __sec_reduce_max(X[i*N:N]);
    min[i] = __sec_reduce_min(X[i*N:N]);
  }

  range[:] = fabs(max[:] - min[:]);
  range[:] += 0.01*range[:];
  qstep[:] = range[:] / nbins;
  
  for(int i=0; i<DIM; i++){
    cilk_spawn quantize(&codes[i*N], &X[i*N], min[i], qstep[i], N);
  }
  cilk_sync;
  
}

void compute_quantization_codes_T(uint* codes, float *X, int N, int nbins){

  float min[DIM], max[DIM], range[DIM], qstep[DIM];

  /* Compute the boundaries */
  cilk_for(int i=0; i<DIM; i++){ 
    max[i] = __sec_reduce_max(X[i:N:DIM]);
    min[i] = __sec_reduce_min(X[i:N:DIM]);
  }
  
  range[:] = fabs(max[:] - min[:]);
  range[:] += 0.01*range[:];
  qstep[:] = range[:] / nbins;
  
  quantizeT(codes, X, min, qstep, N);
  
}


void compute_quantization_codes_TL(uint* codes, float *X, int N, int nbins){

  float min[DIM], max[DIM], range[DIM], qstep[DIM];

  /* Compute the boundaries */
  cilk_for(int i=0; i<DIM; i++){ 
    max[i] = __sec_reduce_max(X[i:N:LDIM]);
    min[i] = __sec_reduce_min(X[i:N:LDIM]);
  }
  
  range[:] = fabs(max[:] - min[:]);
  range[:] += 0.01*range[:];
  qstep[:] = range[:] / nbins;
  
  quantizeTL(codes, X, min, qstep, N);
  
}



