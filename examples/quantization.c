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
#include <iostream>
#include "math.h"
#include "sys/time.h"

#define DIM 3

#ifndef LDIM
#define LDIM 12 // Apparent dimension
#endif

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

void compute_quantization_codes_T(uint* codes, float *X, int N, int nbins){

  float Xmin[DIM] = {0};
  float Xmax[DIM] = {0};
  float X0[DIM];
  for (int b=0; b<N; b++){
    for (int d=0; d<3; d++) {
      Xmin[d] = fmin(X[3*b+d],Xmin[d]);
      Xmax[d] = fmax(X[3*b+d],Xmax[d]);
    }
  }
  for (int d=0; d<3; d++) X0[d] = (Xmax[d] + Xmin[d]) / 2;
  float range = 0;
  for(int d=0; d<DIM; d++) {
    range = fmax(X0[d] - Xmin[d], range);
    range = fmax(Xmax[d] - X0[d], range);
  }
  range *= 1.00001;
  for(int d=0; d<DIM; d++) {
    Xmin[d] = X0[d] - range;
    Xmax[d] = X0[d] + range;
  }
  for(int d=0; d<DIM; d++){
    float qstep = range / nbins;
    cilk_spawn quantize(&codes[d*N], &X[d*N], Xmin[d], qstep, N);
  }
  cilk_sync;

}

void compute_quantization_codes(uint* codes, float *X, int N, int nbins){

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



