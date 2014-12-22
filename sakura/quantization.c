/*
  This file contains the files used for the computation of morton codes 

  author: Nikos Sismanis
  date: Jul 2014

*/

#include "stdio.h"
#include "stdlib.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_max.h>
#include <cilk/reducer_min.h>
#include "cilk/cilk_api.h"
#include "math.h"
#include "sys/time.h"
#include "float.h"
#include "utils.h"
#include <stdint.h>

#define DIM 3
#define LDIM 12
#define NP 128

uint32_t compute_code(float x, float low, float step){
  return floor((x - low) / step);
}

void quantizeTL(uint32_t (*restrict codes), float (*restrict X), 
		float (*restrict low), float step, int N){
  cilk_for(uint64_t i=0; i<N; i++){
    codes[i*DIM:DIM] = compute_code(X[i*LDIM:DIM], low[0:DIM], step); 
  }
}

void quantize_serial(uint32_t (*restrict codes), float (*restrict X), 
		     float (*restrict low), float step, int N){
  for(int i=0; i<N; i++){
    codes[i*DIM:DIM] = compute_code(X[i*LDIM:DIM], low[0:DIM], step); 
  }
}

void compute_quantization_codes_TL(uint32_t (*restrict codes), float (*restrict X), int N, 
				   int nbins, float (*restrict min), 
				   float (*restrict max)){
  float ranges[DIM];
  float qstep;
  ranges[:] = fabs(max[0:DIM] - min[0:DIM]);
  ranges[:] *= 1.00001;
  qstep = __sec_reduce_max(ranges[:]) / (float) nbins;
  int M = (int)ceil((float)N / (float)NP); 
  for(uint64_t i=0; i<NP; i++){
    uint64_t size = ( (i+1)*M < N ) ? M : N - i*M;
    cilk_spawn quantize_serial(&codes[i*M*DIM], &X[i*M*LDIM], min, qstep, size);
  }
}

void space_bounds(float (*restrict min), float (*restrict max), 
		  float (*restrict X), int N){
  cilk::reducer_max<float> max_val[DIM];
  cilk::reducer_min<float> min_val[DIM];
  cilk_for(int i=0; i<N; i++){
    max_val[0]->calc_max(X[i*LDIM]);
    min_val[0]->calc_min(X[i*LDIM]);
    max_val[1]->calc_max(X[i*LDIM + 1]);
    min_val[1]->calc_min(X[i*LDIM + 1]);
    max_val[2]->calc_max(X[i*LDIM + 2]);
    min_val[2]->calc_min(X[i*LDIM + 2]);
  }
  min[0] = min_val[0].get_value();
  min[1] = min_val[1].get_value();
  min[2] = min_val[2].get_value();
  max[0] = max_val[0].get_value();
  max[1] = max_val[1].get_value();
  max[2] = max_val[2].get_value();
}

