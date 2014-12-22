#include <cilk/cilk.h>
#include <math.h>
#include "utils.h"

#define NP 128

void quantize_serial(uint32_t (*restrict codes), float (*restrict X), 
		     float (*restrict low), float step, int N){
  for(int i=0; i<N; i++){
    codes[i*DIM:DIM] = floor((X[i*LDIM:DIM] - low[0:DIM]) / step); 
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
