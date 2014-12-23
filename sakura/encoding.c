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

inline uint64_t splitBy3(uint32_t a){
  uint64_t x = a & 0x1fffff;
  x = (x | x << 32) & 0x1f00000000ffff;
  x = (x | x << 16) & 0x1f0000ff0000ff;
  x = (x | x << 8) & 0x100f00f00f00f00f;
  x = (x | x << 4) & 0x10c30c30c30c30c3;
  x = (x | x << 2) & 0x1249249249249249;
  return x;
}

uint64_t mortonEncode_magicbits(uint32_t x, uint32_t y, uint32_t z){
  uint64_t mcodes = 0;
  mcodes |= splitBy3(x) | splitBy3(y) << 1 | splitBy3(z) << 2;
  return mcodes;
}

inline uint64_t Compact1By2(uint64_t x){
  x &= 0x1249249249249249;
  x = (x ^ (x >>  2)) & 0x10c30c30c30c30c3;
  x = (x ^ (x >>  4)) & 0x100f00f00f00f00f;
  x = (x ^ (x >>  8)) & 0x1f0000ff0000ff;
  x = (x ^ (x >> 16)) & 0x1f00000000ffff;
  x = (x ^ (x >> 32)) & 0x1fffff;
  return x;
}

void decode_morton_code(int *x, int *y, int *z, uint64_t mcode){
  x[0] = Compact1By2(mcode >> 0);
  y[0] = Compact1By2(mcode >> 1);
  z[0] = Compact1By2(mcode >> 2);
}

void morton_encoding_T(uint64_t (*restrict mcodes), uint32_t (*restrict codes), int N){
  cilk_for(uint64_t i=0; i<N; i++){
    mcodes[i] = mortonEncode_magicbits(codes[i*DIM], codes[i*DIM + 1], codes[i*DIM + 2]);
  }
}

void encodeParticles(int N, float * X, float * min, float *max, uint64_t *mcodes, int maxlev) {
  uint32_t *codes = (uint32_t *)sakura_malloc(N, DIM * sizeof(uint32_t), "Hash code array");
  int nbins = (1 << maxlev);
  start_timer();
  compute_quantization_codes_TL(codes, X, N, nbins, min, max);
  stop_timer("Quantization");
  start_timer();
  morton_encoding_T(mcodes, codes, N);
  stop_timer("Morton encoding");
  free(codes);
}
