#include <cilk/cilk.h>
#include "utils.h"

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

void morton_encoding_T(uint64_t (*restrict mcodes), uint32_t (*restrict codes), int N, int max_level){
  cilk_for(uint64_t i=0; i<N; i++){
    mcodes[i] = mortonEncode_magicbits(codes[i*DIM], codes[i*DIM + 1], codes[i*DIM + 2]);
  }
}


