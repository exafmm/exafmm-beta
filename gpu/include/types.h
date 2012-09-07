#ifndef _TYPES_H_
#define _TYPES_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cuda_runtime.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <vector>
#include "cudavec.h"
#include "macros.h"
#include "vec.h"

typedef float real;
typedef vec<3,real> vec3;
typedef vec<4,real> vec4;

const int  P     = 3;
const real EPS2  = 0.0001;
const real THETA = .75;

const int MTERM = P*(P+1)*(P+2)/6;
const int LTERM = (P+1)*(P+2)*(P+3)/6;
typedef vec<MTERM,real> vecM;
typedef vec<LTERM,real> vecL;

namespace {
__host__ __device__
vec3 make_vec3(float3 input) {
  vec3 output;
  output[0] = input.x;
  output[1] = input.y;
  output[2] = input.z;
  return output;
}
__host__ __device__
vec3 make_vec3(real x, real y, real z) {
  vec3 output;
  output[0] = x;
  output[1] = y;
  output[2] = z;
  return output;
}
__host__ __device__
vec3 make_vec3(vec4 input) {
  vec3 output;
  output[0] = input[0];
  output[1] = input[1];
  output[2] = input[2];
  return output;
}
__host__ __device__
vec4 make_vec4(float4 input) {
  vec4 output;
  output[0] = input.x;
  output[1] = input.y;
  output[2] = input.z;
  output[3] = input.w;
  return output;
}
__host__ __device__
vec4 make_vec4(real x, real y, real z, real w) {
  vec4 output;
  output[0] = x;
  output[1] = y;
  output[2] = z;
  output[3] = w;
  return output;
}
}

#endif
