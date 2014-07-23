/*
This file contains the files used for the computation of morton codes

author: Nikos Sismanis
date: Jul 2014

*/

#include "stdio.h"
#include "stdlib.h"
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include "math.h"

void compute_quantization_codes_T(uint* codes, float *X, int N, int nbins) {
  float Xmin[3] = {0};
  float Xmax[3] = {0};
  float X0[3];
  for (int b=0; b<N; b++) {
    for (int d=0; d<3; d++) {
      Xmin[d] = fmin(X[3*b+d],Xmin[d]);
      Xmax[d] = fmax(X[3*b+d],Xmax[d]);
    }
  }
  for (int d=0; d<3; d++) X0[d] = (Xmax[d] + Xmin[d]) / 2;
  float range = 0;
  for(int d=0; d<3; d++) {
    range = fmax(X0[d] - Xmin[d], range);
    range = fmax(Xmax[d] - X0[d], range);
  }
  range *= 1.00001;
  for(int d=0; d<3; d++) {
    Xmin[d] = X0[d] - range;
    Xmax[d] = X0[d] + range;
  }
  float d = range / nbins;
  cilk_for(int i=0; i<N; i++){
    codes[i*3:3] = floor((X[i*3:3] - Xmin[0:3]) / d);
  }
}

