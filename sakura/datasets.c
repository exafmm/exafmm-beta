/*

  This file contains the different distributions used for testing the octree construction

  author: Nikos Sismanis
  date: Jul 2014

*/

#include "stdio.h" 
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include <stdint.h>

#define DIM 3 // Dimension of the space 
#define PI 3.1415927 

#ifndef LDIM
#define LDIM 12 // Leading dimension in memory of the FMM data structure
#endif


void cubeTL(float *X, int N){

  srand(time(NULL));
  for(uint64_t i=0; i<N; i++){
    for(uint64_t j=0; j<DIM; j++){
      X[i*LDIM + j] = (float)rand() / (float) RAND_MAX;
      //X[i*LDIM + j] = ((float)rand() / (float) RAND_MAX) * PI - PI;
    }
    X[i*LDIM + 3] = 1;
    X[i*LDIM + LDIM-1] = i;
  }
  
}

void octant_uniform(float *X, int N){

  srand48(NULL);
  for(uint64_t i=0; i<N; i++){
    float th = (PI/2) * drand48();
    float ph = (PI/2) * drand48();
    X[i*LDIM] = sin(ph);
    X[i*LDIM + 1] = cos(th) * cos(ph);
    X[i*LDIM + 2] = sin(th) * cos(ph);
    X[i*LDIM + 3] = 1;
    X[i*LDIM + LDIM-1] = i;
  }


}


void plummer(float *X, int N){

  srand48(NULL);

  for(uint64_t i=0; i<N; i++){
    float X1 = drand48();
    float X2 = drand48();
    float X3 = drand48();
    float R =  1.0 / sqrt( (pow(X1, -2.0 / 3.0) - 1.0) );
    if(R<100.0){
      float x1 = (1.0 - 2.0 * X2) * R;
      float x2 = sqrt(R * R - X1 * X1) * cos(2.0 * M_PI * X3);
      float x3 = sqrt(R * R - X1 * X1) * sin(2.0 * M_PI * X3);
      float scale = 3.0 * M_PI / 16.0;
      x1 *= scale; x2 *= scale; x3 *= scale;
      X[i*LDIM] = x1;
      X[i*LDIM + 1] = x2;
      X[i*LDIM + 2] = x3;
    }
    X[i*LDIM + 3] = 1;
    X[i*LDIM + LDIM-1] = i;
  }

}

/* Function that creates a dataset for testing 0 (spherical octant) 
   1 (uniform cude) */
void create_dataset_TL(float *X, int N, int dist){

  switch(dist){
  case 0:
    octant_uniform(X, N);
    break;
  case 1:
    cubeTL(X, N);
    break;
  case 2:
    octant_uniform(X, N);
    break;
  case 3: 
    plummer(X, N);
    break;
  default:
    octant_uniform(X, N);
    break;
  }

}
