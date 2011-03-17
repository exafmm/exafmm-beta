#ifndef gpu_h
#define gpu_h
#include "cuprintf.h"

int   *keysDevc;                                                // Keys on device
int   *rangeDevc;                                               // Ranges on device
float *sourceDevc;                                              // Sources on device
float *targetDevc;                                              // Targets on device
__device__ __constant__ float constDevc[1];                     // Constants on device

#endif
