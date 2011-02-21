#ifndef gpu_h
#define gpu_h
#include <sys/time.h>

int const GPUS(4);                                              // Number of GPUs per node
int const THREADS(256);                                         // Number of threads per thread-block

int hostConstant[4];                                            // Constants on host
__device__ __constant__ int deviceConstant[4];                  // Constants on device

double get_gpu_time(void) {
  cudaThreadSynchronize();
  struct timeval tv;                                            // Time value
  gettimeofday(&tv, NULL);                                      // Get time of day in seconds and microseconds
  return double(tv.tv_sec+tv.tv_usec*1e-6);                     // Combine seconds and microseconds and return
}

#endif
