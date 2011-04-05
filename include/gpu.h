#ifndef postgpu_h
#define postgpu_h
#include <cutil.h>

#define CALL_GPU(KERNEL,EVENT)\
void Kernel::KERNEL() {\
  cudaThreadSynchronize();\
  startTimer("cudaMalloc   ");\
  cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );\
  cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );\
  cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(float) );\
  cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(float) );\
  cudaThreadSynchronize();\
  stopTimer("cudaMalloc   ");\
  startTimer("cudaMemcpy   ");\
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),    cudaMemcpyHostToDevice);\
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),   cudaMemcpyHostToDevice);\
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(float),cudaMemcpyHostToDevice);\
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(float),cudaMemcpyHostToDevice);\
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(float));\
  cudaThreadSynchronize();\
  stopTimer("cudaMemcpy   ");\
  cudaThreadSynchronize();\
  startTimer(#EVENT);\
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    KERNEL##_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer(#EVENT);\
  cudaThreadSynchronize();\
  startTimer("cudaMemcpy   ");\
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(float),cudaMemcpyDeviceToHost);\
  cudaThreadSynchronize();\
  stopTimer("cudaMemcpy   ");\
  startTimer("cudaFree     ");\
  cudaFree(keysDevc);\
  cudaFree(rangeDevc);\
  cudaFree(targetDevc);\
  cudaFree(sourceDevc);\
  cudaThreadSynchronize();\
  stopTimer("cudaFree     ");\
}

#endif
