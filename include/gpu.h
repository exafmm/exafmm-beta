#ifndef postgpu_h
#define postgpu_h
#include <cutil.h>

#define CALL_GPU(EQUATION,STAGE,EVENT)\
template<>\
void Kernel<EQUATION>::STAGE() {\
  cudaThreadSynchronize();\
  startTimer("cudaMalloc   ");\
  if( keysHost.size() > keysDevcSize ) {\
    if( keysDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(keysDevc));\
    CUDA_SAFE_CALL(cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) ));\
    keysDevcSize = keysHost.size();\
  }\
  if( rangeHost.size() > rangeDevcSize ) {\
    if( rangeDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(rangeDevc));\
    CUDA_SAFE_CALL(cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) ));\
    rangeDevcSize = rangeHost.size();\
  }\
  if( sourceHost.size() > sourceDevcSize ) {\
    if( sourceDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(sourceDevc));\
    CUDA_SAFE_CALL(cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) ));\
    sourceDevcSize = sourceHost.size();\
  }\
  if( targetHost.size() > targetDevcSize ) {\
    if( targetDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(targetDevc));\
    CUDA_SAFE_CALL(cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) ));\
    targetDevcSize = targetHost.size();\
  }\
  cudaThreadSynchronize();\
  stopTimer("cudaMalloc   ");\
  startTimer("cudaMemcpy   ");\
  CUDA_SAFE_CALL(cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),      cudaMemcpyHostToDevice));\
  CUDA_SAFE_CALL(cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),     cudaMemcpyHostToDevice));\
  CUDA_SAFE_CALL(cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));\
  CUDA_SAFE_CALL(cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice));\
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(gpureal)));\
  cudaThreadSynchronize();\
  stopTimer("cudaMemcpy   ");\
  cudaThreadSynchronize();\
  startTimer(#EVENT);\
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    EQUATION##STAGE##_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  cudaThreadSynchronize();\
  stopTimer(#EVENT);\
  cudaThreadSynchronize();\
  startTimer("cudaMemcpy   ");\
  CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost));\
  cudaThreadSynchronize();\
  stopTimer("cudaMemcpy   ");\
}

#endif
