/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
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
