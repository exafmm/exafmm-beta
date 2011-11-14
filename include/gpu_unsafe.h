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

#define CALL_GPU(KERNEL,EVENT)\
void Kernel::KERNEL() {\
  if( keysHost.size() > keysDevcSize ) {\
    if( keysDevcSize != 0 ) cudaFree(keysDevc);\
    cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );\
    keysDevcSize = keysHost.size();\
  }\
  if( rangeHost.size() > rangeDevcSize ) {\
    if( rangeDevcSize != 0 ) cudaFree(rangeDevc);\
    cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );\
    rangeDevcSize = rangeHost.size();\
  }\
  if( sourceHost.size() > sourceDevcSize ) {\
    if( sourceDevcSize != 0 ) cudaFree(sourceDevc);\
    cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) );\
    sourceDevcSize = sourceHost.size();\
  }\
  if( targetHost.size() > targetDevcSize ) {\
    if( targetDevcSize != 0 ) cudaFree(targetDevc);\
    cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) );\
    targetDevcSize = targetHost.size();\
  }\
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),      cudaMemcpyHostToDevice);\
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),     cudaMemcpyHostToDevice);\
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice);\
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(gpureal),cudaMemcpyHostToDevice);\
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(gpureal));\
  startTimer(#EVENT);\
  int numBlocks = keysHost.size();\
  if( numBlocks != 0 ) {\
    KERNEL##_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  }\
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(gpureal),cudaMemcpyDeviceToHost);\
}

#endif
