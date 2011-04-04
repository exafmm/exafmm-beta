#ifndef postgpu_h
#define postgpu_h
#include <cutil.h>

#define CALL_GPU(KERNEL,EVENT)\
void Kernel::KERNEL() {\
  cudaThreadSynchronize();\
  startTimer("cudaMalloc   ");\
  cudaMalloc( (void**) &keysDevc,   keysHost.size()*sizeof(int) );\
  cudaMalloc( (void**) &rangeDevc,  rangeHost.size()*sizeof(int) );\
  cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(double) );\
  cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(double) );\
  cudaThreadSynchronize();\
  stopTimer("cudaMalloc   ");\
  startTimer("cudaMemcpy   ");\
  cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int),    cudaMemcpyHostToDevice);\
  cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int),   cudaMemcpyHostToDevice);\
  cudaMemcpy(targetDevc,&targetHost[0],targetHost.size()*sizeof(double),cudaMemcpyHostToDevice);\
  cudaMemcpy(sourceDevc,&sourceHost[0],sourceHost.size()*sizeof(double),cudaMemcpyHostToDevice);\
  cudaMemcpyToSymbol(constDevc,&constHost[0],constHost.size()*sizeof(double));\
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
  cudaMemcpy(&targetHost[0],targetDevc,targetHost.size()*sizeof(double),cudaMemcpyDeviceToHost);\
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

#if Laplace
CALL_GPU(LaplaceP2M,P2M GPUkernel);
CALL_GPU(LaplaceM2M,M2M GPUkernel);
CALL_GPU(LaplaceM2L,M2L GPUkernel);
CALL_GPU(LaplaceM2P,M2P GPUkernel);
CALL_GPU(LaplaceP2P,P2P GPUkernel);
CALL_GPU(LaplaceL2L,L2L GPUkernel);
CALL_GPU(LaplaceL2P,L2P GPUkernel);
#elif BiotSavart
CALL_GPU(BiotSavartP2M,P2M GPUkernel);
CALL_GPU(BiotSavartM2M,M2M GPUkernel);
CALL_GPU(BiotSavartM2L,M2L GPUkernel);
CALL_GPU(BiotSavartM2P,M2P GPUkernel);
CALL_GPU(BiotSavartP2P,P2P GPUkernel);
CALL_GPU(BiotSavartL2L,L2L GPUkernel);
CALL_GPU(BiotSavartL2P,L2P GPUkernel);
#elif Stretching
CALL_GPU(StretchingP2M,P2M GPUkernel);
CALL_GPU(StretchingM2M,M2M GPUkernel);
CALL_GPU(StretchingM2L,M2L GPUkernel);
CALL_GPU(StretchingM2P,M2P GPUkernel);
CALL_GPU(StretchingP2P,P2P GPUkernel);
CALL_GPU(StretchingL2L,L2L GPUkernel);
CALL_GPU(StretchingL2P,L2P GPUkernel);
#elif Gaussian
CALL_GPU(GaussianP2M,P2M GPUkernel);
CALL_GPU(GaussianM2M,M2M GPUkernel);
CALL_GPU(GaussianM2L,M2L GPUkernel);
CALL_GPU(GaussianM2P,M2P GPUkernel);
CALL_GPU(GaussianP2P,P2P GPUkernel);
CALL_GPU(GaussianL2L,L2L GPUkernel);
CALL_GPU(GaussianL2P,L2P GPUkernel);
#endif

#endif
