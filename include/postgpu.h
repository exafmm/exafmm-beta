#ifndef postgpu_h
#define postgpu_h

#define CALL_GPU(KERNEL,EVENT)\
void Kernel::KERNEL() {\
  allocGPU();\
  cudaThreadSynchronize();\
  startTimer(#EVENT);\
  int numBlocks = keysHost.size();\
  KERNEL##_GPU<<< numBlocks, THREADS >>>(keysDevc,rangeDevc,targetDevc,sourceDevc);\
  cudaThreadSynchronize();\
  stopTimer(#EVENT);\
  deallocGPU();\
}

CALL_GPU(P2M,P2M GPUkernel);
CALL_GPU(M2M,M2M GPUkernel);
CALL_GPU(M2L,M2L GPUkernel);
CALL_GPU(M2P,M2P GPUkernel);
CALL_GPU(P2P,P2P GPUkernel);
CALL_GPU(L2L,L2L GPUkernel);
CALL_GPU(L2P,L2P GPUkernel);

#endif
