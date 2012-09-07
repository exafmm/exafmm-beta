#ifndef _CUDAVEC_H_
#define _CUDAVEC_H_
#include <cstdio>

#define CU_SAFE_CALL(err)  __checkCudaErrors (err, __FILE__, __LINE__)
inline void __checkCudaErrors(cudaError err, const char *file, const int line ) {
  if(cudaSuccess != err) {
    fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
    exit(-1);
  }
}

template<class T>
class cudaVec {
private:
  int size;
  T *devcPtr;
  T *hostPtr;

public:
  cudaVec() : size(0), devcPtr(NULL), hostPtr(NULL) {}

  ~cudaVec() {
    cudaFree(devcPtr);
#if PINNED
    cudaFreeHost((void*)hostPtr);
#else
    free(hostPtr);
#endif
  }

  void alloc(int n) {
    assert(size == 0);
    size = n;
#if PINNED
    CU_SAFE_CALL(cudaMallocHost((T**)&hostPtr, size*sizeof(T)));
#else
    hostPtr = (T*)malloc(size*sizeof(T));
#endif
    CU_SAFE_CALL(cudaMalloc((T**)&devcPtr, size*sizeof(T)));
  }

  void zeros() {
    CU_SAFE_CALL(cudaMemset((void*)devcPtr, 0, size*sizeof(T)));
  }

  void ones() {
    CU_SAFE_CALL(cudaMemset((void*)devcPtr, 1, size*sizeof(T)));
  }

  void d2h() {
    CU_SAFE_CALL(cudaMemcpy(hostPtr, devcPtr, size*sizeof(T), cudaMemcpyDeviceToHost));
  }

  void d2h(int n) {
    CU_SAFE_CALL(cudaMemcpy(hostPtr, devcPtr, n*sizeof(T),cudaMemcpyDeviceToHost));
  }

  void h2d() {
    CU_SAFE_CALL(cudaMemcpy(devcPtr, hostPtr, size*sizeof(T),cudaMemcpyHostToDevice ));
  }

  void h2d(int n) {
    CU_SAFE_CALL(cudaMemcpy(devcPtr, hostPtr, n*sizeof(T),cudaMemcpyHostToDevice));
  }

  void tex(const char *symbol) {
    const textureReference *texref;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
    CU_SAFE_CALL(cudaGetTextureReference(&texref,symbol));
    CU_SAFE_CALL(cudaBindTexture(0,texref,(void*)devcPtr,&channelDesc,sizeof(T)*size));
  }

  T& operator[] (int i){ return hostPtr[i]; }
  T* devc() {return devcPtr;}
};

#endif
