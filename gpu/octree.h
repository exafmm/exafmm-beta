#ifndef _OCTREE_H_
#define _OCTREE_H_

#include "types.h"

typedef unsigned int uint;

class octree {
 private:
  int numBodies;
  int numLeafs;
  int numSources;
  int numTargets;
  int numLevels;
  float4 *float4buffer;
  cudaVec<uint64_t> bodyKeys;
  cudaVec<int> bodyIndex;
  cudaVec<uint2> bodyRange;
  cudaVec<uint64_t> cellKeys;
  cudaVec<int> cellLevel;
  cudaVec<uint> childRange;
  cudaVec<uint> levelOffset;
  cudaVec<uint2> levelRange;
  cudaVec<uint> validRange;
  cudaVec<uint> compactRange;

  cudaVec<uint> cellIndex;
  cudaVec<uint2> targetRange;

  cudaVec<float> openingAngle;
  cudaVec<float4> targetSizeInfo;
  cudaVec<float4> targetCenterInfo;

  cudaVec<uint>   generalBuffer1;
  float4 corner;

  cudaVec<float3> XMIN;
  cudaVec<float3> XMAX;
  cudaVec<uint>   offset;
  cudaVec<uint>   workToDo;
  
 public:
  cudaVec<float4> bodyPos;
  cudaVec<float4> bodyAcc;
  cudaVec<float4> bodyAcc2;
  cudaVec<float4> cellPos;
  cudaVec<float>  multipole;      

 private:
  bool isPowerOfTwo(const int n) {
    return ((n != 0) && !(n & (n - 1)));
  }
  void gpuCompact(cudaVec<uint> &input, cudaVec<uint> &output, int size);
  void gpuSplit(cudaVec<uint> &input, cudaVec<uint> &output, int size);
  void getBoundaries();
  void getKeys();
  void sortKeys();
  void sortBodies();
  void allocateTreePropMemory();
  void buildTree();
  void linkTree();
  void upward();
  void traverse();

 public:
  octree(const int _n) : numBodies(_n) {
    assert(isPowerOfTwo(NCRIT));
    cudaSetDevice(2);
    bodyPos.alloc(numBodies+1);
    bodyKeys.alloc(numBodies+1);
    bodyIndex.alloc(numBodies+1);
    bodyAcc.alloc(numBodies);
    bodyAcc2.alloc(numBodies);
    bodyRange.alloc(numBodies);
    cellKeys.alloc(numBodies);
    cellLevel.alloc(numBodies);
    childRange.alloc(numBodies);
    levelOffset.alloc(MAXLEVELS*2);
    levelRange.alloc(MAXLEVELS);
    validRange.alloc(2*numBodies);
    compactRange.alloc(2*numBodies);
    bodyAcc.zeros();
    CU_SAFE_CALL(cudaMalloc(&float4buffer, numBodies*sizeof(float4)));

    int treeWalkStackSize = (LMEM_STACK_SIZE * NTHREAD + 2 * NTHREAD) * NBLOCK;
    int sortBufferSize = 4 * ALIGN(numBodies,128) * 128;
    generalBuffer1.alloc(max(treeWalkStackSize,sortBufferSize));

    XMIN.alloc(64);
    XMAX.alloc(64);
    offset.alloc(NBLOCK);
    workToDo.alloc(1);

  }

  double get_time() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return double(tv.tv_sec +1.e-6*tv.tv_usec);
  }

  void iterate();
  void direct(int,int);
};

#endif
