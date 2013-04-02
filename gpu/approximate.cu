#include "octree.h"
#define laneIdx (threadIdx.x & (WARP_SIZE - 1))
#define warpIdx (threadIdx.x >> WARP_SIZE2)
#define IF(x) (-(int)(x))
#define ABS(x) ((int(x) < 0 ) ? -(x) : (x))

__device__ __forceinline__ int inclusiveScanInt(int* prefix, int value) 
{
  prefix[laneIdx] = value;
  for (int i = 0; i < WARP_SIZE2; i++) {
    const int offset = 1 << i;
    const int laneOffset = ABS(laneIdx-offset);
    prefix[laneIdx] += prefix[laneOffset] & IF(laneIdx >= offset);
  }
  return prefix[WARP_SIZE-1];
}

__device__ __forceinline__ int lanemask_lt()
{
  int mask;
  asm("mov.u32 %0, %lanemask_lt;" : "=r" (mask));
  return mask;
}

__device__ int exclusiveScanBit(const bool flag)
{
  const uint flags = __ballot(flag);
  return __popc(flags & lanemask_lt());
}

__device__ int reduceBit(const bool flag)
{
  const uint flags = __ballot(flag);
  return __popc(flags);
}

__device__ __forceinline__ int lanemask_le()
{
  int mask;
  asm("mov.u32 %0, %lanemask_le;" : "=r" (mask));
  return mask;
}

__device__ __forceinline__ int inclusive_segscan_warp(
    int *shmem, const int packed_value, int &dist_block, int &nseg)
{
  const int  flag = packed_value < 0;
  const int  mask = IF(flag);
  const int value = (mask & (-1-packed_value)) + (~mask & 1);
  const int flags = __ballot(flag);

  nseg += __popc(flags) ;
  dist_block = __clz(__brev(flags));

  const int distance = min(__clz(flags & lanemask_le()) + laneIdx - 31, laneIdx);
  shmem[laneIdx] = value;
  for( int i=0; i<WARP_SIZE2; i++ ) {
    const int offset = 1 << i;
    const int laneOffset = ABS(laneIdx-offset);
    shmem[laneIdx] += shmem[laneOffset] & IF(offset <= distance);
  }
  return shmem[WARP_SIZE - 1];
}

__device__ __forceinline__ int inclusive_segscan_array(int *shmem_in, const int N)
{
  int dist, nseg = 0;
  int y = inclusive_segscan_warp(shmem_in, shmem_in[laneIdx], dist, nseg);
  for( int p=WARP_SIZE; p<N; p+=WARP_SIZE ) {
    int *shmem = shmem_in + p;
    int y1 = inclusive_segscan_warp(shmem, shmem[laneIdx], dist, nseg);
    shmem[laneIdx] += y & IF(laneIdx < dist);
    y = y1;
  }
  return nseg;
}

__device__ __forceinline__ int ACCESS(const int i) {
  return (i & (LMEM_STACK_SIZE - 1)) * blockDim.x + threadIdx.x;
}

texture<uint, 1, cudaReadModeElementType> texChildRange;
texture<float, 1, cudaReadModeElementType> texOpening;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;

__device__ __forceinline__ void P2P(
    float4 &acc,  const float4 pos,
    const float4 posj) {
  const float3 dr = make_float3(posj.x - pos.x, posj.y - pos.y, posj.z - pos.z);
  const float r2     = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + EPS2;
  const float rinv   = rsqrtf(r2);
  const float rinv2  = rinv * rinv;
  const float mrinv  = posj.w * rinv;
  const float mrinv3 = mrinv * rinv2;
  acc.w += mrinv;
  acc.x += mrinv3 * dr.x;
  acc.y += mrinv3 * dr.y;
  acc.z += mrinv3 * dr.z;
}

__device__ bool applyMAC(
    const float4 sourceCenter, 
    const float4 targetCenter, 
    const float4 targetSize) {
  float3 dr = make_float3(fabsf(targetCenter.x - sourceCenter.x) - (targetSize.x),
                          fabsf(targetCenter.y - sourceCenter.y) - (targetSize.y),
                          fabsf(targetCenter.z - sourceCenter.z) - (targetSize.z));
  dr.x += fabsf(dr.x); dr.x *= 0.5f;
  dr.y += fabsf(dr.y); dr.y *= 0.5f;
  dr.z += fabsf(dr.z); dr.z *= 0.5f;
  const float ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  return ds2 <= fabsf(sourceCenter.w);
}

__device__ void traverse(
    float4 &pos_i,
    float4 &acc_i,
    float4 targetCenter,
    float4 targetSize,
    uint2 rootRange,
    int *shmem,
    int *lmem) {
  const int stackSize = LMEM_STACK_SIZE << NTHREAD2;
  int *approxSources = lmem + stackSize + 2 * WARP_SIZE * warpIdx;
  int *numDirect = shmem;
  int *stackShrd = numDirect + WARP_SIZE;
  int *directSources = stackShrd + WARP_SIZE;
  float4 *pos_j = (float4*)&directSources[3*WARP_SIZE];
  int *prefix = (int*)&pos_j[WARP_SIZE];

  // stack
  int *stackGlob = lmem;
  // begin tree-walk
  int warpOffsetApprox = 0;
  int warpOffsetDirect = 0;
  for( int root=rootRange.x; root<rootRange.y; root+=WARP_SIZE ) {
    int numSources = min(rootRange.y-root, WARP_SIZE);
    int beginStack = 0;
    int endStack = 1;
    stackGlob[threadIdx.x] = root + laneIdx;
    // walk each level
    while( numSources > 0 ) {
      int numSourcesNew = 0;
      int warpOffsetSplit = 0;
      int numStack = endStack;
      // walk a level
      for( int iStack=beginStack; iStack<endStack; iStack++ ) {
        bool valid = laneIdx < numSources;
        int node = stackGlob[ACCESS(iStack)] & IF(valid);
        numSources -= WARP_SIZE;
        float opening = tex1Dfetch(texOpening, node);
        uint childRange = tex1Dfetch(texChildRange, node);
        float4 sourceCenter = tex1Dfetch(texMultipole, node);
        sourceCenter.w = opening;
        bool split = applyMAC(sourceCenter, targetCenter, targetSize);
        bool leaf = opening <= 0;
        bool flag = split && !leaf && valid;
        int child = childRange & BODYMASK;
        int numChild = ((childRange & INVBMASK) >> LEAFBIT) & IF(flag);
        int sumChild = inclusiveScanInt(prefix, numChild);
        int laneOffset = prefix[laneIdx];
        laneOffset += warpOffsetSplit - numChild;
        for( int i=0; i<numChild; i++ )
          stackShrd[laneOffset+i] = child+i;
        warpOffsetSplit += sumChild;
        while( warpOffsetSplit >= WARP_SIZE ) {
          warpOffsetSplit -= WARP_SIZE;
          stackGlob[ACCESS(numStack)] = stackShrd[warpOffsetSplit+laneIdx];
          numStack++;
          numSourcesNew += WARP_SIZE;
          if( (numStack - iStack) > LMEM_STACK_SIZE ) return;
        }
#if 1   // APPROX
        flag = !split && valid;
        laneOffset = exclusiveScanBit(flag);
        if( flag ) approxSources[warpOffsetApprox+laneOffset] = node;
        warpOffsetApprox += reduceBit(flag);
        if( warpOffsetApprox >= WARP_SIZE ) {
          warpOffsetApprox -= WARP_SIZE;
          node = approxSources[warpOffsetApprox+laneIdx];
          pos_j[laneIdx] = tex1Dfetch(texMultipole, node);
          for( int i=0; i<WARP_SIZE; i++ )
            P2P(acc_i, pos_i, pos_j[i]);
        }
#endif
#if 1   // DIRECT
        flag = split && leaf && valid;
        const int jbody = childRange & BODYMASK;
        int numBodies = (((childRange & INVBMASK) >> LEAFBIT)+1) & IF(flag);
        directSources[laneIdx] = numDirect[laneIdx];
        int sumBodies = inclusiveScanInt(prefix, numBodies);
        laneOffset = prefix[laneIdx];
        if( flag ) prefix[exclusiveScanBit(flag)] = laneIdx;
        numDirect[laneIdx] = laneOffset;
        laneOffset -= numBodies;
        int numFinished = 0;
        while( sumBodies > 0 ) {
          numBodies = min(sumBodies, 3*WARP_SIZE-warpOffsetDirect);
          for( int i=warpOffsetDirect; i<warpOffsetDirect+numBodies; i+=WARP_SIZE )
            directSources[i+laneIdx] = 0;
          if( flag && (numDirect[laneIdx] <= numBodies) && (laneOffset >= 0) )
            directSources[warpOffsetDirect+laneOffset] = -1-jbody;
          numFinished += inclusive_segscan_array(&directSources[warpOffsetDirect], numBodies);
          numBodies = numDirect[prefix[numFinished-1]];
          sumBodies -= numBodies;
          numDirect[laneIdx] -= numBodies;
          laneOffset -= numBodies;
          warpOffsetDirect += numBodies;
          while( warpOffsetDirect >= WARP_SIZE ) {
            warpOffsetDirect -= WARP_SIZE;
            pos_j[laneIdx] = tex1Dfetch(texBody,directSources[warpOffsetDirect+laneIdx]);
            for( int i=0; i<WARP_SIZE; i++ )
              P2P(acc_i, pos_i, pos_j[i]);
          }
        }
        numDirect[laneIdx] = directSources[laneIdx];
#endif
      }

      if( warpOffsetSplit > 0 ) { 
        stackGlob[ACCESS(numStack)] = stackShrd[laneIdx];
        numStack++; 
        numSourcesNew += warpOffsetSplit;
      }
      numSources = numSourcesNew;
      beginStack = endStack;
      endStack = numStack;
    }
  }

  if( warpOffsetApprox > 0 ) {
    if( laneIdx < warpOffsetApprox )  {
      const int node = approxSources[laneIdx];
      pos_j[laneIdx] = tex1Dfetch(texMultipole, node);
    } else {
      pos_j[laneIdx] = make_float4(1.0e10f, 1.0e10f, 1.0e10f, 0.0f);
    }
    for( int i=0; i<WARP_SIZE; i++ )
      P2P(acc_i, pos_i, pos_j[i]);
  }

  if( warpOffsetDirect > 0 ) {
    if( laneIdx < warpOffsetDirect ) {
      const float4 posj = tex1Dfetch(texBody,numDirect[laneIdx]);
      pos_j[laneIdx] = posj;
    } else {
      pos_j[laneIdx] = make_float4(1.0e10f, 1.0e10f, 1.0e10f, 0.0f);
    }
    for( int i=0; i<WARP_SIZE; i++ ) 
      P2P(acc_i, pos_i, pos_j[i]);
  }
}

extern "C" __global__ void traverseKernel(
      const int numTargets,
      uint2 *levelRange,
      float4 *acc,
      float4 *targetSizeInfo,
      float4 *targetCenterInfo,
      int    *MEM_BUF,
      uint   *workToDo) {
  __shared__ int wid[NWARP];
  __shared__ int shmem_pool[10*NTHREAD];
  int *shmem = shmem_pool+10*WARP_SIZE*warpIdx;
  int *lmem = &MEM_BUF[blockIdx.x*(LMEM_STACK_SIZE*NTHREAD+2*NTHREAD)];
  while(true) {
    if( laneIdx == 0 )
      wid[warpIdx] = atomicAdd(workToDo,1);
    if( wid[warpIdx] >= numTargets ) return;
    float4 targetSize = targetSizeInfo[wid[warpIdx]];
    const int targetData = __float_as_int(targetSize.w);
    const uint begin = targetData & CRITMASK;
    const uint numTarget = ((targetData & INVCMASK) >> CRITBIT) + 1;
    float4 targetCenter = targetCenterInfo[wid[warpIdx]];
    uint body_i = begin + laneIdx % numTarget;
    float4 pos_i = tex1Dfetch(texBody,body_i);
    float4 acc_i = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
    traverse(pos_i, acc_i, targetCenter, targetSize, levelRange[2], shmem, lmem);
    if( laneIdx < numTarget )
      acc[body_i] = acc_i;
  }
}

extern "C" __global__ void directKernel(float4 *bodyPos, float4 *bodyAcc, const int N) {
  uint idx = min(blockIdx.x * blockDim.x + threadIdx.x, N-1);
  float4 pos_i = bodyPos[idx];
  float4 acc_i = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
  __shared__ float4 shmem[NTHREAD];
  float4 *pos_j = shmem + WARP_SIZE * warpIdx;
  const int numWarp = ALIGN(N, WARP_SIZE);
  for( int jwarp=0; jwarp<numWarp; jwarp++ ) {
    int jGlob = jwarp*WARP_SIZE+laneIdx;
    pos_j[laneIdx] = bodyPos[min(jGlob,N-1)];
    pos_j[laneIdx].w *= jGlob < N;
    for( int j=0; j<WARP_SIZE; j++ )
      P2P(acc_i, pos_i, pos_j[j]);
  }
  bodyAcc[idx] = acc_i;
}

void octree::traverse() {
  childRange.tex("texChildRange");
  openingAngle.tex("texOpening");
  multipole.tex("texMultipole");
  bodyPos.tex("texBody");
  workToDo.zeros();
  traverseKernel<<<NBLOCK,NTHREAD,0,execStream>>>(
    numTargets,
    levelRange.devc(),
    bodyAcc.devc(),
    targetSizeInfo.devc(),
    targetCenterInfo.devc(),
    (int*)generalBuffer1.devc(),
    workToDo.devc()
  );
}

void octree::iterate() {
  CU_SAFE_CALL(cudaStreamCreate(&execStream));
  double t1 = get_time();
  getBoundaries();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("BOUND : %lf\n",get_time() - t1);;
  t1 = get_time();
  getKeys();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("INDEX : %lf\n",get_time() - t1);;
  t1 = get_time();
  sortKeys();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("KEYS  : %lf\n",get_time() - t1);;
  t1 = get_time();
  sortBodies();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("BODIES: %lf\n",get_time() - t1);;
  t1 = get_time();
  buildTree();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("BUILD : %lf\n",get_time() - t1);;
  t1 = get_time();
  allocateTreePropMemory();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("ALLOC : %lf\n",get_time() - t1);;
  t1 = get_time();
  linkTree();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("LINK  : %lf\n",get_time() - t1);;
  t1 = get_time();
  upward();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("UPWARD: %lf\n",get_time() - t1);;
  t1 = get_time();
  traverse();
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  printf("FMM   : %lf\n",get_time() - t1);;
}

void octree::direct(int numTarget, int numBodies) {
  int blocks = ALIGN(numTarget, NTHREAD);
  directKernel<<<blocks,NTHREAD,0,execStream>>>(bodyPos.devc(),bodyAcc2.devc(),numBodies);
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  CU_SAFE_CALL(cudaStreamDestroy(execStream));
}
