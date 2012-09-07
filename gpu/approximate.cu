#include "octree.h"
#define laneId (threadIdx.x & (WARP_SIZE - 1))
#define warpId (threadIdx.x >> WARP_SIZE2)
#define IF(x) (-(int)(x))
#define ABS(x) ((int(x) < 0 ) ? -(x) : (x))

__device__ __forceinline__ int inclusiveScanInt(int* prefix, int value) 
{
  prefix[laneId] = value;
  for (int i = 0; i < WARP_SIZE2; i++) {
    const int offset = 1 << i;
    const int laneOffset = ABS(laneId-offset);
    prefix[laneId] += prefix[laneOffset] & IF(laneId >= offset);
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

  const int distance = min(__clz(flags & lanemask_le()) + laneId - 31, laneId);
  shmem[laneId] = value;
  for( int i=0; i<WARP_SIZE2; i++ ) {
    const int offset = 1 << i;
    const int laneOffset = ABS(laneId-offset);
    shmem[laneId] += shmem[laneOffset] & IF(offset <= distance);
  }
  return shmem[WARP_SIZE - 1];
}

__device__ __forceinline__ int inclusive_segscan_array(int *shmem_in, const int N)
{
  int dist, nseg = 0;
  int y = inclusive_segscan_warp(shmem_in, shmem_in[laneId], dist, nseg);
  for( int p=WARP_SIZE; p<N; p+=WARP_SIZE ) {
    int *shmem = shmem_in + p;
    int y1 = inclusive_segscan_warp(shmem, shmem[laneId], dist, nseg);
    shmem[laneId] += y & IF(laneId < dist);
    y = y1;
  }
  return nseg;
}

__device__ __forceinline__ int ACCESS(const int i) {
  return (i & (LMEM_STACK_SIZE - 1)) * blockDim.x + threadIdx.x;
}

texture<uint, 1, cudaReadModeElementType> texNodeChild;
texture<float, 1, cudaReadModeElementType> texOpening;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;

__device__ __forceinline__ void P2P(
    float4 &acc,  const float4 pos,
    const float4 posj) {
  const float3 dr = make_float3(posj.x - pos.x, posj.y - pos.y, posj.z - pos.z);
  const float r2     = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + EPS2;
  const float rinv   = rsqrtf(r2);
  const float rinv2  = rinv*rinv;
  const float mrinv  = posj.w * rinv;
  const float mrinv3 = mrinv * rinv2;
  acc.w -= mrinv;
  acc.x += mrinv3 * dr.x;
  acc.y += mrinv3 * dr.y;
  acc.z += mrinv3 * dr.z;
}

__device__ bool applyMAC(
    const float4 sourceCenter, 
    const float4 groupCenter, 
    const float4 groupSize) {
  float3 dr = make_float3(fabsf(groupCenter.x - sourceCenter.x) - (groupSize.x),
                          fabsf(groupCenter.y - sourceCenter.y) - (groupSize.y),
                          fabsf(groupCenter.z - sourceCenter.z) - (groupSize.z));
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
  int *approxNodes = lmem + stackSize + 2 * WARP_SIZE * warpId;
  int *numDirect = shmem;
  int *stackShrd = numDirect + WARP_SIZE;
  int *directNodes = stackShrd + WARP_SIZE;
  float4 *pos_j = (float4*)&directNodes[3*WARP_SIZE];
  int *prefix = (int*)&pos_j[WARP_SIZE];

  // stack
  int *stackGlob = lmem;
  // begin tree-walk
  int warpOffsetApprox = 0;
  int warpOffsetDirect = 0;
  for( int root=rootRange.x; root<rootRange.y; root+=WARP_SIZE ) {
    int numNodes = min(rootRange.y-root, WARP_SIZE);
    int beginStack = 0;
    int endStack = 1;
    stackGlob[threadIdx.x] = root + laneId;
    // walk each level
    while( numNodes > 0 ) {
      int numNodesNew = 0;
      int warpOffsetSplit = 0;
      int numStack = endStack;
      // walk a level
      for( int iStack=beginStack; iStack<endStack; iStack++ ) {
        bool valid = laneId < numNodes;
        int node = stackGlob[ACCESS(iStack)] & IF(valid);
        numNodes -= WARP_SIZE;
        float opening = tex1Dfetch(texOpening, node);
        uint sourceData = tex1Dfetch(texNodeChild, node);
        float4 sourceCenter = tex1Dfetch(texMultipole, node);
        sourceCenter.w = opening;
        bool split = applyMAC(sourceCenter, targetCenter, targetSize);
        bool leaf = opening <= 0;
        bool flag = split && !leaf && valid;
        int child = sourceData & 0x0FFFFFFF;
        int numChild = ((sourceData & 0xF0000000) >> 28) & IF(flag);
        int sumChild = inclusiveScanInt(prefix, numChild);
        int laneOffset = prefix[laneId];
        laneOffset += warpOffsetSplit - numChild;
        for( int i=0; i<numChild; i++ )
          stackShrd[laneOffset+i] = child+i;
        warpOffsetSplit += sumChild;
        while( warpOffsetSplit >= WARP_SIZE ) {
          warpOffsetSplit -= WARP_SIZE;
          stackGlob[ACCESS(numStack)] = stackShrd[warpOffsetSplit+laneId];
          numStack++;
          numNodesNew += WARP_SIZE;
          if( (numStack - iStack) > LMEM_STACK_SIZE ) return;
        }
#if 1   // APPROX
        flag = !split && valid;
        laneOffset = exclusiveScanBit(flag);
        if( flag ) approxNodes[warpOffsetApprox+laneOffset] = node;
        warpOffsetApprox += reduceBit(flag);
        if( warpOffsetApprox >= WARP_SIZE ) {
          warpOffsetApprox -= WARP_SIZE;
          node = approxNodes[warpOffsetApprox+laneId];
          pos_j[laneId] = tex1Dfetch(texMultipole, node);
          for( int i=0; i<WARP_SIZE; i++ )
            P2P(acc_i, pos_i, pos_j[i]);
        }
#endif
#if 1   // DIRECT
        flag = split && leaf && valid;
        const int jbody = sourceData & BODYMASK;
        int numBodies = (((sourceData & INVBMASK) >> LEAFBIT)+1) & IF(flag);
        directNodes[laneId] = numDirect[laneId];

        int sumBodies = inclusiveScanInt(prefix, numBodies);
        laneOffset = prefix[laneId];
        if( flag ) prefix[exclusiveScanBit(flag)] = laneId;
        numDirect[laneId] = laneOffset;
        laneOffset -= numBodies;
        int numFinished = 0;
        while( sumBodies > 0 ) {
          numBodies = min(sumBodies, 3*WARP_SIZE-warpOffsetDirect);
          for( int i=warpOffsetDirect; i<warpOffsetDirect+numBodies; i+=WARP_SIZE )
            directNodes[i+laneId] = 0;
          if( flag && (numDirect[laneId] <= numBodies) && (laneOffset >= 0) )
            directNodes[warpOffsetDirect+laneOffset] = -1-jbody;
          numFinished += inclusive_segscan_array(&directNodes[warpOffsetDirect], numBodies);
          numBodies = numDirect[prefix[numFinished-1]];
          sumBodies -= numBodies;
          numDirect[laneId] -= numBodies;
          laneOffset -= numBodies;
          warpOffsetDirect += numBodies;
          while( warpOffsetDirect >= WARP_SIZE ) {
            warpOffsetDirect -= WARP_SIZE;
            pos_j[laneId] = tex1Dfetch(texBody,directNodes[warpOffsetDirect+laneId]);
            for( int i=0; i<WARP_SIZE; i++ )
              P2P(acc_i, pos_i, pos_j[i]);
          }
        }
        numDirect[laneId] = directNodes[laneId];
#endif
      }

      if( warpOffsetSplit > 0 ) { 
        stackGlob[ACCESS(numStack)] = stackShrd[laneId];
        numStack++; 
        numNodesNew += warpOffsetSplit;
      }
      numNodes = numNodesNew;
      beginStack = endStack;
      endStack = numStack;
    }
  }

  if( warpOffsetApprox > 0 ) {
    if( laneId < warpOffsetApprox )  {
      const int node = approxNodes[laneId];
      pos_j[laneId] = tex1Dfetch(texMultipole, node);
    } else {
      pos_j[laneId] = make_float4(1.0e10f, 1.0e10f, 1.0e10f, 0.0f);
    }
    for( int i=0; i<WARP_SIZE; i++ )
      P2P(acc_i, pos_i, pos_j[i]);
  }

  if( warpOffsetDirect > 0 ) {
    if( laneId < warpOffsetDirect ) {
      const float4 posj = tex1Dfetch(texBody,numDirect[laneId]);
      pos_j[laneId] = posj;
    } else {
      pos_j[laneId] = make_float4(1.0e10f, 1.0e10f, 1.0e10f, 0.0f);
    }
    for( int i=0; i<WARP_SIZE; i++ ) 
      P2P(acc_i, pos_i, pos_j[i]);
  }
}

extern "C" __global__ void
  traverseKernel(
      const int numGroups,
      uint2 *levelRange,
      float4 *acc,
      float4 *groupSizeInfo,
      float4 *groupCenterInfo,
      int    *MEM_BUF,
      uint   *workToDo) {
  __shared__ int wid[4];
  __shared__ int shmem_pool[10*NTHREAD];
  int *shmem = shmem_pool+10*WARP_SIZE*warpId;
  int *lmem = &MEM_BUF[blockIdx.x*(LMEM_STACK_SIZE*NTHREAD+2*NTHREAD)];
  while(true) {
    if( laneId == 0 )
      wid[warpId] = atomicAdd(workToDo,1);
    if( wid[warpId] >= numGroups ) return;
    float4 groupSize = groupSizeInfo[wid[warpId]];
    const int groupData = __float_as_int(groupSize.w);
    const uint begin = groupData & CRITMASK;
    const uint numGroup = ((groupData & INVCMASK) >> CRITBIT) + 1;
    float4 groupCenter = groupCenterInfo[wid[warpId]];
    uint body_i = begin + laneId % numGroup;
    float4 pos_i = tex1Dfetch(texBody,body_i);
    float4 acc_i = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

    traverse(pos_i, acc_i, groupCenter, groupSize, levelRange[2], shmem, lmem);
    if( laneId < numGroup )
      acc[body_i] = acc_i;
  }
}

extern "C" __global__ void directKernel(float4 *bodyPos, float4 *bodyAcc, const int N) {
  uint idx = min(blockIdx.x * blockDim.x + threadIdx.x, N-1);
  float4 pos_i = bodyPos[idx];
  float4 acc_i = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
  __shared__ float4 shmem[NTHREAD];
  float4 *pos_j = shmem + WARP_SIZE * warpId;
  const int numWarp = ALIGN(N, WARP_SIZE);
  for( int jwarp=0; jwarp<numWarp; jwarp++ ) {
    int jGlob = jwarp*WARP_SIZE+laneId;
    pos_j[laneId] = bodyPos[min(jGlob,N-1)];
    pos_j[laneId].w *= jGlob < N;
    for( int i=0; i<WARP_SIZE; i++ )
      P2P(acc_i, pos_i, pos_j[i]);
  }
  bodyAcc[idx] = acc_i;
}

void octree::traverse() {
  nodeChild.tex("texNodeChild");
  openingAngle.tex("texOpening");
  multipole.tex("texMultipole");
  bodyPos.tex("texBody");
  workToDo.zeros();
  traverseKernel<<<NBLOCK,NTHREAD,0,execStream>>>(
    numGroups,
    levelRange.devc(),
    bodyAcc.devc(),
    groupSizeInfo.devc(),
    groupCenterInfo.devc(),
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

void octree::direct() {
  int blocks = ALIGN(numBodies/100, NTHREAD);
  directKernel<<<blocks,NTHREAD,0,execStream>>>(bodyPos.devc(),bodyAcc2.devc(),numBodies);
  CU_SAFE_CALL(cudaStreamSynchronize(execStream));
  CU_SAFE_CALL(cudaStreamDestroy(execStream));
}
