#include "octree.h"

//Warp based summation
static __device__ int inexclusive_scan_warp(volatile int *ptr,bool inclusive, const uint idx, int value) {
  const uint lane = idx & 31;
  
  if (lane >=  1) ptr[idx] = value = ptr[idx -  1]   + value;
  if (lane >=  2) ptr[idx] = value = ptr[idx -  2]   + value;
  if (lane >=  4) ptr[idx] = value = ptr[idx -  4]   + value;
  if (lane >=  8) ptr[idx] = value = ptr[idx -  8]   + value;
  if (lane >= 16) ptr[idx] = value = ptr[idx -  16]  + value;

  if(inclusive)  
    return value;    //Inclusive
  else
    return(lane > 0) ? ptr[idx-1] : 0;  //Exclusive
}

//N is number of previous blocks in the count call
extern "C" __global__ void exclusive_scan_block(uint *ptr, const int size, uint *count)  {
 if (*count == 0) return;
  const uint idx = threadIdx.x;
  const uint lane   = idx & 31;
  const uint warpid = idx >> 5;
  __shared__ int shmemESB[512];

  int value;

  //Read the data in shmem
  if(idx < size + 1)
    shmemESB[idx] = value = ptr[idx];
  else
    shmemESB[idx] = value = 0;
  
  __syncthreads();

  // step 1: Intra-warp scan in each warp
  int val = inexclusive_scan_warp(&shmemESB[0], false, idx, value);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) shmemESB[warpid] = shmemESB[idx];
  __syncthreads();
  value = shmemESB[idx];
  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inexclusive_scan_warp(&shmemESB[0], false, idx, value);
  __syncthreads();
  
  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = shmemESB[warpid - 1] + val;
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();


  // ptr[blockDim.x - 1] + lastValue; //count
  if(idx == 0)//Thread 0 saves the total count value
    *count = ptr[blockDim.x - 1];
}

static __device__  __forceinline__ int hillisSteele5(volatile uint tmp[], int &count, uint val, const int idx)
{
  tmp[idx-16] = 0;
  tmp[idx] = val;

  //Since we set half the array to 0 we don't need ifs!
  tmp[idx] = val = tmp[idx -  1]  + val;
  tmp[idx] = val = tmp[idx -  2]  + val;
  tmp[idx] = val = tmp[idx -  4]  + val;
  tmp[idx] = val = tmp[idx -  8]  + val;
  tmp[idx] = val = tmp[idx -  16] + val;

  //Inclusive sum/count
  count = tmp[blockDim.x-1];

  //Exclusive index
  return (idx > 0) ? tmp[idx-1] : 0;
}


static __device__ void reduce_block2(int tid, volatile int *shmem, int val)
{
  //Reduce the 32 block
  if(tid < 16){
    shmem[tid] = val = val + shmem[tid+16];
    shmem[tid] = val = val + shmem[tid+8];
    shmem[tid] = val = val + shmem[tid+4];
    shmem[tid] = val = val + shmem[tid+2];
    shmem[tid] = val = val + shmem[tid+1];  
  }
}


//Count the number of valid elements in this BLOCK
extern "C" __global__ void compact_count(volatile uint *values,
                                         uint *compactOffset,
                                         const int size,
                                         const uint *workToDo) {
  if (*workToDo == 0) return;
  volatile __shared__ int shmemCC[32];
  int jobs = (size / blockDim.x) / gridDim.x;
  int blocksWithExtraJobs = (size / blockDim.x) % gridDim.x;
  int extraElements = size % blockDim.x;
  int extraOffset = size - extraElements;

  //Determine the parameters and loop over the particles
  int jobSize, offset, count = 0;

  jobSize = jobs; 
  if(blockIdx.x < blocksWithExtraJobs)
    jobSize++;
  
  if(blockIdx.x <= blocksWithExtraJobs)
    offset = (jobs+1) * 32 * blockIdx.x;
  else
    offset = blocksWithExtraJobs * (jobs+1) * 32
           + (blockIdx.x - blocksWithExtraJobs) * jobs * 32;

  for(int i=0; i < jobSize; i++) {
    count  += (values[offset + threadIdx.x]   >> 31);
    offset += blockDim.x;
  }

  //Reduce to get the count of this block
  shmemCC[threadIdx.x] = count;
  __syncthreads();
  reduce_block2(threadIdx.x, &shmemCC[0], count);
  
  //Save the values / count of the current block
  if(threadIdx.x == 0)  
    compactOffset[blockIdx.x] = shmemCC[0];

  if(blockIdx.x == gridDim.x-1) {
    //Here i use single element reads for ease of boundary conditions and steps
    count = 0;
    int offset = extraOffset + threadIdx.x;   

    for( int i=offset; i<size; i+=blockDim.x )
      count += values[i] >> 31;

    //Reduce
    shmemCC[threadIdx.x] = count;
    __syncthreads();
    reduce_block2(threadIdx.x, &shmemCC[0], count);
  
    //Save the count
    if(threadIdx.x == 0)  
      compactOffset[gridDim.x] = shmemCC[0];

  }
}

//The kernel that actually moves the data
extern "C" __global__ void compact_move(uint *values,
                                        uint *output,
                                        uint *compactOffset,
                                        const int size,
                                        const uint *workToDo) {
  if(*workToDo == 0) return;
  //Determine the parameters and loop over the particles
  int jobSize, offset;
  volatile __shared__ uint shmemCM[48];
  int jobs = (size / blockDim.x) / gridDim.x;
  int blocksWithExtraJobs = (size / blockDim.x) % gridDim.x;
  int extraElements = size % blockDim.x;
  int extraOffset = size - extraElements;

  jobSize = jobs; 
  if(blockIdx.x < blocksWithExtraJobs)
    jobSize++;
  
  if(blockIdx.x <= blocksWithExtraJobs)
    offset = (jobs+1) * 32 * blockIdx.x;
  else
  {
    offset = blocksWithExtraJobs*(jobs+1)*32;
    offset += (blockIdx.x - blocksWithExtraJobs) * jobs * 32;
  }

  int outputOffset = compactOffset[blockIdx.x];
  int curCount;

  //Do per step the prefix scan to determine the output locations
  for(int i=0; i < jobSize; i++)
  {     
    uint validBase = values[offset + threadIdx.x];
    int value = (validBase >> 31);  

    int idx = hillisSteele5(&shmemCM[16], curCount, value, threadIdx.x);

    if(value)
      output[idx + outputOffset] = validBase & 0x7FFFFFFF;

    outputOffset += curCount;
    offset       += blockDim.x;
  }
  
  //Block 0 handles any extra elements that couldn't be divided equally
  if(blockIdx.x == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    offset       = extraOffset;   
    outputOffset = compactOffset[gridDim.x]; 

    for(int i=0; i < extraElements;  i += blockDim.x)
    {
      if( (offset + i +  threadIdx.x) < size ) {
        uint validBase = values[offset + i +  threadIdx.x];     
        int value = validBase >> 31;
        int idx  = hillisSteele5(&shmemCM[16], curCount, value, threadIdx.x);
        if(value)
          output[idx + outputOffset] = validBase & 0x7FFFFFFF; 
      }
      outputOffset += curCount;
    }
  }//end if bid==0 
}//end compact_move

//The kernel that actually moves/splits the data
extern "C" __global__ void split_move(uint *valid,
                                      uint *output,
                                      uint *compactOffset,
                                      const int size) {
  //Determine the parameters and loop over the particles
  int jobSize, offset;
  volatile __shared__ uint shmemSM[48];
  int jobs = (size / blockDim.x) / gridDim.x;
  int blocksWithExtraJobs = (size / blockDim.x) % gridDim.x;
  int extraElements = size % blockDim.x;
  int extraOffset = size - extraElements;

  jobSize = jobs; 
  if(blockIdx.x < blocksWithExtraJobs)
    jobSize++;
  
  if(blockIdx.x <= blocksWithExtraJobs)
    offset = (jobs+1) * 32 * blockIdx.x;
  else
  {
    offset = blocksWithExtraJobs * (jobs+1) * 32;
    offset += (blockIdx.x - blocksWithExtraJobs) * jobs * 32;
  }

  int outputOffset = compactOffset[blockIdx.x];

  int rightOutputOffset = compactOffset[gridDim.x+1];
  rightOutputOffset     = rightOutputOffset + offset - outputOffset;

  int curCount;
  int idx, ridx;

  //Do per step the prefix scan to determine the output locations
  for(int i=0; i < jobSize; i++)
  {     
    uint validBase = valid[offset + threadIdx.x];
    int value = (validBase >> 31);  
    idx  = hillisSteele5(&shmemSM[16], curCount, value, threadIdx.x);
    ridx = threadIdx.x - idx;
    
    if((validBase >> 31))
    {
      output[idx + outputOffset] = validBase & 0x7FFFFFFF;
      idx++;
    }
    else
    {
      output[ridx + rightOutputOffset] = validBase & 0x7FFFFFFF;      
      ridx++;
    }

    outputOffset      += curCount;
    rightOutputOffset += 32 - curCount; //64 (32*2) since we do 2 items a time
    offset            += blockDim.x;    //Step to the next N threads
  }
  
  //Block 0 handles any extra elements that couldn't be divided equally
  if(blockIdx.x == 0)
  {
    //Here i use single element reads for ease of boundary conditions and steps
    offset              = extraOffset;   
    outputOffset        = compactOffset[gridDim.x];    
    rightOutputOffset   = compactOffset[gridDim.x+1];
    rightOutputOffset   = rightOutputOffset + offset - outputOffset;

    uint* valid2 = (uint*) valid;
    
    for(int i=0; i < extraElements;  i += blockDim.x)
    {
      uint value = 0;    
      if((offset + i +  threadIdx.x) < size)  //Make sure we dont read more than there are items      
        value = valid2[offset + i +  threadIdx.x];     

      idx  = hillisSteele5(&shmemSM[16], curCount, value >> 31, threadIdx.x);
      ridx = threadIdx.x - idx;

      if((offset + i +  threadIdx.x) < size)  
        if(value >> 31)
          output[idx + outputOffset]       = value & 0x7FFFFFFF; 
        else
          output[ridx + rightOutputOffset] = value & 0x7FFFFFFF;

      outputOffset += curCount;
      rightOutputOffset += 32 - curCount;
    }
  }//end if bid==0 
}//end split_move

void octree::gpuCompact(cudaVec<uint> &input, cudaVec<uint> &output, int size) {
  int blocks = NBLOCK-2;
  compact_count<<<blocks,32,0,execStream>>>(input.devc(),offset.devc(),size,workToDo.devc());
  exclusive_scan_block<<<1,NBLOCK,0,execStream>>>(offset.devc(),blocks,workToDo.devc());
  compact_move<<<blocks,32,0,execStream>>>(input.devc(),output.devc(),offset.devc(),size,workToDo.devc());
}

void octree::gpuSplit(cudaVec<uint> &input, cudaVec<uint> &output, int size) {
  int blocks = NBLOCK-2;
  compact_count<<<blocks,32,0,execStream>>>(input.devc(),offset.devc(),size,workToDo.devc());
  exclusive_scan_block<<<1,NBLOCK,0,execStream>>>(offset.devc(),blocks,workToDo.devc());
  split_move<<<blocks,32,0,execStream>>>(input.devc(),output.devc(),offset.devc(),size);
}
