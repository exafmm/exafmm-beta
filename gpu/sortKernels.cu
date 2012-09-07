#include <b40c/radix_sort/enactor.cuh>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/device_vector.h> 
#include <thrust/iterator/transform_iterator.h>
#include "octree.h"

Sort90::Sort90(uint size, uint* generalBuffer) 
{
  double_buffer = new b40c::util::DoubleBuffer<uint, uint>;
  sort_enactor = new b40c::radix_sort::Enactor;
  int offset = 0;
  int stride = ALIGN(size,128) * 128;
  double_buffer->d_keys[0] = generalBuffer+offset;
  offset += stride;
  double_buffer->d_keys[1] = generalBuffer+offset;
  offset += stride;
  double_buffer->d_values[0] = generalBuffer+offset;
  offset += stride;
  double_buffer->d_values[1] = generalBuffer+offset;
}

Sort90::~Sort90() 
{
  delete double_buffer;
  delete sort_enactor;
}

template<typename KeyPtr, typename PermutationPtr, typename OutputPtr>
void apply_permutation(KeyPtr& thrust_in,
                       PermutationPtr& thrustPermutation,
                       OutputPtr& thrust_out,
                       int size)
{
  // permute the keys into out vector
  thrust::gather(thrustPermutation, thrustPermutation + size, thrust_in, thrust_out);
}


// Extract 32-bit word from uint4
template<int keyIdx>
struct ExtractBits: public thrust::unary_function<uint4, uint>
{
  __host__ __device__ __forceinline__ uint operator()(uint4 key) const
  {
    if      (keyIdx == 0) return key.x;
    else if (keyIdx == 1) return key.y;
    else                  return key.z;
  }
};


template<int keyIdx, typename KeyPtr>
void update_permutation(KeyPtr& thrustInput, 
                        int size,
                        b40c::util::DoubleBuffer<uint, uint> &double_buffer,
		                      b40c::radix_sort::Enactor &sort_enactor)
{
  thrust::device_ptr<uint> thrustPermutation = 
    thrust::device_pointer_cast(double_buffer.d_values[double_buffer.selector]);

  // thrust ptr to temporary 32-bit keys
  thrust::device_ptr<uint> thust_32bit_temp = 
    thrust::device_pointer_cast(double_buffer.d_keys[double_buffer.selector]);

  // gather into temporary keys with the current reordering
  thrust::gather(thrustPermutation,
                 thrustPermutation + size,
                 thrust::make_transform_iterator(thrustInput, ExtractBits<keyIdx>()),
                 thust_32bit_temp);

  // Stable-sort the top 30 bits of the temp keys (and
  // associated thrustPermutation values)
  sort_enactor.Sort<30, 0>(double_buffer, size);
}


  // Back40 90-bit sorting: sorts the lower 30 bits in uint4's key
void Sort90::sort(uint4 *input,
                  cudaVec<uint4> &output,
                  int size)
{
  // thrust ptr to input
  thrust::device_ptr<uint4> thrustInput = 
    thrust::device_pointer_cast(input);

  // thrust ptr to output
  thrust::device_ptr<uint4> thrustOutput = 
    thrust::device_pointer_cast(output.devc());

  // thrust ptr to permutation buffer
  thrust::device_ptr<uint> thrustPermutation = 
    thrust::device_pointer_cast(double_buffer->d_values[double_buffer->selector]);

  // initialize values (thrustPermutation) to [0, 1, 2, ... ,size-1]
  thrust::sequence(thrustPermutation, thrustPermutation + size);

  // sort z, y, x
  // careful: note 2, 1, 0 key word order, NOT 0, 1, 2.
  update_permutation<2>(thrustInput, size, *double_buffer, *sort_enactor);
  update_permutation<1>(thrustInput, size, *double_buffer, *sort_enactor);
  update_permutation<0>(thrustInput, size, *double_buffer, *sort_enactor);

  // refresh thrust ptr to permutation buffer (may have changed inside ping-pong)
  thrustPermutation = 
    thrust::device_pointer_cast(double_buffer->d_values[double_buffer->selector]);

  // Note: thrustPermutation now maps unsorted keys to sorted order
  apply_permutation(thrustInput, thrustPermutation, thrustOutput, size);
}
