#ifndef sort_h
#define sort_h
#include "types.h"

namespace exafmm {
  //! Custom radix sort for body and structures
  class Sort {
  private:
    Bodies output;                                              //!< Output buffer

  private:
    //! Radixsorts the values using the keys
    void radixsort(int * key, int * value, int size) {
      const int bitStride = 8;                                  // Number of bits in one stride
      const int stride = 1 << bitStride;                        // Size of stride in decimal
      const int mask = stride - 1;                              // Mask the bits in one stride
      int maxKey = 0;                                           // Maximum value of key
      int bucket[stride];                                       // Bucket
      int * buffer = new int [size];                            // Buffer for both key and value
      int * permutation = new int [size];                       // Permutation index
      for (int i=0; i<size; i++)                                // Loop over keys
	if (key[i] > maxKey)                                    //  If key is larger than maxKey
	  maxKey = key[i];                                      //   Update maxKey per thread
      while (maxKey > 0) {                                      // While there are bits in maxKey to process
	for (int i=0; i<stride; i++)                            //   Loop over strides
	  bucket[i] = 0;                                        //    Initialize bucket
	for (int i=0; i<size; i++)                              //   Loop over keys
	  bucket[key[i] & mask]++;                              //    Increment bucket
	for (int i=1; i<stride; i++)                            //   Loop over strides
	  bucket[i] += bucket[i-1];                             //    Scan bucket over strides
	for (int i=size-1; i>=0; i--)                           //   Loop over keys backwards
	  permutation[i] = --bucket[key[i] & mask];             //    Reverse scan bucket to get permutation
	for (int i=0; i<size; i++)                              //   Loop over values
	  buffer[permutation[i]] = value[i];                    //    Sort into buffer
	for (int i=0; i<size; i++)                              //   Loop over values
	  value[i] = buffer[i];                                 //    Copy back from buffer
	for (int i=0; i<size; i++)                              //   Loop over keys
	  buffer[permutation[i]] = key[i];                      //    Sort into buffer
	for (int i=0; i<size; i++)                              //   Loop over keys
	  key[i] = buffer[i] >> bitStride;                      //    Copy back from buffer and bit shift keys
	maxKey >>= bitStride;                                   //   Bit shift maxKey
      }                                                         //  End while for bits in maxKey
      delete[] buffer;                                          // Deallocate buffer
      delete[] permutation;                                     // Deallocate permutation index
    }

  public:
    //! Sort input accoring to ibody
    Bodies ibody(Bodies & input) {
      const int size = input.size();                            // Size of bodies vector
      int * key = new int [size];                               // Allocate key array
      int * index = new int [size];                             // Allocate index array
      for (B_iter B=input.begin(); B!=input.end(); B++) {       // Loop over input bodies
	int i = B-input.begin();                                //  Body index
	key[i] = B->IBODY;                                      //  Copy IBODY to key array
	index[i] = i;                                           //  Initialize index array
      }                                                         // End loop over input bodies
      radixsort(key,index,size);                                // Radix sort index according to key
      output.resize(size);                                      // Resize output buffer
      for (B_iter B=output.begin(); B!=output.end(); B++) {     // Loop over output boides
	int i = B-output.begin();                               //  Body index
	*B = input[index[i]];                                   //  Permute according to index
      }                                                         // End loop over output bodies
      delete[] key;                                             // Deallocate key array
      delete[] index;                                           // Deallocate index array
      return output;                                            // Return output
    }

    //! Sort input accoring to irank
    Bodies irank(Bodies & input) {
      const int size = input.size();                            // Size of bodies vector
      int * key = new int [size];                               // Allocate key array
      int * index = new int [size];                             // Allocate index array
      for (B_iter B=input.begin(); B!=input.end(); B++) {       // Loop over input bodies
	int i = B-input.begin();                                //  Body index
	key[i] = B->IRANK;                                      //  Copy IRANK to key array
	index[i] = i;                                           //  Initialize index array
      }                                                         // End loop over input bodies
      radixsort(key,index,size);                                // Radix sort index according to key
      output.resize(size);                                      // Resize output buffer
      for (B_iter B=output.begin(); B!=output.end(); B++) {     // Loop over output boides
	int i = B-output.begin();                               //  Body index
	*B = input[index[i]];                                   //  Permute according to index
      }                                                         // End loop over output bodies
      delete[] key;                                             // Deallocate key array
      delete[] index;                                           // Deallocate index array
      return output;                                            // Return output
    }

    //! Sort bodies back to original order
    Bodies unsort(Bodies & bodies) {
      bodies = ibody(bodies);                                   // Sort bodies
      return bodies;
    }
  };
}
#endif
