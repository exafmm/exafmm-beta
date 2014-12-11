#ifndef sort_h
#define sort_h
#include "types.h"
#ifndef _OPENMP
int omp_get_num_threads() {return 1;}
int omp_get_thread_num() {return 0;}
#else
#include <omp.h>
#endif

//! Custom radix sort for body and structures
class Sort {
private:
  Bodies output;                                                //!< Output buffer

private:
  //! Radixsorts the values using the keys
  void radixsort(int * key, int * value, int size) {
    const int bitStride = 8;                                    // Number of bits in one stride
    const int stride = 1 << bitStride;                          // Size of stride in decimal
    const int mask = stride - 1;                                // Mask the bits in one stride
    int numThreads;                                             // Number of OpenMP threads
    int maxKey = 0;                                             // Maximum value of key
    int (*bucketPerThread)[stride];                             // Bucket for each thread
    int * maxKeyPerThread;                                      // MaxKey for each thread
    int * buffer = new int [size];                              // Buffer for both key and value
    int * permutation = new int [size];                         // Permutation index
#pragma omp parallel
    {                                                           // Start OpenMP parallel clause
      numThreads = omp_get_num_threads();                       //  Get number of OpenMP threads
#pragma omp single
      {                                                         //  Start serial clause
	bucketPerThread = new int [numThreads][stride]();       //  Allocate bucket per thread
	maxKeyPerThread = new int [numThreads];                 //  Allocate maxKey per thread
	for (int i=0; i<numThreads; i++)                        //  Loop over threads
	  maxKeyPerThread[i] = 0;                               //   Initialize maxKey per thread
      }                                                         //  End serial clause
#pragma omp for
      for( int i=0; i<size; i++ )                               //  Loop over keys
	if( key[i] > maxKeyPerThread[omp_get_thread_num()] )    //   If key is larger than maxKey
	  maxKeyPerThread[omp_get_thread_num()] = key[i];       //    Update maxKey per thread
#pragma omp single
      for( int i=0; i<numThreads; i++ )                         //  Loop over threads
	if( maxKeyPerThread[i] > maxKey ) maxKey = maxKeyPerThread[i];// Update maxKey
      while( maxKey > 0 ) {                                     //  While there are bits in maxKey to process
	int bucket[stride] = {0};                               //   Initialize bucket
#pragma omp single
	for( int t=0; t<numThreads; t++ )                       //   Loop over threads
	  for( int i=0; i<stride; i++ )                         //    Loop over strides
	    bucketPerThread[t][i] = 0;                          //     Initialize bucket per thread
#pragma omp for
	for( int i=0; i<size; i++ )                             //   Loop over keys
	  bucketPerThread[omp_get_thread_num()][key[i] & mask]++;//   Increment bucket
#pragma omp single
	{                                                       //   Start serial clause
	  for( int t=0; t<numThreads; t++ )                     //    Loop over threads
	    for( int i=0; i<stride; i++ )                       //     Loop over strides
	      bucket[i] += bucketPerThread[t][i];               //      Update bucket from all threads
	  for( int i=1; i<stride; i++ )                         //    Loop over strides
	    bucket[i] += bucket[i-1];                           //     Scan bucket over strides
	  for( int i=size-1; i>=0; i-- )                        //    Loop over keys backwards
	    permutation[i] = --bucket[key[i] & mask];           //     Reverse scan bucket to get permutation
	}                                                       //   End serial clause
#pragma omp for
	for( int i=0; i<size; i++ )                             //   Loop over values
	  buffer[permutation[i]] = value[i];                    //    Sort into buffer
#pragma omp for
	for( int i=0; i<size; i++ )                             //   Loop over values
	  value[i] = buffer[i];                                 //    Copy back from buffer
#pragma omp for
	for( int i=0; i<size; i++ )                             //   Loop over keys
	  buffer[permutation[i]] = key[i];                      //    Sort into buffer
#pragma omp for
	for( int i=0; i<size; i++ )                             //   Loop over keys
	  key[i] = buffer[i] >> bitStride;                      //    Copy back from buffer and bit shift keys
#pragma omp single
	maxKey >>= bitStride;                                   //   Bit shift maxKey
      }                                                         //  End while for bits to process
    }                                                           // End OpenMP parallel clause
    delete[] bucketPerThread;                                   // Deallocate bucket per thread
    delete[] maxKeyPerThread;                                   // Deallocate maxKey per thread
    delete[] buffer;                                            // Deallocate buffer
    delete[] permutation;                                       // Deallocate permutation index
  }

public:
  //! Sort input accoring to ibody
  Bodies ibody(Bodies & input) {
    const int size = input.size();                              // Size of bodies vector
    int * key = new int [size];                                 // Allocate key array
    int * index = new int [size];                               // Allocate index array
    for (B_iter B=input.begin(); B!=input.end(); B++) {         // Loop over input bodies
      int i = B-input.begin();                                  //  Body index
      key[i] = B->IBODY;                                        //  Copy IBODY to key array
      index[i] = i;                                             //  Initialize index array
    }                                                           // End loop over input bodies
    radixsort(key,index,size);                                  // Radix sort index according to key
    output.resize(size);                                        // Resize output buffer
    for (B_iter B=output.begin(); B!=output.end(); B++) {       // Loop over output boides
      int i = B-output.begin();                                 //  Body index
      *B = input[index[i]];                                     //  Permute according to index
    }                                                           // End loop over output bodies
    delete[] key;                                               // Deallocate key array
    delete[] index;                                             // Deallocate index array
    return output;                                              // Return output
  }

  //! Sort input accoring to irank
  Bodies irank(Bodies & input) {
    const int size = input.size();                              // Size of bodies vector
    int * key = new int [size];                                 // Allocate key array
    int * index = new int [size];                               // Allocate index array
    for (B_iter B=input.begin(); B!=input.end(); B++) {         // Loop over input bodies
      int i = B-input.begin();                                  //  Body index
      key[i] = B->IRANK;                                        //  Copy IRANK to key array
      index[i] = i;                                             //  Initialize index array
    }                                                           // End loop over input bodies
    radixsort(key,index,size);                                  // Radix sort index according to key
    output.resize(size);                                        // Resize output buffer
    for (B_iter B=output.begin(); B!=output.end(); B++) {       // Loop over output boides
      int i = B-output.begin();                                 //  Body index
      *B = input[index[i]];                                     //  Permute according to index
    }                                                           // End loop over output bodies
    delete[] key;                                               // Deallocate key array
    delete[] index;                                             // Deallocate index array
    return output;                                              // Return output
  }

  //! Sort bodies back to original order
  Bodies unsort(Bodies & bodies) {
    bodies = ibody(bodies);                                     // Sort bodies
    return bodies;
  }
};
#endif
