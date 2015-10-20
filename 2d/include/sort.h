#ifndef sort_h
#define sort_h
#include "types.h"

//! Custom bucket sort for body and structures
class Sort {
 private:
  typedef Bodies::reverse_iterator B_ritr;                      //!< Reverse iterator for Bodies
  std::vector<int> bucket;                                      //!< Bucket for sorting
  Bodies output;                                                //!< Output buffer

 public:
//! Sort input accoring to cell index
  Bodies sortBodies(Bodies &input) {
    if (input.size() == 0) return input;                        // Return if input is empty
    int Imin = input[0].ICELL;                                  // Initialize minimum index
    int Imax = input[0].ICELL;                                  // Initialize maximum index
    for (B_iter B=input.begin(); B!=input.end(); B++) {         // Loop over vector
      if      (B->ICELL < Imin) Imin = B->ICELL;                //  Set minimum index
      else if (B->ICELL > Imax) Imax = B->ICELL;                //  Set maximum index
    }                                                           // End loop over vector
    int numBucket = Imax - Imin + 1;                            // Use range of indices as bucket size
    if( numBucket > int(bucket.size()) ) {                      // If bucket size needs to be enlarged
      bucket.resize(numBucket);                                 //  Resize bucket vector
    }                                                           // End if for resize
    if( input.size() > output.capacity() ) {                    // If buffer size need to be enlarged
      output.reserve(input.size());                             //  Resize output buffer
    }                                                           // End if for resize
    output = input;                                             // Resize output
    for (int i=0; i<numBucket; i++) bucket[i] = 0;              // Initialize bucket
    for (B_iter B=input.begin(); B!=input.end(); B++) bucket[B->ICELL-Imin]++;// Fill bucket
    for (int i=1; i<numBucket; i++) bucket[i] += bucket[i-1];   // Scan bucket
    for (B_ritr B=input.rbegin(); B!=input.rend(); B++) {       // Loop over data backwards
      bucket[B->ICELL-Imin]--;                                  //  Empty bucket
      int inew = bucket[B->ICELL-Imin];                         //  Permutation index
      output[inew] = *B;                                        //  Fill output
    }                                                           // End loop over data
    return output;                                              // Return output
  }
};
#endif
