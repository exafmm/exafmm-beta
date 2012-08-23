#ifndef sort_h
#define sort_h
#include "types.h"

//! Custom bucket sort for body and cell structures
class Sort {
private:
  std::vector<int> bucket;                                      //!< Bucket for sorting

//! Get bucket size for sorting
  template<typename T>
  void getBucketSize(T &values, int begin, int end, int &Imin, int &numBucket) {
    typename T::iterator V0 = values.begin()+begin;             // Get begin iterator
    typename T::iterator VN = values.begin()+end;               // Get end iterator
    Imin = V0->ICELL;                                           // Initialize minimum index
    int Imax = V0->ICELL;                                       // Initialize maximum index
    for (typename T::iterator V=V0; V!=VN; V++) {               // Loop over vector
      if      (V->ICELL < Imin) Imin = V->ICELL;                //  Set minimum index
      else if (V->ICELL > Imax) Imax = V->ICELL;                //  Set maximum index
    }                                                           // End loop over vector
    numBucket = Imax - Imin + 1;                                // Use range of indices as bucket size
    if( numBucket > int(bucket.size()) ) {                      // If bucket size needs to be enlarged
      bucket.resize(numBucket);                                 //  Resize bucket vector
    }                                                           // Endif for resize
  }

//! Bucket sort for small indices
  template<typename T>
  void sortICELL(T &values, T &buffer, int Imin,
                 int numBucket, bool ascend, int begin, int end) {
    for (int i=0; i<numBucket; i++) bucket[i] = 0;              // Initialize bucket
    for (int i=begin; i<end; i++) bucket[values[i].ICELL-Imin]++;// Fill bucket
    for (int i=1; i<numBucket; i++) bucket[i] += bucket[i-1];   // Scan bucket
    for (int i=end-1; i>=begin; i--) {                          // Loop over data backwards
      bucket[values[i].ICELL-Imin]--;                           //  Empty bucket
      int inew = bucket[values[i].ICELL-Imin]+begin;            //  Permutation index
      buffer[inew] = values[i];                                 //  Fill buffer
    }                                                           // End loop over data
    if (ascend) {                                               // If sorting in ascending order
      for (int i=begin; i<end; i++) values[i] = buffer[i];      //  Copy back bodiess in order
    } else {                                                    // If sorting in descending order
      for (int i=begin; i<end; i++) values[end-i+begin-1] = buffer[i];// Copy back bodiess in reverse order
    }                                                           // Endif for sorting order
  }

public:
//! Sort bodies accoring to cell index
  void sortBodies(Bodies &bodies, Bodies &buffer, bool ascend=true, int begin=0, int end=0) {
    if (bodies.size() == 0) return;                             // Don't do anything if vector is empty
    if (end == 0) end = bodies.size();                          // Default range is the whole vector
    int numBucket = 0;                                          // Initialize bucket size
    int Imin = 0;                                               // Initialize minimum index
    getBucketSize(bodies,begin,end,Imin,numBucket);             // Get bucket size for sorting
    sortICELL(bodies,buffer,Imin,numBucket,ascend,begin,end);   // Call bucket sort for small indices
  }
};

#endif
