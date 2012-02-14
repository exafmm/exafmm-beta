/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef sort_h
#define sort_h
#include "../include/logger.h"

//! Custom bucket sort for body and cell structures
class Sort : public Logger {
private:
  std::vector<int> bucket;                                      //!< Bucket for sorting

//! Get bucket size for sorting
  template<typename T>
  void getBucketSize(T &values, int begin, int end, bigint &Imin, int &numBucket) {
    typename T::iterator V0 = values.begin()+begin;             // Get begin iterator
    typename T::iterator VN = values.begin()+end;               // Get end iterator
    Imin = V0->ICELL;                                           // Initialize minimum index
    bigint Imax = V0->ICELL;                                    // Initialize maximum index
    for( typename T::iterator V=V0; V!=VN; ++V ) {              // Loop over vector
      if     ( V->ICELL < Imin ) Imin = V->ICELL;               //  Set minimum index
      else if( V->ICELL > Imax ) Imax = V->ICELL;               //  Set maximum index
    }                                                           // End loop over vector
    numBucket = Imax - Imin + 1;                                // Use range of indices as bucket size
    if( numBucket > int(bucket.size()) ) {                      // If bucket size needs to be enlarged
      bucket.resize(numBucket);                                 //  Resize bucket vector
    }                                                           // Endif for resize
  }

//! Bucket sort for small indices
  template<typename T>
  void sortICELL(T &values, T &buffer, bigint Imin,
                 int numBucket, bool ascend, int begin, int end) {
    startTimer("Fill bucket  ");                                // Start timer
    for( int i=0; i!=numBucket; ++i ) bucket[i] = 0;            // Initialize bucket
    for( int i=begin; i!=end; ++i ) bucket[values[i].ICELL-Imin]++;// Fill bucket
    for( int i=1; i!=numBucket; ++i ) bucket[i] += bucket[i-1]; // Scan bucket
    stopTimer("Fill bucket  ");                                 // Stop timer
    startTimer("Empty bucket ");                                // Start timer
    for( int i=end-1; i>=begin; --i ) {                         // Loop over data backwards
      bucket[values[i].ICELL-Imin]--;                           //  Empty bucket
      int inew = bucket[values[i].ICELL-Imin]+begin;            //  Permutation index
      buffer[inew] = values[i];                                 //  Fill buffer
    }                                                           // End loop over data
    stopTimer("Empty bucket ");                                 // Stop timer
    startTimer("Copy value   ");                                // Start timer
    if( ascend ) {                                              // If sorting in ascending order
#pragma omp parallel for
      for( int i=begin; i<end; ++i ) values[i] = buffer[i];     //  Copy back bodiess in order
    } else {                                                    // If sorting in descending order
#pragma omp parallel for
      for( int i=begin; i<end; ++i ) values[end-i+begin-1] = buffer[i];// Copy back bodiess in reverse order
    }                                                           // Endif for sorting order
    stopTimer("Copy value   ");                                 // Stop timer
  }

public:
//! Sort bodies accoring to cell index
  void sortBodies(Bodies &bodies, Bodies &buffer, bool ascend=true, int begin=0, int end=0) {
    startTimer("Sort bodies  ");                                // Start timer
    if( bodies.size() == 0 ) return;                            // Don't do anything if vector is empty
    if( end == 0 ) end = bodies.size();                         // Default range is the whole vector
    int numBucket = 0;                                          // Initialize bucket size
    bigint Imin = 0;                                            // Initialize minimum index
    getBucketSize(bodies,begin,end,Imin,numBucket);             // Get bucket size for sorting
    stopTimer("Sort bodies  ");                                 // Stop timer
    sortICELL(bodies,buffer,Imin,numBucket,ascend,begin,end);   // Call bucket sort for small indices
  }

//! Sort cells according to cell index
  void sortCells(Cells &cells, Cells &buffer, bool ascend=true, int begin=0, int end=0) {
    startTimer("Sort cells   ");                                // Start timer
    if( cells.size() == 0 ) return;                             // Don't do anything if vector is empty
    if( end == 0 ) end = cells.size();                          // Default rage is the whole vector
    int numBucket = 0;                                          // Initialize bucket size
    bigint Imin = 0;                                            // Initialize minimum index
    getBucketSize(cells,begin,end,Imin,numBucket);              // Get bucket size for sorting
    stopTimer("Sort cells   ");                                 // Stop timer
    if( buffer.size() < cells.size() ) buffer.resize(cells.size());// Resize sort buffer if necessary
    sortICELL(cells,buffer,Imin,numBucket,ascend,begin,end);    // Call bucket sort for small indices
  }
};

#endif
