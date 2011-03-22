#ifndef sort_h
#define sort_h

class Sort {                                                    // Custom bucket sort for body and cell structures
private:
  std::vector<int> bucket,bucket1,bucket2,permut1,permut2;      // Bucket and permutation for sorting

  template<typename T>
  void sortSmall(T &value, T &vbuffer, bigint Imin,             // Bucket sort for small indices
                 int numBucket, bool ascend, int begin, int end) {
    if( end == 0 ) end = value.size();                          // Default range is the whole vector
    for( int i=0; i!=numBucket; ++i ) bucket[i] = 0;            // Initialize bucket
    for( int i=begin; i!=end; ++i ) bucket[value[i].I-Imin]++;  // Fill bucket
    for( int i=1; i!=numBucket; ++i ) bucket[i] += bucket[i-1]; // Scan bucket
    for( int i=end-1; i>=begin; --i ) {                         // Loop over data backwards
      bucket[value[i].I-Imin]--;                                //  Empty bucket
      int inew = bucket[value[i].I-Imin]+begin;                 //  Permutation index
      vbuffer[inew] = value[i];                                 //  Fill buffer
    }                                                           // End loop over data
    if( ascend ) {                                              // If sorting in ascending order
      for( int i=begin; i!=end; ++i ) value[i] = vbuffer[i];    //  Copy back values in order
    } else {                                                    // If sorting in descending order
      for( int i=begin; i!=end; ++i ) value[end-i+begin-1] = vbuffer[i];// Copy back values in reverse order
    }                                                           // Endif for sorting order
  }

  template<typename T>
  void sortLarge(T &value, T &vbuffer, bigint Imin,             // Bucket sort for large indices
                 int numBucket, bool ascend, int begin, int end) {
    if( end == 0 ) end = value.size();                          // Default range is the whole vector
    for( int i=0; i!=numBucket; ++i )                           // Loop over bucket1
      bucket1[i] = 0;                                           //  Initialize bucket
    for( int i=begin; i!=end; ++i )                             // Loop over data
      bucket1[(value[i].I - Imin) / numBucket]++;               //  Fill bucket1
    for( int i=1; i!=numBucket; ++i )                           // Loop over bucket1
      bucket1[i] += bucket1[i-1];                               //  Scan bucket1
    for( int i=end-1; i>=begin; --i ) {                         // Loop over data backwards
      bucket1[(value[i].I - Imin) / numBucket]--;               //  Empty bucket1
      int inew = bucket1[(value[i].I - Imin) / numBucket]+begin;//  Permuation index
      permut1[inew] = i;                                        //  Store permut1
    }                                                           // End loop over data
    for( int i=0; i!=numBucket; ++i )                           // Loop over bucket1
      bucket1[i] = 0;                                           //  Initialize bucket1
    for( int i=begin; i!=end; ++i )                             // Loop over data
      bucket1[(value[permut1[i]].I - Imin) / numBucket]++;      //  Fill bucket1
    int offset = 0;                                             // Initialize offset
    for( int j=0; j!=numBucket; ++j ) {                         // Loop over bucket1
      if( bucket1[j] > 0 ) {                                    //  If bucket1 is not empty
        for( int i=0; i!=numBucket; ++i )                       //   Loop over bucket2
          bucket2[i] = 0;                                       //    Initialize bucket2
        for( int i=0; i!=bucket1[j]; ++i )                      //   Loop over data
          bucket2[(value[permut1[i+offset]].I - Imin) % numBucket]++;// Fill bucket2
        for( int i=1; i!=numBucket; ++i )                       //   Loop over bucket2
          bucket2[i] += bucket2[i-1];                           //    Scan bucket2
        for( int i=bucket1[j]-1; i>=0; --i ) {                  //   Loop over data backwards
          bucket2[(value[permut1[i+offset]].I - Imin) % numBucket]--;// Empty bucket2
          int inew = bucket2[(value[permut1[i+offset]].I - Imin) % numBucket] + begin;// Permutation index
          permut2[inew+offset] = i+offset;                      //    Store permut2
        }                                                       //   End loop over data
        offset += bucket1[j];                                   //   Increment offset
      }                                                         //  Endif for bucket1
    }                                                           // End loop over bucket1
    for( int i=begin; i!=end; ++i ) {                           // Loop over data
      vbuffer[i] = value[permut1[i]];                           //  Fill buffer
    }                                                           // End loop over data
    if( ascend ) {                                              // If sorting in ascending order
      for( int i=begin; i!=end; ++i ) {                         //  Loop over data
        value[i] = vbuffer[permut2[i]];                         //   Copy back values in order
      }                                                         //  End loop over data
    } else {                                                    // If sorting in descending order
      for( int i=begin; i!=end; ++i ) {                         //  Loop over data
        value[end-i+begin-1] = vbuffer[permut2[i]];             //   Copy back values in reverse order
      }                                                         //  End loop over data
    }                                                           // Endif for sorting order
  }

public:
  template<typename T>
  void sort(T &value, T &vbuffer, bool ascend=true, int begin=0, int end=0) {// Interface for bucket sort
    const int N = value.size();                                 // Set size of sort vector
    if( N == 0 ) return;                                        // Don't do anything if vector is empty
    const int threshold = 100000000;                            // If index is too large use sortLarge()
    int numBucket = 0;                                          // Initialize bucket size
    if( end == 0 ) end = N;                                     // Default range is the whole vector
    typename T::iterator V0 = value.begin()+begin;              // Get begin iterator
    typename T::iterator VN = value.begin()+end;                // Get end iterator
    bigint Imin = V0->I, Imax = V0->I;                          // Initialize min and max of index
    for(typename T::iterator V=V0; V!=VN; ++V ) {               // Loop over data
      if     ( V->I < Imin ) Imin = V->I;                       //  Set minimum index
      else if( V->I > Imax ) Imax = V->I;                       //  Set maximum index
    }                                                           // End loop over data
    bigint Isize = Imax - Imin + 1;                             // Range of indices
    if( Isize < threshold ) {                                   // If range of indices is less than threshold
      numBucket = Isize;                                        //  Use range of indices as bucket size
      if( numBucket > int(bucket.size()) ) {                    //  If bucket size needs to be enlarged
        bucket.resize(numBucket);                               //   Resize bucket vector
      }                                                         //  Endif for resize
      sortSmall(value,vbuffer,Imin,numBucket,ascend,begin,end); //  Call bucket sort for small indices
    } else {                                                    // If range of indices is larger than threshold
      numBucket = threshold;                                    //  Use threshold as bucket size
      if( N > int(permut1.size()) ) {                           //  If sizes of permut and bucket need to be enlarged
        bucket1.resize(numBucket);                              //   Resize bucket1 vector
        bucket2.resize(numBucket);                              //   Resize bucket2 vector
        permut1.resize(N);                                      //   Resize permut1 vector
        permut2.resize(N);                                      //   Resize permut2 vector
      }                                                         //  Endif for resize
      sortLarge(value,vbuffer,Imin,numBucket,ascend,begin,end); //  Cell bucket sort for large indices
    }                                                           // Endif for range of indices
  }

};

#endif
