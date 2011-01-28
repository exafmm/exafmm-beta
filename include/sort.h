#ifndef sort_h
#define sort_h

class Sort {
private:
  std::vector<int> bucket,bucket1,bucket2,permut1,permut2;      // Bucket and permutation for sorting
protected:
  std::vector<bigint> ibuffer;                                  // Index buffer for sorting
public:
  template<typename T>
  void sortSmall(Bigints &index, T &value, T &vbuffer, bigint Imin,
                 int numBucket, bool ascend, int begin, int end) {
    int const N = value.size();
    if( end == 0 ) end = N;
    for( int i=0; i!=numBucket; ++i ) bucket[i] = 0;
    for( int i=begin; i!=end; ++i ) bucket[index[i]-Imin]++;
    for( int i=1; i!=numBucket; ++i ) bucket[i] += bucket[i-1];
    for( int i=end-1; i>=begin; --i ) {
      bucket[index[i]-Imin]--;
      int inew = bucket[index[i]-Imin]+begin;
      ibuffer[inew] = index[i];
      vbuffer[inew] = value[i];
    }
    if( ascend ) {
      for( int i=begin; i!=end; ++i ) index[i] = ibuffer[i];
      for( int i=begin; i!=end; ++i ) value[i] = vbuffer[i];
    } else {
      for( int i=begin; i!=end; ++i ) index[end-i+begin-1] = ibuffer[i];
      for( int i=begin; i!=end; ++i ) value[end-i+begin-1] = vbuffer[i];
    }
  }

  template<typename T>
  void sortLarge(Bigints &index, T &value, T &vbuffer, bigint Imin,
                 int numBucket, bool ascend, int begin, int end) {
    int const N = value.size();
    if( end == 0 ) end = N;
    for( int i=0; i!=numBucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(index[i] - Imin) / numBucket]++;
    for( int i=1; i!=numBucket; ++i )
      bucket1[i] += bucket1[i-1];
    for( int i=end-1; i>=begin; --i ) {
      bucket1[(index[i] - Imin) / numBucket]--;
      int inew = bucket1[(index[i] - Imin) / numBucket]+begin;
      permut1[inew] = i;
    }
    for( int i=0; i!=numBucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(index[permut1[i]] - Imin) / numBucket]++;
    int offset(0);
    for( int j=0; j!=numBucket; ++j ) {
      if( bucket1[j] > 0 ) {
        for( int i=0; i!=numBucket; ++i )
          bucket2[i] = 0;
        for( int i=0; i!=bucket1[j]; ++i )
          bucket2[(index[permut1[i+offset]] - Imin) % numBucket]++;
        for( int i=1; i!=numBucket; ++i )
          bucket2[i] += bucket2[i-1];
        for( int i=bucket1[j]-1; i>=0; --i ) {
          bucket2[(index[permut1[i+offset]] - Imin) % numBucket]--;
          int inew = bucket2[(index[permut1[i+offset]] - Imin) % numBucket] + begin;
          permut2[inew+offset] = i+offset;
        }
        offset += bucket1[j];
      }
    }
    for( int i=begin; i!=end; ++i ) {
      ibuffer[i] = index[permut1[i]];
      vbuffer[i] = value[permut1[i]];
    }
    if( ascend ) {
      for( int i=begin; i!=end; ++i ) {
        index[i] = ibuffer[permut2[i]];
        value[i] = vbuffer[permut2[i]];
      }
    } else {
      for( int i=begin; i!=end; ++i ) {
        index[end-i+begin-1] = ibuffer[permut2[i]];
        value[end-i+begin-1] = vbuffer[permut2[i]];
      }
    }
  }

  template<typename T>
  void sort(Bigints &index, T &value, T &vbuffer, bool ascend=true, int begin=0, int end=0) {
    int const N = value.size();
    int const threshold(100000000);
    int numBucket(0);
    if( end == 0 ) end = N;
    BI_iter BI0 = index.begin()+begin;
    BI_iter BIN = index.begin()+end;
    bigint Imin  = *(std::min_element(BI0,BIN));
    bigint Imax  = *(std::max_element(BI0,BIN));
    bigint Isize = Imax - Imin + 1;
    if( N > int(ibuffer.size()) ) {
      ibuffer.resize(N);
    }
    if( Isize < threshold ) {
      numBucket = Isize;
      if( numBucket > int(bucket.size()) ) {
        bucket.resize(numBucket);
      }
      sortSmall(index,value,vbuffer,Imin,numBucket,ascend,begin,end);
    } else {
      numBucket = threshold;
      if( N > int(permut1.size()) ) {
        bucket1.resize(numBucket);
        bucket2.resize(numBucket);
        permut1.resize(N);
        permut2.resize(N);
      }
      sortLarge(index,value,vbuffer,Imin,numBucket,ascend,begin,end);
    }
  }
};

#endif
