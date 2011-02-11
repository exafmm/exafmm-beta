#ifndef sort_h
#define sort_h

class Sort {
private:
  std::vector<int> bucket,bucket1,bucket2,permut1,permut2;      // Bucket and permutation for sorting
public:
  template<typename T>
  void sortSmall(T &value, T &vbuffer, bigint Imin,
                 int numBucket, bool ascend, int begin, int end) {
    int const N = value.size();
    if( end == 0 ) end = N;
    for( int i=0; i!=numBucket; ++i ) bucket[i] = 0;
    for( int i=begin; i!=end; ++i ) bucket[value[i].I-Imin]++;
    for( int i=1; i!=numBucket; ++i ) bucket[i] += bucket[i-1];
    for( int i=end-1; i>=begin; --i ) {
      bucket[value[i].I-Imin]--;
      int inew = bucket[value[i].I-Imin]+begin;
      vbuffer[inew] = value[i];
    }
    if( ascend ) {
      for( int i=begin; i!=end; ++i ) value[i] = vbuffer[i];
    } else {
      for( int i=begin; i!=end; ++i ) value[end-i+begin-1] = vbuffer[i];
    }
  }

  template<typename T>
  void sortLarge(T &value, T &vbuffer, bigint Imin,
                 int numBucket, bool ascend, int begin, int end) {
    int const N = value.size();
    if( end == 0 ) end = N;
    for( int i=0; i!=numBucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(value[i].I - Imin) / numBucket]++;
    for( int i=1; i!=numBucket; ++i )
      bucket1[i] += bucket1[i-1];
    for( int i=end-1; i>=begin; --i ) {
      bucket1[(value[i].I - Imin) / numBucket]--;
      int inew = bucket1[(value[i].I - Imin) / numBucket]+begin;
      permut1[inew] = i;
    }
    for( int i=0; i!=numBucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(value[permut1[i]].I - Imin) / numBucket]++;
    int offset = 0;
    for( int j=0; j!=numBucket; ++j ) {
      if( bucket1[j] > 0 ) {
        for( int i=0; i!=numBucket; ++i )
          bucket2[i] = 0;
        for( int i=0; i!=bucket1[j]; ++i )
          bucket2[(value[permut1[i+offset]].I - Imin) % numBucket]++;
        for( int i=1; i!=numBucket; ++i )
          bucket2[i] += bucket2[i-1];
        for( int i=bucket1[j]-1; i>=0; --i ) {
          bucket2[(value[permut1[i+offset]].I - Imin) % numBucket]--;
          int inew = bucket2[(value[permut1[i+offset]].I - Imin) % numBucket] + begin;
          permut2[inew+offset] = i+offset;
        }
        offset += bucket1[j];
      }
    }
    for( int i=begin; i!=end; ++i ) {
      vbuffer[i] = value[permut1[i]];
    }
    if( ascend ) {
      for( int i=begin; i!=end; ++i ) {
        value[i] = vbuffer[permut2[i]];
      }
    } else {
      for( int i=begin; i!=end; ++i ) {
        value[end-i+begin-1] = vbuffer[permut2[i]];
      }
    }
  }

  template<typename T>
  void sort(T &value, T &vbuffer, bool ascend=true, int begin=0, int end=0) {
    int const N = value.size();
    if( N == 0 ) return;
    int const threshold = 100000000;
    int numBucket = 0;
    if( end == 0 ) end = N;
    typename T::iterator V0 = value.begin()+begin;
    typename T::iterator VN = value.begin()+end;
    bigint Imin = V0->I, Imax = V0->I;
    for(typename T::iterator V=V0; V!=VN; ++V ) {
      if     ( V->I < Imin ) Imin = V->I;
      else if( V->I > Imax ) Imax = V->I;
    }
    bigint Isize = Imax - Imin + 1;
    if( Isize < threshold ) {
      numBucket = Isize;
      if( numBucket > int(bucket.size()) ) {
        bucket.resize(numBucket);
      }
      sortSmall(value,vbuffer,Imin,numBucket,ascend,begin,end);
    } else {
      numBucket = threshold;
      if( N > int(permut1.size()) ) {
        bucket1.resize(numBucket);
        bucket2.resize(numBucket);
        permut1.resize(N);
        permut2.resize(N);
      }
      sortLarge(value,vbuffer,Imin,numBucket,ascend,begin,end);
    }
  }
};

#endif
