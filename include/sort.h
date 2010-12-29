#ifndef sort_h
#define sort_h

class Sort {
protected:
  int sortAlloc,smallAlloc,largeAlloc;                          // Size of allocation for sort variables
  int *bucket,*bucket1,*bucket2,*permut1,*permut2;              // Bucket and permutation for sorting
  bigint *ibuffer;                                              // Index buffer for sorting

public:
  Sort() : sortAlloc(0),smallAlloc(0),largeAlloc(0) {}          // Constructor
  ~Sort() {}                                                    // Destructor
  template<typename T>
  void sortSmall(bigint *index, T &value, T &vbuffer, bigint Imin,
                 int Nbucket, bool ascend, int begin, int end) {
    int const N = value.size();
    int inew;
    if( end == 0 ) end = N;
    for( int i=0; i!=Nbucket; ++i ) bucket[i] = 0;
    for( int i=begin; i!=end; ++i ) bucket[index[i]-Imin]++;
    for( int i=1; i!=Nbucket; ++i ) bucket[i] += bucket[i-1];
    for( int i=end-1; i>=begin; i-- ) {
      bucket[index[i]-Imin]--;
      inew = bucket[index[i]-Imin]+begin;
      ibuffer[inew] = index[i];
      vbuffer[inew] = value[i];
    }
    if( ascend ) {
      for( int i=begin; i!=end; ++i ) index[i] = ibuffer[i];
      for( int i=begin; i!=end; ++i ) value[i] = vbuffer[i];
    } else {
      for( int i=begin; i!=end; ++i ) index[end-i-1] = ibuffer[i];
      for( int i=begin; i!=end; ++i ) value[end-i-1] = vbuffer[i];
    }
  }

  template<typename T>
  void sortLarge(bigint *index, T &value, T &vbuffer, bigint Imin,
                 int Nbucket, bool ascend, int begin, int end) {
    int const N = value.size();
    int inew;
    if( end == 0 ) end = N;
    for( int i=0; i!=Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(index[i]-Imin)/Nbucket]++;
    for( int i=1; i!=Nbucket; ++i )
      bucket1[i] += bucket1[i-1];
    for( int i=end-1; i>=begin; i-- ) {
      bucket1[(index[i]-Imin)/Nbucket]--;
      inew = bucket1[(index[i]-Imin)/Nbucket]+begin;
      permut1[inew] = i;
    }
    for( int i=0; i!=Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[(index[permut1[i]]-Imin)/Nbucket]++;
    int offset(0);
    for( int j=0; j!=Nbucket; ++j ) {
      if( bucket1[j] > 0 ) {
        for( int i=0; i!=Nbucket; ++i )
          bucket2[i] = 0;
        for( int i=0; i!=bucket1[j]; ++i )
          bucket2[(index[permut1[i+offset]]-Imin)%Nbucket]++;
        for( int i=1; i!=Nbucket; ++i )
          bucket2[i] += bucket2[i-1];
        for( int i=bucket1[j]-1; i>=0; --i ) {
          bucket2[(index[permut1[i+offset]]-Imin)%Nbucket]--;
          inew = bucket2[(index[permut1[i+offset]]-Imin)%Nbucket]+begin;
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
        index[end-i-1] = ibuffer[permut2[i]];
        value[end-i-1] = vbuffer[permut2[i]];
      }
    }
  }

  template<typename T>
  void sort(bigint *index, T &value, T &vbuffer, bool ascend=true, int begin=0, int end=0) {
    int const N = value.size();
    int const threshold = 100000000;
    int Nbucket(0);
    bigint Imin(index[begin]),Imax(index[begin]),Isize(0);
    if( end == 0 ) end = N;
    for( int i=begin; i!=end; ++i ) {
      if     (index[i] < Imin) Imin = index[i];
      else if(index[i] > Imax) Imax = index[i];
    }
    Isize = Imax-Imin+1;
    if( N > sortAlloc ) {
      if( sortAlloc != 0 ) delete[] ibuffer;
      ibuffer = new bigint [N];
      sortAlloc = N;
    }
    if( Isize < threshold ) {
      Nbucket = Isize;
      if( Nbucket > smallAlloc ) {
        if( smallAlloc != 0 ) delete[] bucket;
        smallAlloc = 1 << int(log(1.+Nbucket)/M_LN2/3+1)*3;
        bucket = new int [smallAlloc];
      }
      sortSmall(index,value,vbuffer,Imin,Nbucket,ascend,begin,end);
    } else {
      Nbucket = threshold;
      if( N > largeAlloc ) {
        if( largeAlloc != 0 ) {
          delete[] bucket1;
          delete[] bucket2;
          delete[] permut1;
          delete[] permut2;
        }
        bucket1 = new int [Nbucket];
        bucket2 = new int [Nbucket];
        permut1 = new int [N];
        permut2 = new int [N];
        largeAlloc = N;
      }
      sortLarge(index,value,vbuffer,Imin,Nbucket,ascend,begin,end);
    }
  }

  void sortDealloc() {
    if( sortAlloc ) delete[] ibuffer;
    if( smallAlloc ) delete[] bucket;
    if( largeAlloc ) {
      delete[] bucket1;
      delete[] bucket2;
      delete[] permut1;
      delete[] permut2;
    }
  }
};

#endif
