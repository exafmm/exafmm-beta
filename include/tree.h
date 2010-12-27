#ifndef tree_h
#define tree_h

class Tree {
  bodies &B;                                                    // Bodies in the tree
  vect X0;                                                      // Center of root cell
  real R0;                                                      // Radius of root cell
  int *bucket,*bucket1,*bucket2,*permut1,*permut2;              // Bucket and permutation for sorting
  bigint *ibuffer;                                              // Index buffer for sorting
  bool sortAlloc,smallAlloc,largeAlloc;                         // Are sort variables allocated?
public:
  Tree(bodies &b) : B(b),X0(0),R0(0),sortAlloc(false),
                    smallAlloc(false),largeAlloc(false) {}      // Constructor
  ~Tree() {}                                                    // Destructor

  vect getX0() {return X0;}
  real getR0() {return R0;}

  void setDomain() {
    vect xmin,xmax;                                             // Min,Max of domain
    B.begin();                                                  // Reset bodies counter
    xmin = xmax = B.pos();                                      // Initialize xmin,xmax
    for( B=B.begin(); B!=B.end(); ++B ) {                       // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B.pos()[d] < xmin[d]) xmin[d] = B.pos()[d];     //   Determine xmin
        else if(B.pos()[d] < xmax[d]) xmax[d] = B.pos()[d];     //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B.pos();                                            //  Sum positions
    }                                                           // End loop over all bodies
    X0 /= B.size();                                             // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d]-X0[d],R0);                          //  Calculate max distance from center
      R0 = std::max(X0[d]-xmin[d],R0);                          //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1.+log(R0)/M_LN2));                         // Add some leeway to root radius
  }

  int getMaxLevel() {
    int const Nbody = B.size();                                 // Number of bodies
    int level;                                                  // Define max level
    level = Nbody >= NCRIT ? 1+log(Nbody/NCRIT)/M_LN2/3 : 0;    // Decide max level from Nbody/Ncrit
    return level;                                               // Return max level
  }

  void getMorton(bigint *index, int level=0, int begin=0, int end=0 ) {
    bigint i;                                                   // Temporary Morton index
    if( level == 0 ) level = getMaxLevel();                     // Decide max level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // Define 3-D index
    if( end == 0 ) end = B.end();
    for( B=begin; B!=end; ++B ) {                               // Loop over all bodies
      i = 0;                                                    //  Initialize Morton index
      for( int d=0; d!=3; ++d )                                 //  Loop over dimension
        nx[d] = int( ( B.pos()[d]-(X0[d]-R0) )/r );             //   3-D index
      for( int l=0; l!=level; ++l ) {                            //  Loop over all levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d]%2 << (3*l+d);                              //    Accumulate Morton index
          nx[d] >>= 1;                                          //    Bitshift 3-D index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      index[B] = i;                                             //  Store results
    }                                                           // End loop over all bodies
  }

  template<typename T>
  void sortSmall(bigint *index, T &value, T &vbuffer, int Nbucket, int begin=0, int end=0) {
    int const N = value.size();
    int inew;
    if( end == 0 ) end = N;
    for( int i=0; i!=Nbucket; ++i ) bucket[i] = 0;
    for( int i=begin; i!=end; ++i ) bucket[index[i]]++;
    for( int i=1; i!=Nbucket; ++i ) bucket[i] += bucket[i-1];
    for( int i=end-1; i>=begin; i-- ) {
      bucket[index[i]]--;
      inew = bucket[index[i]]+begin;
      ibuffer[inew] = index[i];
      vbuffer[inew] = value[i];
    }
    for( int i=begin; i!=end; ++i ) index[i] = ibuffer[i];
    for( int i=begin; i!=end; ++i ) value[i] = vbuffer[i];
  }

  template<typename T>
  void sortLarge(bigint *index, T &value, T &vbuffer, int Nbucket, int begin=0, int end=0) {
    int const N = value.size();
    int inew;
    if( end == 0 ) end = N;
    for( int i=0; i!=Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[index[i]/Nbucket]++;
    for( int i=1; i!=Nbucket; ++i )
      bucket1[i] += bucket1[i-1];
    for( int i=end-1; i>=begin; i-- ) {
      bucket1[index[i]/Nbucket]--;
      inew = bucket1[index[i]/Nbucket]+begin;
      permut1[inew] = i;
    }
    for( int i=0; i!=Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=begin; i!=end; ++i )
      bucket1[index[permut1[i]]/Nbucket]++;
    int offset(0);
    for( int j=0; j!=Nbucket; ++j ) {
      if( bucket1[j] > 0 ) {
        for( int i=0; i!=Nbucket; ++i )
          bucket2[i] = 0;
        for( int i=0; i!=bucket1[j]; ++i )
          bucket2[index[permut1[i+offset]]%Nbucket]++;
        for( int i=1; i!=Nbucket; ++i )
          bucket2[i] += bucket2[i-1];
        for( int i=bucket1[j]-1; i>=0; --i ) {
          bucket2[index[permut1[i+offset]]%Nbucket]--;
          inew = bucket2[index[permut1[i+offset]]%Nbucket]+begin;
          permut2[inew+offset] = i+offset;
        }
        offset += bucket1[j];
      }
    }
    for( int i=begin; i!=end; ++i ) {
      ibuffer[i] = index[permut1[i]];
      vbuffer[i] = value[permut1[i]];
    }
    for( int i=begin; i!=end; ++i ) {
      index[i] = ibuffer[permut2[i]];
      value[i] = vbuffer[permut2[i]];
    }
  }

  template<typename T>
  void sort(bigint *index, T &value, T &vbuffer, int begin=0, int end=0) {
    int const N = value.size();
    int const threshold = 100000000;
    int Nbucket(0);
    bigint maxIndex = index[0];
    if( end == 0 ) end = N;
    for( int i=begin; i!=end; ++i ) maxIndex = std::max(maxIndex,index[i]);
    if( !sortAlloc ) {
      ibuffer = new bigint [N];
      sortAlloc = true;
    }
    if( maxIndex < threshold ) {
      Nbucket = maxIndex+1;
      if( !smallAlloc ) {
        bucket = new int [Nbucket];
        smallAlloc = true;
      }
      sortSmall(index,value,vbuffer,Nbucket,begin,end);
    } else {
      Nbucket = threshold;
      if( !largeAlloc ) {
        bucket1 = new int [Nbucket];
        bucket2 = new int [Nbucket];
        permut1 = new int [N];
        permut2 = new int [N];
        largeAlloc = true;
      }
      sortLarge(index,value,vbuffer,Nbucket,begin,end);
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

  void grow(bigint *index, bodies &buffer, int level=0, int begin=0, int end=0) {
    int icell(index[0]),offset(0),size(0),Ncell(1);
    if( level == 0 ) level = getMaxLevel();
    if( end == 0 ) end = B.end();
    for( B=begin; B!=end; ++B ) {
      if( index[B] != icell ) {
        if( size > NCRIT ) {
          level++;
          getMorton(index,level,offset,offset+size);
          sort(index,B,buffer,offset,offset+size);
          grow(index,buffer,level,offset,offset+size);
          level--;
        }
        offset = B;
        size = 0;
        icell = index[B];
        Ncell++;
      }
      size++;
    }
    if( size > NCRIT ) {
      level++;
      getMorton(index,level,offset,offset+size);
      sort(index,B,buffer,offset,offset+size);
      grow(index,buffer,level,offset,offset+size);
      level--;
    }
  }

};

#endif
