#ifndef tree_h
#define tree_h

class Tree {
  bodies &B;                                                    // Bodies in the tree
  vect X0;                                                      // Center of root cell
  real R0;                                                      // Radius of root cell
public:
  Tree(bodies &b) : B(b),X0(0),R0(0) {}                         // Constructor
  ~Tree() {}                                                    // Destructor

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

  int get_max_level() {
    int const Nbody = B.size();                                 // Number of bodies
    int level;                                                  // Define max level
    level = Nbody >= NCRIT ? 1+log(Nbody/NCRIT)/M_LN2/3 : 0;    // Decide max level from Nbody/Ncrit
    return level;                                               // Return max level
  }

  void get_morton(bigint *index) {
    bigint i;                                                   // Temporary Morton index
    int level = get_max_level();                                // Decide max level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // Define 3-D index
    for( B=B.begin(); B!=B.end(); ++B ) {                       // Loop over all bodies
      i = 0;                                                    //  Initialize Morton index
      for( int d=0; d!=3; ++d )                                 //  Loop over dimension
        nx[d] = int( ( B.pos()[d]-(X0[d]-R0) )/r );             //   3-D index
      for( int l=0; l<level; ++l ) {                            //  Loop over all levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d]%2 << (3*l+d);                              //    Accumulate Morton index
          nx[d] >>= 1;                                          //    Bitshift 3-D index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      index[B] = i;                                             //  Store results
    }                                                           // End loop over all bodies
  }

  template<typename T>
  void sort(bigint *index, T &value, int Nbucket) {
    int const N = value.size();
    int inew,*bucket;
    bigint *ibuffer;
    bucket = new int [Nbucket];
    ibuffer = new bigint [N];
    T vbuffer(N,N);
    for( int i=0; i<Nbucket; ++i ) bucket[i] = 0;
    for( int i=0; i<N; ++i ) bucket[index[i]]++;
    for( int i=1; i<Nbucket; ++i ) bucket[i] += bucket[i-1];
    for( int i=N-1; i>=0; i-- ) {
      bucket[index[i]]--;
      inew = bucket[index[i]];
      ibuffer[inew] = index[i];
      vbuffer[inew] = value[i];
    }
    for( int i=0; i<N; ++i ) index[i] = ibuffer[i];
    for( int i=0; i<N; ++i ) value[i] = vbuffer[i];
    delete[] bucket;
    delete[] ibuffer;
  }

  template<typename T>
  void sortLarge(bigint *index, T &value, int Nbucket) {
    int const N = value.size();
    int inew,*bucket1,*bucket2,*permut1,*permut2;
    bigint *ibuffer;
    bucket1 = new int [Nbucket];
    bucket2 = new int [Nbucket];
    permut1 = new int [N];
    permut2 = new int [N];
    ibuffer = new bigint [N];
    T vbuffer(N,N);
    for( int i=0; i<Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=0; i<N; ++i )
      bucket1[index[i]/Nbucket]++;
    for( int i=1; i<Nbucket; ++i )
      bucket1[i] += bucket1[i-1];
    for( int i=N-1; i>=0; i-- ) {
      bucket1[index[i]/Nbucket]--;
      inew = bucket1[index[i]/Nbucket];
      permut1[inew] = i;
    }
    for( int i=0; i<Nbucket; ++i )
      bucket1[i] = 0;
    for( int i=0; i<N; ++i )
      bucket1[index[permut1[i]]/Nbucket]++;
    int offset(0);
    for( int j=0; j<Nbucket; ++j ) {
      if( bucket1[j] > 0 ) {
        for( int i=0; i<Nbucket; ++i )
          bucket2[i] = 0;
        for( int i=0; i<bucket1[j]; ++i )
          bucket2[index[permut1[i+offset]]%Nbucket]++;
        for( int i=1; i<Nbucket; ++i )
          bucket2[i] += bucket2[i-1];
        for( int i=bucket1[j]-1; i>=0; --i ) {
          bucket2[index[permut1[i+offset]]%Nbucket]--;
          inew = bucket2[index[permut1[i+offset]]%Nbucket];
          permut2[inew+offset] = i+offset;
        }
        offset += bucket1[j];
      }
    }
    for( int i=0; i<N; ++i ) {
      ibuffer[i] = index[permut1[i]];
      vbuffer[i] = value[permut1[i]];
    }
    for( int i=0; i<N; ++i ) {
      index[i] = ibuffer[permut2[i]];
      value[i] = vbuffer[permut2[i]];
    }
    delete[] bucket1;
    delete[] bucket2;
    delete[] permut1;
    delete[] permut2;
    delete[] ibuffer;
  }
};

#endif
