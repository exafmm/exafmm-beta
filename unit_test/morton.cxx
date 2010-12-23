#include "body.h"

struct node {
  int NLEAF;
  unsigned NCHILD;
  u64 I;
  bodies *LEAF;

  void init() {
    NLEAF = 0;
    NCHILD = 0;
  }
};

int get_max_level(bodies &B) {
  int const Nbody = B.size();                                   // Number of bodies
  int level;                                                    // Define max level
  level = Nbody >= NCRIT ? 1+log(Nbody/NCRIT)/M_LN2/3 : 0;      // Decide max level from Nbody/Ncrit
  return level;                                                 // Return max level
}

void get_morton(bodies &B, vect x0, real r0, u64 *index) {
  u64 i;                                                        // Temporary Morton index
  int level = get_max_level(B);                                 // Decide max level
  float r = r0 / (1 << (level-1));                              // Radius at finest level
  vec<3,int> nx;                                                // Define 3-D index
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    i = 0;                                                      //  Initialize Morton index
    for( int d=0; d!=3; ++d )                                   //  Loop over dimension
      nx[d] = int( ( B.pos()[d]-(x0[d]-r0) )/r );               //   3-D index
    for( int l=0; l<level; ++l ) {                              //  Loop over all levels of tree
      for( int d=0; d!=3; ++d ) {                               //   Loop over dimension
        i += nx[d]%2 << (3*l+d);                                //    Accumulate Morton index
        nx[d] >>= 1;                                            //    Bitshift 3-D index
      }                                                         //   End loop over dimension
    }                                                           //  End loop over levels
    index[B] = i;                                               //  Store results
  }                                                             // End loop over all bodies
}

vec<3,int> un_morton(u64 i) {
  int d(0),l(0);                                                // Define dimension and level counter
  vec<3,int> nx(0);                                             // Define 3-D index
  while( i != 0 ) {                                             // Loop while working index is nozero
    nx[d] += (i%2)*(1 << l);                                    //  Accumulate 3-D index
    i >>= 1;                                                    //  Bitshift working index
    d = (d+1)%3;                                                //  Update dimension counter
    if( d == 0 ) ++l;                                           //  Update level counter
  }                                                             // End while loop
  return nx;                                                    // Return 3-D index
}

void sort(u64 *index, bodies &B) {
  int const Nbody = B.size();
  int const level = get_max_level(B);
  int const Nbucket = 1 << 3*level;
  unsigned inew,*bucket;
  bucket = new unsigned [Nbucket];
  bodies buffer(Nbody,Nbody);
  for( int i=0; i<Nbucket; ++i ) bucket[i] = 0;
  for( int i=0; i<int(Nbody); ++i ) bucket[index[i]]++;
  for( int i=1; i<Nbucket; ++i ) bucket[i] += bucket[i-1];
  for( int i=Nbody-1; i>=0; i-- ) {
    bucket[index[i]]--;
    inew = bucket[index[i]];
    buffer[inew] = B[i];
  }
  for( int i=0; i<int(Nbody); ++i ) B[i] = buffer[i];
  delete[] bucket;
}

void sort_large(u64 *index, bodies &B) {
  int const Nbody = B.size();
  int const level = get_max_level(B);
  int const digit = floor(log10(1 << 3*level)/2)+1;
  int const Nbucket = pow(10,digit);
  unsigned inew,*bucket1,*bucket2,*permut1,*permut2;
  bucket1 = new unsigned [Nbucket];
  bucket2 = new unsigned [Nbucket];
  permut1 = new unsigned [Nbody];
  permut2 = new unsigned [Nbody];
  bodies buffer(Nbody,Nbody);
  for( int i=0; i<Nbucket; ++i )
    bucket1[i] = 0;
  for( int i=0; i<int(Nbody); ++i )
    bucket1[index[i]/Nbucket]++;
  for( int i=1; i<Nbucket; ++i )
    bucket1[i] += bucket1[i-1];
  for( int i=Nbody-1; i>=0; i-- ) {
    bucket1[index[i]/Nbucket]--;
    inew = bucket1[index[i]/Nbucket];
    permut1[inew] = i;
  }
  for( int i=0; i<Nbucket; ++i )
    bucket1[i] = 0;
  for( int i=0; i<int(Nbody); ++i )
    bucket1[index[permut1[i]]/Nbucket]++;
  int offset(0);
  for( int j=0; j<Nbucket; ++j ) {
    if( bucket1[j] > 0 ) {
      for( int i=0; i<Nbucket; ++i )
        bucket2[i] = 0;
      for( int i=0; i<int(bucket1[j]); ++i )
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
  for( int i=0; i<int(Nbody); ++i ) buffer[i] = B[permut1[i]];
  for( int i=0; i<int(Nbody); ++i ) B[i] = buffer[permut2[i]];
  delete[] bucket1;
  delete[] bucket2;
  delete[] permut1;
  delete[] permut2;
}

int main()
{
  double tic,toc;
  unsigned const Nbody=50;
  bodies B(Nbody,Nbody);

  tic = get_time();
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] = rand()/(1.+RAND_MAX)*2-1;                    //   Initialize positions
    B.scal() = 1./B.size();                                     //  Initialize source value
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  real r0(0);                                                   // Root radius
  vect xmin,xmax,x0(0);                                         // Min,Max,Center of domain
  B.begin();                                                    // Reset bodies counter
  xmin = xmax = B.pos();                                        // Initialize xmin,xmax
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over each dimension
      if     (B.pos()[d] < xmin[d]) xmin[d] = B.pos()[d];       //   Determine xmin
      else if(B.pos()[d] < xmax[d]) xmax[d] = B.pos()[d];       //   Determine xmax
    }                                                           //  End loop over each dimension
    x0 += B.pos();                                              //  Sum positions
  }                                                             // End loop over all bodies
  x0 /= B.size();                                               // Calculate average position
  for( int d=0; d!=3; ++d ) {                                   // Loop over each dimension
    x0[d] = int(x0[d]+.5);                                      //  Shift center to nearest integer
    r0 = std::max(xmax[d]-x0[d],r0);                            //  Calculate max distance from center
    r0 = std::max(x0[d]-xmin[d],r0);                            //  Calculate max distance from center
  }                                                             // End loop over each dimension
  r0 = pow(2.,int(1.+log(r0)/M_LN2));                           // Add some leeway to root radius
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  u64 *index;
  node *N;
  index = new u64 [Nbody];
  get_morton(B,x0,r0,index);
  sort(index,B);
  toc = get_time();
  std::cout << "Construct tree: " << toc-tic << std::endl;
  get_morton(B,x0,r0,index);
  for( B=B.begin(); B!=B.end(); ++B ) {
    std::cout << B << " " << index[B] << std::endl;
  }
  N = new node [Nbody];
  delete[] N;
  delete[] index;
}
