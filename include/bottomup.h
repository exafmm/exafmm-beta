#ifndef bottomup_h
#define bottomup_h
#include "tree.h"
#include "sort.h"

class BottomUpTreeConstructor : public TreeStructure, public Sort {
public:
  BottomUpTreeConstructor(bodies &b) : TreeStructure(b){}       // Constructor
  ~BottomUpTreeConstructor() {}                                 // Destructor

  int getMaxLevel() {
    int const Nbody = B.size();                                 // Number of bodies
    int level;                                                  // Define max level
    level = Nbody >= NCRIT ? 1+log(Nbody/NCRIT)/M_LN2/3 : 0;    // Decide max level from Nbody/Ncrit
    return level;                                               // Return max level
  }

  void setMorton(bigint *index, int level=0, int begin=0, int end=0 ) {
    bigint i;                                                   // Levelwise Morton index
    if( level == 0 ) level = getMaxLevel();                     // Decide max level
    bigint offset = ((1 << 3*level) - 1)/7;                     // Offset for each level
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
      index[B] = i+offset;                                      //  Store results with levelwise offset
    }                                                           // End loop over all bodies
  }

  int getMortonLevel(bigint index) {
    int level(-1);
    while( index >= 0 ) {
      level++;
      index -= 1 << 3*level;
    }
    return level;
  }

  void prune(bigint *index) {
    int icell,begin,size;
    int maxLevel = getMaxLevel();
    int level;
    bigint cOffset,pOffset,pIndex;
    for( int l=maxLevel; l>0; --l ) {
      level = getMortonLevel(index[0]);
      cOffset = ((1 << 3*level) - 1)/7;
      pOffset = ((1 << 3*(l-1)) - 1)/7;
      icell = ((index[0]-cOffset) >> 3*(level-l+1)) + pOffset;
      begin = 0;
      size = 0;
      for( B=B.begin(); B!=B.end(); ++B ) {
        level = getMortonLevel(index[B]);
        cOffset = ((1 << 3*level) - 1)/7;
        pIndex = ((index[B]-cOffset) >> 3*(level-l+1)) + pOffset;
        if( pIndex != icell ) {
          if( size <= NCRIT ) {
            for( int i=begin; i!=begin+size; ++i ) {
              index[i] = icell;
            }
          }
          begin = B;
          size = 0;
          icell = pIndex;
        }
        size++;
      }
      if( size <= NCRIT ) {
        for( int i=begin; i!=begin+size; ++i ) {
          index[i] = icell;
        }
      }
    }
  }

  void grow(bigint *index, bodies &buffer, int level=0, int begin=0, int end=0) {
    int icell(index[begin]),offset(begin),size(0);
    if( level == 0 ) level = getMaxLevel();
    if( end == 0 ) end = B.end();
    for( B=begin; B!=end; ++B ) {
      if( index[B] != icell ) {
        if( size > NCRIT ) {
          level++;
          setMorton(index,level,offset,offset+size);
          sort(index,B,buffer,true,offset,offset+size);
          grow(index,buffer,level,offset,offset+size);
          level--;
        }
        offset = B;
        size = 0;
        icell = index[B];
      }
      size++;
    }
    if( size > NCRIT ) {
      level++;
      setMorton(index,level,offset,offset+size);
      sort(index,B,buffer,true,offset,offset+size);
      grow(index,buffer,level,offset,offset+size);
      level--;
    }
  }

  void link(bigint *index) {
    CN = C0;
    int icell(index[0]),begin(0),size(0),Ncell(0);
    for( B=B.begin(); B!=B.end(); ++B ) {
      if( index[B] != icell ) {
        CN->NLEAF = size;
        CN->NCHILD = 0;
        CN->I = icell;
        CN->LEAF = &B[begin];
        ++CN;
        begin = B;
        size = 0;
        icell = index[B];
        Ncell++;
      }
      size++;
    }
    CN->NLEAF = size;
    CN->NCHILD = 0;
    CN->I = icell;
    CN->LEAF = &B[begin];
    ++CN;

    for( C=C0; C!=CN; ++C ) {
      std::cout << C->I << std::endl;
    }
  }
};

#endif
