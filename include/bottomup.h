#ifndef bottomup_h
#define bottomup_h
#include "tree.h"
#include "sort.h"

class BottomUpTreeConstructor : public TreeStructure, public Sort {
public:
  BottomUpTreeConstructor(Bodies &b) : TreeStructure(b){}       // Constructor
  ~BottomUpTreeConstructor() {}                                 // Destructor

  int getMaxLevel() {
    int const N = bodies.size();                                // Number of bodies
    int level;                                                  // Define max level
    level = N >= NCRIT ? 1+log(N/NCRIT)/M_LN2/3 : 0;            // Decide max level from N/Ncrit
    return level;                                               // Return max level
  }

  void setMorton(int level=0, int begin=0, int end=0 ) {
    bigint i;                                                   // Levelwise Morton index
    if( level == 0 ) level = getMaxLevel();                     // Decide max level
    bigint offset = ((1 << 3*level) - 1)/7;                     // Offset for each level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // Define 3-D index
    if( end == 0 ) end = bodies.size();
    for( int b=begin; b!=end; ++b ) {                           // Loop over all bodies
      for( int d=0; d!=3; ++d )                                 //  Loop over dimension
        nx[d] = int( ( bodies.at(b).pos[d]-(X0[d]-R0) )/r );    //   3-D index
      i = 0;                                                    //  Initialize Morton index
      for( int l=0; l!=level; ++l ) {                           //  Loop over all levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d]%2 << (3*l+d);                              //    Accumulate Morton index
          nx[d] >>= 1;                                          //    Bitshift 3-D index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      Ibody[b] = i+offset;                                      //  Store results with levelwise offset
    }                                                           // End loop over all bodies
  }

  int getLevel(bigint index) {
    int level(-1);
    while( index >= 0 ) {
      level++;
      index -= 1 << 3*level;
    }
    return level;
  }

  bigint getParent(bigint index) {
    int level = getLevel(index);
    int cOffset = ((1 << 3*level) - 1)/7;
    int pOffset = ((1 << 3*(level-1)) - 1)/7;
    bigint i = ((index-cOffset) >> 3) + pOffset;
    return i;
  }

  void sortCells(Cells &buffer) {
    int begin = CC0-C0;
    int end = CCN-C0;
    int c = begin;
    for( C=CC0; C!=CCN; ++C,++c ) Icell[c] = C->I;
    sort(Icell,cells,buffer,false,begin,end);
  }

  void linkParent(Cells &buffer) {
    CCN = CN;
    sortCells(buffer);
    CN->I = getParent(CC0->I);
    CN->LEAF = CC0->LEAF;
    for( C=CC0; C!=CCN; ++C ) {
      if( getParent(C->I) != CN->I ) {
        ++CN;
        CN->I = getParent(C->I);
        CN->LEAF = C->LEAF;
      }
      C->PARENT = CN;
      CN->NLEAF += C->NLEAF;
      CN->CHILD[CN->NCHILD] = C;
      CN->NCHILD++;
    }
    ++CN;
    CC0 = CCN;
  }

  void prune() {
    int icell,begin,size;
    int maxLevel = getMaxLevel();
    int level;
    bigint cOffset,pOffset,pIndex;
    for( int l=maxLevel; l>0; --l ) {
      level = getLevel(Ibody[0]);
      cOffset = ((1 << 3*level) - 1)/7;
      pOffset = ((1 << 3*(l-1)) - 1)/7;
      icell = ((Ibody[0]-cOffset) >> 3*(level-l+1)) + pOffset;
      begin = 0;
      size = 0;
      int b = 0;
      for( B=bodies.begin(); B!=bodies.end(); ++B,++b ) {
        level = getLevel(Ibody[b]);
        cOffset = ((1 << 3*level) - 1)/7;
        pIndex = ((Ibody[b]-cOffset) >> 3*(level-l+1)) + pOffset;
        if( pIndex != icell ) {
          if( size <= NCRIT ) {
            for( int i=begin; i!=begin+size; ++i ) {
              Ibody[i] = icell;
            }
          }
          begin = b;
          size = 0;
          icell = pIndex;
        }
        size++;
      }
      if( size <= NCRIT ) {
        for( int i=begin; i!=begin+size; ++i ) {
          Ibody[i] = icell;
        }
      }
    }
  }

  void grow(Bodies &buffer, int level=0, int begin=0, int end=0) {
    int icell(Ibody[begin]),offset(begin),size(0);
    if( level == 0 ) level = getMaxLevel();
    if( end == 0 ) end = bodies.size();
    for( int b=begin; b!=end; ++b ) {
      if( Ibody[b] != icell ) {
        if( size > NCRIT ) {
          level++;
          setMorton(level,offset,offset+size);
          sort(Ibody,bodies,buffer,true,offset,offset+size);
          grow(buffer,level,offset,offset+size);
          level--;
        }
        offset = b;
        size = 0;
        icell = Ibody[b];
      }
      size++;
    }
    if( size > NCRIT ) {
      level++;
      setMorton(level,offset,offset+size);
      sort(Ibody,bodies,buffer,true,offset,offset+size);
      grow(buffer,level,offset,offset+size);
      level--;
    }
  }

  void link() {
    CCN = CC0 = CN = C0 = cells.begin();
    int icell(Ibody[0]),size(0),level(getLevel(icell));
    B_iter begin(bodies.begin());
    Cells buffer(cells.size());
    int b = 0;
    for( B=bodies.begin(); B!=bodies.end(); ++B,++b ) {
      if( Ibody[b] != icell ) {
        CN->NLEAF = size;
        CN->NCHILD = 0;
        CN->I = icell;
        CN->LEAF = begin;
        ++CN;
        if( getLevel(Ibody[b]) != level ) {
          linkParent(buffer);
          level = getLevel(Ibody[b]);
        }
        begin = B;
        size = 0;
        icell = Ibody[b];
      }
      size++;
    }
    CN->NLEAF = size;
    CN->NCHILD = 0;
    CN->I = icell;
    CN->LEAF = begin;
    ++CN;
    CC0 = CCN;
    for( int l=level; l>0; --l ) {
      linkParent(buffer);
    }

    for( C=C0; C!=CN; ++C ) {
      std::cout << C->I << std::endl;
    }
  }
};

#endif
