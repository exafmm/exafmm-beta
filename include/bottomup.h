#ifndef bottomup_h
#define bottomup_h
#include "tree.h"

class BottomUpTreeConstructor : virtual public TreeStructure {
public:
  BottomUpTreeConstructor(Bodies &b) : TreeStructure(b){}       // Constructor
  ~BottomUpTreeConstructor() {}                                 // Destructor

  int getMaxLevel() {                                           // Max level for bottom up tree build
    int const N = bodies.size();                                // Number of bodies
    int level;                                                  // Define max level
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;   // Decide max level from N/Ncrit
    return level;                                               // Return max level
  }

  void setIndex(int level=0, int begin=0, int end=0 ) {         // Set cell index of all bodies
    bigint i;                                                   // Levelwise cell index
    if( level == 0 ) level = getMaxLevel();                     // Decide max level
    bigint off = ((1 << 3*level) - 1) / 7;                      // Offset for each level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // Define 3-D cell index
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over all bodies
      for( int d=0; d!=3; ++d )                                 //  Loop over dimension
        nx[d] = int( ( bodies.at(b).pos[d] - (X0[d]-R0) ) / r );//   3-D cell index
      i = 0;                                                    //  Initialize cell index
      for( int l=0; l!=level; ++l ) {                           //  Loop over all levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d] % 2 << (3 * l + d);                        //    Accumulate cell index
          nx[d] >>= 1;                                          //    Bitshift 3-D cell index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      Ibody[b] = i+off;                                         //  Store results with levelwise offset
    }                                                           // End loop over all bodies
  }

  void prune() {                                                // Prune tree by merging cells
    int maxLevel = getMaxLevel();                               // Max level for bottom up tree build
    for( int l=maxLevel; l>0; --l ) {                           // Loop upwards from bottom level
      int level = getLevel(Ibody[0]);                           //  Current level
      bigint cOff = ((1 << 3 * level) - 1) / 7;                 //  Current ce;; index offset
      bigint pOff = ((1 << 3 * (l-1)) - 1) / 7;                 //  Parent cell index offset
      bigint icell = ((Ibody[0]-cOff) >> 3*(level-l+1)) + pOff; //  Current cell index
      int begin = 0;                                            //  Begin cell index for bodies in cell
      int size = 0;                                             //  Number of bodies in cell
      int b = 0;                                                //  Current body index
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {//  Loop over all bodies
        level = getLevel(Ibody[b]);                             //   Level of twig
        cOff = ((1 << 3*level) - 1) / 7;                        //   Offset of twig
        bigint p = ((Ibody[b]-cOff) >> 3*(level-l+1)) + pOff;   //   Cell index of parent cell
        if( p != icell ) {                                      //   If it's a new parent cell
          if( size < NCRIT ) {                                  //    If parent cell has few enough bodies
            for( int i=begin; i!=begin+size; ++i ) {            //     Loop over all bodies in that cell
              Ibody[i] = icell;                                 //      Renumber cell index to parent cell
            }                                                   //     End loop over all bodies in cell
          }                                                     //    Endif for merging
          begin = b;                                            //    Set new begin index
          size = 0;                                             //    Reset number of bodies
          icell = p;                                            //    Set new parent cell
        }                                                       //   Endif for new cell
        size++;                                                 //   Increment body counter
      }                                                         //  End loop over all bodies
      if( size < NCRIT ) {                                      //  If last parent cell has few enough bodies
        for( int i=begin; i!=begin+size; ++i ) {                //   Loop over all bodies in that cell
          Ibody[i] = icell;                                     //    Renumber cell index to parent cell
        }                                                       //   End loop over all bodies in cell
      }                                                         //  Endif for merging
    }                                                           // End loop over levels
  }

  void grow(Bodies &buffer, int level=0, int begin=0, int end=0) { // Grow tree by splitting cells
    bigint icell(Ibody[begin]);                                 // Initialize cell index
    int off(begin), size(0);                                    // Initialize offset, and size
    if( level == 0 ) level = getMaxLevel();                     // Max level for bottom up tree build
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over all bodies under consideration
      if( Ibody[b] != icell ) {                                 //  If it's a new cell
        if( size >= NCRIT ) {                                   //   If the cell has too many bodies
          level++;                                              //    Increment level
          setIndex(level,off,off+size);                         //    Set new cell index considering new level
          sort(Ibody,bodies,buffer,true,off,off+size);          //    Sort new cell index
          grow(buffer,level,off,off+size);                      //    Recursively grow tree
          level--;                                              //    Go back to previous level
        }                                                       //   Endif for splitting
        off = b;                                                //   Set new offset
        size = 0;                                               //   Reset number of bodies
        icell = Ibody[b];                                       //   Set new cell
      }                                                         //  Endif for new cell
      size++;                                                   //  Increment body counter
    }                                                           // End loop over bodies
    if( size >= NCRIT ) {                                       // If last cell has too many bodies
      level++;                                                  //  Increment level
      setIndex(level,off,off+size);                             //  Set new cell index considering new level
      sort(Ibody,bodies,buffer,true,off,off+size);              //  Sort new cell index
      grow(buffer,level,off,off+size);                          //  Recursively grow tree
      level--;                                                  //  Go back to previous level
    }                                                           // Endif for splitting
  }
};

#endif
