/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef bottomup_h
#define bottomup_h
#include "topdown.h"

//! Bottomup tree constructor
template<Equation equation>
class BottomUp : public TopDown<equation> {
public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using TreeStructure<equation>::buffer;                        //!< Buffer for MPI communication & sorting
  using TreeStructure<equation>::getLevel;                      //!< Get level from cell index

protected:
//! Max level for bottom up tree build
  int getMaxLevel(Bodies &bodies) {
    const long N = bodies.size() * MPISIZE;                     // Number of bodies
    int level;                                                  // Max level
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;   // Decide max level from N/Ncrit
    int MPIlevel = int(log(MPISIZE-1) / M_LN2 / 3) + 1;         // Level of local root cell
    if( MPISIZE == 1 ) MPIlevel = 0;                            // For serial execution local root cell is root cell
    if( MPIlevel > level ) {                                    // If process hierarchy is deeper than tree
//      std::cout << "Process hierarchy is deeper than tree @ rank" << MPIRANK << std::endl;
      level = MPIlevel;
    }
    return level;                                               // Return max level
  }

public:
//! Constructor
  BottomUp() : TopDown<equation>() {}
//! Destructor
  ~BottomUp() {}

//! Set cell index of all bodies
  void setIndex(Bodies &bodies, int level=-1, int begin=0, int end=0, bool update=false) {
    startTimer("Set index    ");                                // Start timer
    bigint i;                                                   // Levelwise cell index
    if( level == -1 ) level = getMaxLevel(bodies);              // Decide max level
    bigint off = ((1 << 3*level) - 1) / 7;                      // Offset for each level
    real r = R0 / (1 << (level-1));                             // Radius at finest level
    vec<3,int> nx;                                              // 3-D cell index
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        nx[d] = int( ( bodies[b].X[d] - (X0[d]-R0) ) / r );     //   3-D cell index
      }                                                         //  End loop over dimension
      i = 0;                                                    //  Initialize cell index
      for( int l=0; l!=level; ++l ) {                           //  Loop over levels of tree
        for( int d=0; d!=3; ++d ) {                             //   Loop over dimension
          i += nx[d] % 2 << (3 * l + d);                        //    Accumulate cell index
          nx[d] >>= 1;                                          //    Bitshift 3-D cell index
        }                                                       //   End loop over dimension
      }                                                         //  End loop over levels
      if( !update ) {                                           //  If this not an update
        bodies[b].ICELL = i+off;                                //   Store index in bodies
      } else if( i+off > bodies[b].ICELL ) {                    //  If the new cell index is larger
        bodies[b].ICELL = i+off;                                //   Store index in bodies
      }                                                         //  Endif for update
    }                                                           // End loop over bodies
    stopTimer("Set index    ",printNow);                        // Stop timer
  }

//! Prune tree by merging cells
  void prune(Bodies &bodies) {
    startTimer("Prune tree   ");                                // Start timer
    int maxLevel = getMaxLevel(bodies);                         // Max level for bottom up tree build
    for( int l=maxLevel; l>0; --l ) {                           // Loop upwards from bottom level
      int level = getLevel(bodies[0].ICELL);                    //  Current level
      bigint cOff = ((1 << 3 * level) - 1) / 7;                 //  Current ce;; index offset
      bigint pOff = ((1 << 3 * (l-1)) - 1) / 7;                 //  Parent cell index offset
      bigint index = ((bodies[0].ICELL-cOff) >> 3*(level-l+1)) + pOff;// Current cell index
      int begin = 0;                                            //  Begin cell index for bodies in cell
      int size = 0;                                             //  Number of bodies in cell
      int b = 0;                                                //  Current body index
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {//  Loop over bodies
        level = getLevel(B->ICELL);                             //   Level of twig
        cOff = ((1 << 3*level) - 1) / 7;                        //   Offset of twig
        bigint p = ((B->ICELL-cOff) >> 3*(level-l+1)) + pOff;   //   Cell index of parent cell
        if( p != index ) {                                      //   If it's a new parent cell
          if( size < NCRIT ) {                                  //    If parent cell has few enough bodies
            for( int i=begin; i!=begin+size; ++i ) {            //     Loop over bodies in that cell
              bodies[i].ICELL = index;                          //      Renumber cell index to parent cell
            }                                                   //     End loop over bodies in cell
          }                                                     //    Endif for merging
          begin = b;                                            //    Set new begin index
          size = 0;                                             //    Reset number of bodies
          index = p;                                            //    Set new parent cell
        }                                                       //   Endif for new cell
        size++;                                                 //   Increment body counter
      }                                                         //  End loop over bodies
      if( size < NCRIT ) {                                      //  If last parent cell has few enough bodies
        for( int i=begin; i!=begin+size; ++i ) {                //   Loop over bodies in that cell
          bodies[i].ICELL = index;                              //    Renumber cell index to parent cell
        }                                                       //   End loop over bodies in cell
      }                                                         //  Endif for merging
    }                                                           // End loop over levels
    stopTimer("Prune tree   ",printNow);                        // Stop timer
  }

//! Grow tree by splitting cells
  void grow(Bodies &bodies, int level=0, int begin=0, int end=0) {
    bigint index = bodies[begin].ICELL;                         // Initialize cell index
    int off=begin, size=0;                                      // Initialize offset, and size
    if( level == 0 ) level = getMaxLevel(bodies);               // Max level for bottom up tree build
    if( end == 0 ) end = bodies.size();                         // Default size is all bodies
    for( int b=begin; b!=end; ++b ) {                           // Loop over bodies under consideration
      if( bodies[b].ICELL != index ) {                          //  If it's a new cell
        if( size >= NCRIT ) {                                   //   If the cell has too many bodies
          level++;                                              //    Increment level
          setIndex(bodies,level,off,off+size);                  //    Set new cell index considering new level
          sortBodies(bodies,buffer,false,off,off+size);         //    Sort new cell index
          grow(bodies,level,off,off+size);                      //    Recursively grow tree
          level--;                                              //    Go back to previous level
        }                                                       //   Endif for splitting
        off = b;                                                //   Set new offset
        size = 0;                                               //   Reset number of bodies
        index = bodies[b].ICELL;                                //   Set new cell
      }                                                         //  Endif for new cell
      size++;                                                   //  Increment body counter
    }                                                           // End loop over bodies
    if( size >= NCRIT ) {                                       // If last cell has too many bodies
      level++;                                                  //  Increment level
      setIndex(bodies,level,off,off+size);                      //  Set new cell index considering new level
      sortBodies(bodies,buffer,false,off,off+size);             //  Sort new cell index
      grow(bodies,level,off,off+size);                          //  Recursively grow tree
      level--;                                                  //  Go back to previous level
    }                                                           // Endif for splitting
  }
};

#endif
