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

template<Equation equation>
class BottomUp : public TopDown<equation> {
public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Evaluator<equation>::MAXLEVEL;                          //!< Max depth of tree
  using TopDown<equation>::getMorton;                           //!< Get Morton index

private:
  inline int getMaxLevel(Bodies &bodies) const {
    const long N = bodies.size();
    int level;
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;
    return level;
  }

  inline void initCell(Cell &cell, int child, B_iter LEAF, int level, real diameter) const {
    cell.CHILD  = child;
    cell.LEAF   = LEAF;
    int ix = int((LEAF->X[0] + R0 - X0[0]) / diameter);
    int iy = int((LEAF->X[1] + R0 - X0[1]) / diameter);
    int iz = int((LEAF->X[2] + R0 - X0[2]) / diameter);
    cell.ICELL  = getMorton(ix,iy,iz,level);
    cell.X[0]   = diameter * (ix + .5) + X0[0] - R0;
    cell.X[1]   = diameter * (iy + .5) + X0[1] - R0;
    cell.X[2]   = diameter * (iz + .5) + X0[2] - R0;
    cell.R      = diameter * .5;
    cell.RCRIT  = cell.R / THETA;
  }

  inline void buildBottom(Bodies &bodies, Cells &cells) const {
    int I = -1;
    C_iter C;
    cells.clear();
    cells.reserve(1 << (3 * MAXLEVEL));
    float diameter = 2 * R0 / (1 << MAXLEVEL);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int IC = B->ICELL;
      if( IC != I ) {
        Cell cell;
        initCell(cell,0,B,MAXLEVEL,diameter);
        cells.push_back(cell);
        C = cells.end()-1;
        I = IC;
      }
      C->NCLEAF++;
      C->NDLEAF++;
    }
  }

  inline void twigs2cells(Cells &cells) const {
    int begin = 0, end = cells.size();
    float diameter = 2 * R0 / (1 << MAXLEVEL);
    for( int level=MAXLEVEL-1; level>=0; --level ) {
      int I = -1;
      int p = end - 1;
      int offset = ((1 << 3 * (level + 1)) - 1) / 7;
      diameter *= 2;
      for( int c=begin; c!=end; ++c ) {
        int IC = (cells[c].ICELL - offset) / 8;
        if( IC != I ) {
          Cell cell;
          initCell(cell,c,cells[c].LEAF,level,diameter);
          cells.push_back(cell);
          p++;
          I = IC;
        }
        cells[p].NCHILD++;
        cells[p].NDLEAF += cells[c].NDLEAF;
        cells[c].PARENT = p;
      }
      begin = end;
      end = cells.size();
    }
  }

protected:
  inline void setIndex(Bodies &bodies) const {
    float diameter = 2 * R0 / (1 << MAXLEVEL);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int ix = int((B->X[0] + R0 - X0[0]) / diameter);
      int iy = int((B->X[1] + R0 - X0[1]) / diameter);
      int iz = int((B->X[2] + R0 - X0[2]) / diameter);
      B->ICELL = getMorton(ix,iy,iz,MAXLEVEL);
    }
  }

  void setDomain(Bodies &bodies) {
    startTimer("Set domain");
    MAXLEVEL = getMaxLevel(bodies);
    vect xmin, xmax;
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];
      }
      X0 += B->X;
    }
    X0 /= bodies.size();
    for( int d=0; d!=3; ++d ) {
      X0[d] = int(X0[d]+.5);
      R0 = std::max(xmax[d] - X0[d], R0);
      R0 = std::max(X0[d] - xmin[d], R0);
    }
    R0 *= 1.000001;
    if( IMAGES != 0 ) {
      X0 = 0;
      R0 = M_PI;
    }
    stopTimer("Set domain",printNow);
  }

  void buildTree(Bodies &bodies, Cells &cells) {
    startTimer("Morton index");
    setIndex(bodies);
    Bodies buffer = bodies;
    stopTimer("Morton index",printNow);
    startTimer("Sort bodies");
    sortBodies(bodies,buffer,false);
    stopTimer("Sort bodies",printNow);
    startTimer("Build bottom");
    buildBottom(bodies,cells);
    stopTimer("Build bottom",printNow);
  }

  void linkTree(Cells &cells) {
    startTimer("Link tree");
    twigs2cells(cells);
    stopTimer("Link tree",printNow);
  }

};

#endif
