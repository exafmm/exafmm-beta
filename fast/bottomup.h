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

class BottomUp : public TopDown {
private:
  int getMaxLevel(Bodies &bodies) {
    const long N = bodies.size();
    int level;
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;
    return level;
  }

  inline void getIndex(Bodies &bodies) {
    float d = 2 * R0 / (1 << MAXLEVEL);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int ix = int((B->X[0] + R0 - X0[0]) / d);
      int iy = int((B->X[1] + R0 - X0[1]) / d);
      int iz = int((B->X[2] + R0 - X0[2]) / d);
      int id = 0;
      for( int l=0; l!=MAXLEVEL; ++l ) {
        id += ix % 2 << (3 * l);
        id += iy % 2 << (3 * l + 1);
        id += iz % 2 << (3 * l + 2);
        ix >>= 1;
        iy >>= 1;
        iz >>= 1;
      }
      B->ICELL = id;
    }
  }

  inline void initCell(Cell &cell, int child, B_iter LEAF, real diameter) {
    cell.NCHILD = 0;
    cell.NCLEAF = 0;
    cell.NDLEAF = 0;
    cell.CHILD  = child;
    cell.LEAF   = LEAF;
    int ix = int((LEAF->X[0] + R0 - X0[0]) / diameter);
    int iy = int((LEAF->X[1] + R0 - X0[1]) / diameter);
    int iz = int((LEAF->X[2] + R0 - X0[2]) / diameter);
    cell.X[0]   = diameter * (ix + .5) + X0[0] - R0;
    cell.X[1]   = diameter * (iy + .5) + X0[1] - R0;
    cell.X[2]   = diameter * (iz + .5) + X0[2] - R0;
    cell.R      = diameter * .5;
  }

  void buildBottom(Bodies &bodies, Cells &cells) {
    int I = -1;
    C_iter C;
    cells.clear();
    cells.reserve(1 << (3 * MAXLEVEL));
    float d = 2 * R0 / (1 << MAXLEVEL);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int IC = B->ICELL;
      if( IC != I ) {
        Cell cell;
        initCell(cell,0,B,d);
        cells.push_back(cell);
        C = cells.end()-1;
        I = IC;
      }
      C->NCLEAF++;
      C->NDLEAF++;
    }
  }

  void twigs2cells(Cells &cells) {
    int begin = 0, end = cells.size();
    float d = 2 * R0 / (1 << MAXLEVEL);
    for( int l=0; l!=MAXLEVEL; ++l ) {
      int div = (8 << (3 * l));
      int I = -1;
      int p = end - 1;
      d *= 2;
      for( int c=begin; c!=end; ++c ) {
        B_iter B = cells[c].LEAF;
        int IC = B->ICELL / div;
        if( IC != I ) {
          Cell cell;
          initCell(cell,c,cells[c].LEAF,d);
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
  void setDomain(Bodies &bodies) {
    startTimer("Set domain   ");
    MAXLEVEL = getMaxLevel(bodies);
    vect xmin, xmax;
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->X;                                               //  Sum positions
    }                                                           // End loop over bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 *= 1.000001;                                             // Add some leeway to root radius
    stopTimer("Set domain   ",printNow);
  }

  void buildTree(Bodies &bodies, Cells &cells) {
    startTimer("Morton index ");
    getIndex(bodies);
    Bodies buffer = bodies;
    stopTimer("Morton index ",printNow);
    startTimer("Sort bodies  ");
    sortBodies(bodies,buffer);
    stopTimer("Sort bodies  ",printNow);
    startTimer("Build bottom ");
    buildBottom(bodies,cells);
    stopTimer("Build bottom ",printNow);
  }

  void linkTree(Cells &cells) {
    startTimer("Link tree    ");
    twigs2cells(cells);
    stopTimer("Link tree    ",printNow);
  }

  void upwardPass(Cells &cells) {
    startTimer("Upward pass  ");
    setRootCell(cells);
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      real Rmax = 0;
      setCenter(C);
      P2M(C,Rmax);
      M2M(C,Rmax);
    }
#if Cartesian
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass  ",printNow);
  }

  void downwardPass(Cells &cells) const {
    for( C_iter C=cells.end()-2; C!=cells.begin()-1; --C ) {
      L2L(C);
      L2P(C);
    }
  }

};

#endif
