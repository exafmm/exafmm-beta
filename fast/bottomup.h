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
      int ix = (B->X[0] + R0 - X0[0]) / d;
      int iy = (B->X[1] + R0 - X0[1]) / d;
      int iz = (B->X[2] + R0 - X0[2]) / d;
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
    int ix = (LEAF->X[0] + R0 - X0[0]) / diameter;
    int iy = (LEAF->X[1] + R0 - X0[1]) / diameter;
    int iz = (LEAF->X[2] + R0 - X0[2]) / diameter;
    cell.X[0]   = diameter * (ix + .5) + X0[0] - R0;
    cell.X[1]   = diameter * (iy + .5) + X0[1] - R0;
    cell.X[2]   = diameter * (iz + .5) + X0[2] - R0;
    cell.R      = diameter * .5;
  }

  void buildBottom(Bodies &bodies) {
    int I = -1;
    C_iter C;
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

  void twigs2cells() {
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
  }

  void build(Bodies &bodies) {
    startTimer("Morton index ");
    getIndex(bodies);
    Bodies buffer = bodies;
    stopTimer("Morton index ",printNow);
    startTimer("Sort bodies  ");
    sortBodies(bodies,buffer);
    stopTimer("Sort bodies  ",printNow);
    startTimer("Build bottom ");
    buildBottom(bodies);
    stopTimer("Build bottom ",printNow);
  }

  void link() {
    startTimer("Link tree    ");
    twigs2cells();
    NCELL = cells.size();
    C0 = cells.begin();
    CN = cells.end()-1;
    stopTimer("Link tree    ",printNow);
  }

  void upward() {
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      setCenter(C);
      P2M(C);
      M2M(C);
    }
#if CART
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M[1] /= C->M[0];
      C->M[2] /= C->M[0];
      C->M[3] /= C->M[0];
      C->M[4] /= C->M[0];
      C->M[5] /= C->M[0];
      C->M[6] /= C->M[0];
    }
#endif
    set_rcrit();
  }

public:
  BottomUp() {}
  ~BottomUp() {}
};

#endif
