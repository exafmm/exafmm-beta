#ifndef tree_h
#define tree_h
#include "types.h"

class TreeStructure {
protected:
  Bodies &bodies;                                               // Bodies in the tree
  B_iter B,Bi,Bj;                                               // Body iterators
  Cells  cells;                                                 // Cells in the tree
  C_iter C,C0,CN,CC0,CCN;                                       // Cell iterators
  vect X0;                                                      // Center of root cell
  real R0;                                                      // Radius of root cell
public:
  bigint *Ibody;                                                // Morton index of body
  bigint *Icell;                                                // Morton index of cell

  TreeStructure(Bodies &b) : bodies(b),X0(0),R0(0) {            // Constructor
    int const N = bodies.size();
    Ibody = new bigint [N];
    Icell = new bigint [N];
    cells.resize(N);
    for( C=cells.begin(); C!=cells.end(); ++C ) {
      C->NLEAF = 0;
      C->NCHILD = 0;
    }
  }
  ~TreeStructure() {                                            // Destructor
    delete[] Ibody;
    delete[] Icell;
  }

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain() {
    vect xmin,xmax;                                             // Min,Max of domain
    B = bodies.begin();                                         // Reset body iterator
    xmin = xmax = B->pos;                                       // Initialize xmin,xmax
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < xmin[d]) xmin[d] = B->pos[d];       //   Determine xmin
        else if(B->pos[d] > xmax[d]) xmax[d] = B->pos[d];       //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->pos;                                             //  Sum positions
    }                                                           // End loop over all bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d]-X0[d],R0);                          //  Calculate max distance from center
      R0 = std::max(X0[d]-xmin[d],R0);                          //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1.+log(R0)/M_LN2));                         // Add some leeway to root radius
  }
};

#endif
