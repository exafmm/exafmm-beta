#ifndef tree_h
#define tree_h
#include "body.h"

class TreeStructure {
protected:
  bodies &B;                                                    // Bodies in the tree
  vect X0;                                                      // Center of root cell
  real R0;                                                      // Radius of root cell
  Cell *C0,*CN,*C;                                              // Pointer to cells
public:
  TreeStructure(bodies &b) : B(b),X0(0),R0(0) {}                // Constructor
  ~TreeStructure() {}                                           // Destructor

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain() {
    vect xmin,xmax;                                             // Min,Max of domain
    B.begin();                                                  // Reset bodies counter
    xmin = xmax = B.pos();                                      // Initialize xmin,xmax
    for( B=B.begin(); B!=B.end(); ++B ) {                       // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B.pos()[d] < xmin[d]) xmin[d] = B.pos()[d];     //   Determine xmin
        else if(B.pos()[d] > xmax[d]) xmax[d] = B.pos()[d];     //   Determine xmax
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

  void allocateCells() {
    C0 = new Cell [B.size()];
  }

  void deallocateCells() {
    delete[] C0;
  }
};

#endif
