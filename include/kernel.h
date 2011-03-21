#ifndef kernel_h
#define kernel_h
#include "logger.h"
#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)

class Kernel : public Logger {                                  // Unified CPU/GPU kernel class
protected:
  B_iter BI0;                                                   // Target bodies begin iterator
  B_iter BIN;                                                   // Target bodies end iterator
  B_iter BJ0;                                                   // Source bodies begin iterator
  B_iter BJN;                                                   // Source bodies end iterator
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator
  vect   X0;                                                    // Center of root cell
  real   R0;                                                    // Radius of root cell
  vect   Xperiodic;                                             // Coordinate offset of periodic image

  std::vector<int>   keysHost;                                  // Offsets for rangeHost
  std::vector<int>   rangeHost;                                 // Offsets for sourceHost
  std::vector<float> constHost;                                 // Constants on host
  std::vector<float> sourceHost;                                // Sources on host
  std::vector<float> targetHost;                                // Targets on host
  Map                sourceBegin;                               // Define map for offset of source cells
  Map                sourceSize;                                // Define map for size of source cells
  Map                targetBegin;                               // Define map for offset of target cells

  int getLevel(bigint index) {                                  // Get level from cell index
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

public:
  Kernel() : X0(0), R0(0) {}                                    // Constructor
  ~Kernel() {}                                                  // Destructor

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain(Bodies &bodies) {                              // Set center and size of root cell
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    X0 = R0 = 0;                                                // Initialize center and size of root cell
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
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
    R0 = pow(2.,int(1. + log(R0) / M_LN2));                     // Add some leeway to root radius
  }

  void initialize();

  void precalculate();

  void allocGPU();

  void deallocGPU();

  void P2M();

  void M2M();

  void M2M_CPU();

  void M2L();

  void M2P();

  void P2P();

  void L2L();

  void L2P();

  void postcalculate();

  void finalize();
};

#endif
