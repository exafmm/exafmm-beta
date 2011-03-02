#ifndef kernel_h
#define kernel_h
#define KERNEL
#include "types.h"
#undef KERNEL

class Kernel {
protected:
  B_iter BI0;                                                   // Target bodies begin iterator
  B_iter BIN;                                                   // Target bodies end iterator
  B_iter BJ0;                                                   // Source bodies begin iterator
  B_iter BJN;                                                   // Source bodies end iterator
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator

  std::vector<int> keysHost;                                    // Keys on host
  std::vector<int> rangeHost;                                   // Ranges on host
  std::vector<real> sourceHost;                                 // Sources on host
  std::vector<real> targetHost;                                 // Targets on host

  int getLevel(bigint index) {                                  // Get level from cell index
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

public:
  void initialize();

  void P2M();

  void M2M();

  void M2L();

  void M2P();

  void P2P();

  void L2L();

  void L2P();

  void finalize();
};

#endif
